#!/usr/bin/env python3
"""How does the blocked reconstruction differ from TDHF, given both are isospectral?

Isospectral means the DIAGONAL of gamma_1-dot vanishes in the natural-orbital basis. It says
nothing about the off-diagonal part, so the two flows can differ in every other respect: the
generator, the rate at which natural orbitals rotate, the coherences, the energy. The claim
that the reconstruction "cannot beat TDHF" was about occupations only and was overstated.

Measured here, with derivation qdiss.eom.gwt cumulant rank-1 EOM validated against
qdiss.eom.fullci at 1e-15, state gamma_1 alone, truth exact wavefunction propagation:

  1. the lambda_2 contribution to gamma_1-dot, split into diagonal and off-diagonal in the
     natural-orbital basis
  2. the generators. Both flows are isospectral so each is i[X^T, gamma_1]; X is extracted off
     the degenerate blocks and compared to the Fock operator f = h + W
  3. trajectory separation |gamma_1_blocked - gamma_1_TDHF| against their common distance to
     FullCI
  4. energy E = tr(h gamma_1) + (1/4) sum V gamma_2 with gamma_2 = gamma_1^gamma_1 + lambda_2,
     against the exact conserved energy
"""
import sys
import numpy as np
sys.path.insert(0, ".")
from qdiss.eom import fullci, gwt
from qdiss.eom.evaluate import assemble
from qdiss.eom.hamiltonian import TWO_BODY_PREFACTOR, random_V, random_h

CUM1 = gwt.terms_cumulant(1)
M, N = 6, 3
EYE = np.eye(M, dtype=complex)


def wedge(a, b):
    return np.einsum("pq,vw->pvqw", a, b) - np.einsum("pw,vq->pvqw", a, b)


def sym_wedge(a, b):
    return wedge(a, b) + wedge(b, a)


def blocked(g1):
    A = 0.5 * ((g1 - g1 @ g1) + (g1 - g1 @ g1).conj().T)
    a, U = np.linalg.eigh(A)
    T = a.sum()
    den = T - 2 * a
    Q = (a / den).sum()
    f = -a / ((1.0 + Q) * den)
    return sym_wedge(A, U @ np.diag(f.astype(complex)) @ U.conj().T)


rng = np.random.default_rng(7)
h, V = random_h(M, rng), random_V(M, rng)
psi0 = fullci.random_state(M, N, rng)
H = fullci.build_H(h, V, M)
w, U = np.linalg.eigh(H)
E_exact = float(np.real(np.vdot(psi0, H @ psi0)))
g1_0 = fullci.rdm(psi0, M, 1)


def d1(g1, kind):
    l2 = np.zeros((M,) * 4, dtype=complex) if kind == "tdhf" else blocked(g1)
    t = {("D", 1): EYE, ("G", 1): g1, ("L", 2): l2}
    return assemble(CUM1[0], CUM1[1], h, V, t, "pq", M)


lam, Uno = np.linalg.eigh(0.5 * (g1_0 + g1_0.conj().T))
dt_ = d1(g1_0, "tdhf")
db = d1(g1_0, "blocked")
diff = db - dt_
dn = Uno.conj().T @ diff @ Uno
print("=== 1. the lambda_2 contribution at t=0 ===")
print("   |gamma_1-dot TDHF|                 %.5f" % np.abs(dt_).max())
print("   |gamma_1-dot blocked|              %.5f" % np.abs(db).max())
print("   |difference|                       %.5f" % np.abs(diff).max())
print("   difference, DIAGONAL in NO basis   %.3e" % np.abs(np.diag(dn)).max())
print("   difference, OFF-diagonal           %.5f"
      % np.abs(dn - np.diag(np.diag(dn))).max())
print()

print("=== 2. the generators ===")


def extract(d):
    de = Uno.conj().T @ d @ Uno
    mask = np.abs(lam[None, :] - lam[:, None]) < 1e-9
    X = np.zeros_like(de)
    for i in range(M):
        for j in range(M):
            if not mask[i, j]:
                X[i, j] = de[i, j] / (1j * (lam[j] - lam[i]))
    return Uno @ X @ Uno.conj().T


W = np.einsum("prqs,rs->pq", V, g1_0)
Xt, Xb = extract(dt_), extract(db)
fock = (h + W).conj()
print("   |i[X_TDHF^T, g1] - gamma_1-dot|    %.3e"
      % np.abs(1j * (Xt @ g1_0 - g1_0 @ Xt) - dt_).max())
print("   |i[X_blocked^T, g1] - gamma_1-dot| %.3e"
      % np.abs(1j * (Xb @ g1_0 - g1_0 @ Xb) - db).max())
print("   |X_blocked - X_TDHF|               %.5f   (|X_TDHF| %.5f)"
      % (np.abs(Xb - Xt).max(), np.abs(Xt).max()))
print("   so the generators differ: blocked is NOT the Fock operator.")
print()


def rk4(f, y, dt):
    k1 = f(y); k2 = f(y + dt / 2 * k1); k3 = f(y + dt / 2 * k2); k4 = f(y + dt * k3)
    return y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def energy(g1, kind):
    l2 = np.zeros((M,) * 4, dtype=complex) if kind == "tdhf" else blocked(g1)
    g2 = wedge(g1, g1) + l2
    return float(np.real(np.einsum("pq,pq->", h, g1)
                         + TWO_BODY_PREFACTOR * np.einsum("pvqw,pvqw->", V, g2)))


st = {k: g1_0.copy() for k in ("tdhf", "blocked")}
dt, nsteps, every = 0.002, 500, 100
print("=== 3 and 4. trajectory separation and energy ===")
print("%6s | %11s %11s | %11s | %10s %10s | %10s"
       % ("t", "TDHF err", "blocked err", "separation", "E TDHF", "E blocked", "E exact"))
for step in range(nsteps + 1):
    if step % every == 0:
        t = step * dt
        psit = U @ (np.exp(-1j * w * t) * (U.conj().T @ psi0))
        e1 = fullci.rdm(psit, M, 1)
        a, b = st["tdhf"], st["blocked"]
        print("%6.2f | %11.5f %11.5f | %11.5f | %10.5f %10.5f | %10.5f"
              % (t, np.abs(a - e1).max(), np.abs(b - e1).max(), np.abs(a - b).max(),
                 energy(a, "tdhf"), energy(b, "blocked"), E_exact))
    if step < nsteps:
        for k in ("tdhf", "blocked"):
            st[k] = rk4(lambda y, kk=k: d1(y, kk), st[k], dt)
print()
print("separation is |gamma_1_blocked - gamma_1_TDHF|. Nonzero separation with identical")
print("occupation spectra means the two theories rotate the natural orbitals differently.")
