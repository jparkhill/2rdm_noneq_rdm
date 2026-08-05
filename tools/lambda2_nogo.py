#!/usr/bin/env python3
"""Why a reconstruction lambda_2[gamma_1] cannot differ from TDHF in the occupations.

Both candidates satisfying the partial-trace identity propagated with drift 0.0000: the
occupation spectrum never moved. The suspected mechanism is structural rather than accidental.

If lambda_2 is a matrix FUNCTION of gamma_1 -- built from A = gamma_1 - gamma_1^2 and any
g(gamma_1) -- then every factor commutes with gamma_1 and is diagonal in the natural-orbital
basis. A wedge of two such matrices is nonzero only when the index pairs coincide as sets, so
lambda_2 has no coherences in the natural-orbital pair basis. Its contribution to the DIAGONAL
of gamma_1-dot then involves only diagonal elements of V, which are real, and the overall
factor of i makes the diagonal contribution imaginary; since gamma_1-dot must be Hermitian,
that contribution has to vanish. So ndot = 0 and the flow is isospectral.

Tested three ways: measure the diagonal of gamma_1-dot in the natural-orbital basis directly;
confirm the corrected minimum-norm form satisfies the identity; and add a V-DEPENDENT
trace-free piece, which should restore occupation motion if the diagnosis is right.

Trace-free projection. For any antisymmetric W with partial trace S = sum_x W[px,qx],
    sum_x (a^b)[px,qx] = a tr(b) - a b
    sum_x (S^I + I^S)  = (M-2) S + tr(S) I ,   sum_x (I^I) = (M-1) I
so W - (1/(M-2))(S^I + I^S) + [tr S/((M-2)(M-1))](I^I) has zero partial trace.

Derivation: qdiss.eom.gwt cumulant rank-1 EOM, validated against qdiss.eom.fullci at 1e-15.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
from qdiss.eom import fullci, gwt
from qdiss.eom.evaluate import assemble
from qdiss.eom.hamiltonian import random_V, random_h

CUM1 = gwt.terms_cumulant(1)


def wedge(a, b):
    return np.einsum("pq,vw->pvqw", a, b) - np.einsum("pw,vq->pvqw", a, b)


def sym_wedge(a, b):
    return wedge(a, b) + wedge(b, a)


def trace_free(W, M):
    S = np.einsum("pxqx->pq", W)
    I = np.eye(M, dtype=complex)
    return (W - (1.0 / (M - 2)) * sym_wedge(S, I)
            + (np.trace(S) / ((M - 2) * (M - 1))) * wedge(I, I))


def blocked(g1, M):
    A = 0.5 * ((g1 - g1 @ g1) + (g1 - g1 @ g1).conj().T)
    a, U = np.linalg.eigh(A)
    T = a.sum()
    den = T - 2 * a
    Q = (a / den).sum()
    f = -a / ((1.0 + Q) * den)
    F = U @ np.diag(f.astype(complex)) @ U.conj().T
    return sym_wedge(A, F)


def minnorm(g1, M):
    A = 0.5 * ((g1 - g1 @ g1) + (g1 - g1 @ g1).conj().T)
    I = np.eye(M, dtype=complex)
    return (-(1.0 / (M - 2)) * sym_wedge(A, I)
            + (np.trace(A) / ((M - 2) * (M - 1))) * wedge(I, I))


M, N = 6, 3
EYE = np.eye(M, dtype=complex)
rng = np.random.default_rng(7)
h, V = random_h(M, rng), random_V(M, rng)
psi0 = fullci.random_state(M, N, rng)
H = fullci.build_H(h, V, M)
w, U = np.linalg.eigh(H)
g1_0 = fullci.rdm(psi0, M, 1)
A0 = g1_0 - g1_0 @ g1_0
Vpert = trace_free(V - np.transpose(V, (1, 0, 2, 3)), M)

print("trace-identity residuals at t=0")
for nm, fn in (("minnorm (corrected)", minnorm), ("blocked", blocked)):
    print("   %-20s %.3e" % (nm, np.abs(np.einsum("pxqx->pq", fn(g1_0, M)) + A0).max()))
print("   %-20s %.3e (must be zero)"
      % ("V-perturbation", np.abs(np.einsum("pxqx->pq", Vpert)).max()))
print()


def recon(g1, kind, c=0.0):
    if kind == "tdhf":
        return np.zeros((M,) * 4, dtype=complex)
    base = minnorm(g1, M) if kind == "minnorm" else blocked(g1, M)
    return base + c * Vpert


def d1(g1, kind, c=0.0):
    t = {("D", 1): EYE, ("G", 1): g1, ("L", 2): recon(g1, kind, c)}
    return assemble(CUM1[0], CUM1[1], h, V, t, "pq", M)


lam0, Uno = np.linalg.eigh(0.5 * (g1_0 + g1_0.conj().T))
print("diagonal of gamma_1-dot in the natural-orbital basis at t=0")
for nm, kind, c in (("TDHF", "tdhf", 0.0), ("minnorm", "minnorm", 0.0),
                    ("blocked", "blocked", 0.0), ("blocked + 0.3 V", "blocked", 0.3),
                    ("exact lambda_2", "exact", 0.0)):
    if kind == "exact":
        t = {("D", 1): EYE, ("G", 1): g1_0,
             ("L", 2): fullci.rdm(psi0, M, 2) - wedge(g1_0, g1_0)}
        d = assemble(CUM1[0], CUM1[1], h, V, t, "pq", M)
    else:
        d = d1(g1_0, kind, c)
    print("   %-16s max|ndot| = %.3e" % (nm, np.abs(np.real(np.diag(Uno.conj().T @ d @ Uno))).max()))
print()


def rk4(f, y, dt):
    k1 = f(y); k2 = f(y + dt / 2 * k1); k3 = f(y + dt / 2 * k2); k4 = f(y + dt * k3)
    return y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


runs = [("blocked", "blocked", 0.0), ("blocked+0.1V", "blocked", 0.1),
        ("blocked+0.3V", "blocked", 0.3)]
st = {nm: g1_0.copy() for nm, _, _ in runs}
occ0 = np.sort(lam0)
dt, nsteps, every = 0.002, 500, 125
print("%6s | %s | %8s" % ("t", " | ".join("%-22s" % nm for nm, _, _ in runs), "exact"))
print("%6s | %s | %8s"
      % ("", " | ".join("%7s %7s %6s" % ("nmin", "drift", "err") for _ in runs), "nmin"))
for step in range(nsteps + 1):
    if step % every == 0:
        t = step * dt
        psit = U @ (np.exp(-1j * w * t) * (U.conj().T @ psi0))
        e1 = fullci.rdm(psit, M, 1)
        cells = []
        for nm, kind, c in runs:
            g = st[nm]
            o = np.sort(np.linalg.eigvalsh(0.5 * (g + g.conj().T)))
            cells.append("%7.4f %7.4f %6.3f"
                         % (o.min(), np.abs(o - occ0).max(), np.abs(g - e1).max()))
        print("%6.2f | %s | %8.4f"
              % (t, " | ".join(cells),
                 np.linalg.eigvalsh(0.5 * (e1 + e1.conj().T)).min()))
    if step < nsteps:
        for nm, kind, c in runs:
            st[nm] = rk4(lambda y, kk=kind, cc=c: d1(y, kk, cc), st[nm], dt)
