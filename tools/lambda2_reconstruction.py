#!/usr/bin/env python3
"""Reconstruct lambda_2 from gamma_1 and propagate the MKNO 1-RDM EOM.

Target: a 1-RDM EOM that preserves 0 <= gamma_1 <= 1 and differs meaningfully from TDHF.

Derivation: qdiss.eom.gwt cumulant form of the rank-1 EOM, unapproximated -- 2 one-body and
12 two-body terms over gamma_1 and lambda_2, of which 4 carry lambda_2. Validated against
qdiss.eom.fullci at 1e-15. Propagated state: gamma_1 alone, with lambda_2 := R[gamma_1] a
reconstruction. Ground truth: exact wavefunction propagation.

The reconstruction must satisfy the partial-trace identity

    sum_x lambda_2[px,qx] = -(gamma_1 - gamma_1^2) = -A ,   A = gamma_1 eta_1

which is affine, so there is a family. Two privileged members are compared.

  minnorm   the minimum Frobenius norm solution. By orbital covariance it lies in the span of
            A wedge I and I wedge I, giving
                lambda_2 = -(1/(M-2)) (A^I + I^A) + [tr A / ((M-1)(M-2))] (I^I)
            Canonical and unique, but its coefficient decays as 1/(M-2), so it vanishes in a
            large basis and collapses to TDHF.

  blocked   demand that EVERY index carry the blocking factor rather than only the trace.
            Take lambda_2 = A^f + f^A with f = g(gamma_1), so A and f commute and
                tr(f) A + tr(A) f - 2 A f = -A
            In the natural-orbital basis with a_i = n_i (1 - n_i), T = sum_j a_j,
                f_i = -a_i / [ (1+Q) (T - 2 a_i) ] ,   Q = sum_j a_j / (T - 2 a_j)
            Both factors are proportional to a_i, so each vanishes as an orbital empties or
            fills. This is derived from the constraint, not inserted.

Baselines: TDHF (lambda_2 = 0) and exact.
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


def recon(g1, M, kind):
    if kind == "tdhf":
        return np.zeros((M,) * 4, dtype=complex)
    A = g1 - g1 @ g1
    A = 0.5 * (A + A.conj().T)
    I = np.eye(M, dtype=complex)
    if kind == "minnorm":
        return (-1.0 / (M - 2)) * sym_wedge(A, I) + (np.trace(A) / ((M - 1) * (M - 2))) * wedge(I, I) * 2
    a, U = np.linalg.eigh(A)
    T = a.sum()
    den = T - 2 * a
    if np.abs(den).min() < 1e-10:
        raise RuntimeError("degenerate denominator in blocked reconstruction")
    Q = (a / den).sum()
    f = -a / ((1.0 + Q) * den)
    F = U @ np.diag(f.astype(complex)) @ U.conj().T
    return sym_wedge(A, F)


def rk4(f, y, dt):
    k1 = f(y); k2 = f(y + dt / 2 * k1); k3 = f(y + dt / 2 * k2); k4 = f(y + dt * k3)
    return y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def run(M, N, seed=7, dt=0.002, nsteps=500, every=100):
    EYE = np.eye(M, dtype=complex)
    rng = np.random.default_rng(seed)
    h, V = random_h(M, rng), random_V(M, rng)
    psi0 = fullci.random_state(M, N, rng)
    H = fullci.build_H(h, V, M)
    w, U = np.linalg.eigh(H)
    g1_0 = fullci.rdm(psi0, M, 1)
    A0 = g1_0 - g1_0 @ g1_0

    print("=== M=%d N=%d ===" % (M, N))
    for kind in ("minnorm", "blocked"):
        l2 = recon(g1_0, M, kind)
        err = np.abs(np.einsum("pxqx->pq", l2) + A0).max()
        print("   %-8s trace-identity residual %.3e   |lambda_2| %.5f  (exact |lambda_2| %.5f)"
              % (kind, err, np.abs(l2).max(),
                 np.abs(fullci.rdm(psi0, M, 2) - wedge(g1_0, g1_0)).max()))

    def d1(g1, kind):
        t = {("D", 1): EYE, ("G", 1): g1, ("L", 2): recon(g1, M, kind)}
        return assemble(CUM1[0], CUM1[1], h, V, t, "pq", M)

    st = {k: g1_0.copy() for k in ("tdhf", "minnorm", "blocked")}
    occ0 = np.sort(np.linalg.eigvalsh(0.5 * (g1_0 + g1_0.conj().T)))
    print()
    print("   %6s | %-24s | %-24s | %-24s | %8s"
          % ("t", "TDHF", "minnorm", "blocked", "exact"))
    print("   %6s | %7s %7s %8s | %7s %7s %8s | %7s %7s %8s | %8s"
          % ("", "nmin", "drift", "err", "nmin", "drift", "err", "nmin", "drift", "err", "nmin"))
    for step in range(nsteps + 1):
        if step % every == 0:
            t = step * dt
            psit = U @ (np.exp(-1j * w * t) * (U.conj().T @ psi0))
            e1 = fullci.rdm(psit, M, 1)
            cells = []
            for k in ("tdhf", "minnorm", "blocked"):
                g = st[k]
                o = np.sort(np.linalg.eigvalsh(0.5 * (g + g.conj().T)))
                cells.append((o.min(), np.abs(o - occ0).max(), np.abs(g - e1).max()))
            print("   %6.2f | %7.4f %7.4f %8.4f | %7.4f %7.4f %8.4f | %7.4f %7.4f %8.4f | %8.4f"
                  % (t, cells[0][0], cells[0][1], cells[0][2],
                     cells[1][0], cells[1][1], cells[1][2],
                     cells[2][0], cells[2][1], cells[2][2],
                     np.linalg.eigvalsh(0.5 * (e1 + e1.conj().T)).min()))
        if step < nsteps:
            for k in ("tdhf", "minnorm", "blocked"):
                st[k] = rk4(lambda y, kk=k: d1(y, kk), st[k], dt)
    print()


run(6, 3)
run(8, 4)
print("drift is the largest change in the sorted occupation spectrum since t=0; TDHF is")
print("isospectral so its drift is integrator noise. A reconstruction that differs")
print("meaningfully from TDHF must show nonzero drift while keeping nmin >= 0.")
