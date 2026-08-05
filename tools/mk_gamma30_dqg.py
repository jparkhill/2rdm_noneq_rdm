#!/usr/bin/env python3
"""DQG conditions along MK_gamma30 propagation.

Derivation: qdiss.eom.gwt, MKNO cumulant form (12 + 140 terms over gamma_1, lambda_2,
lambda_3), validated against qdiss.eom.fullci at 1e-15. Propagated state (gamma_1, lambda_2)
with gamma_1-dot exact and lambda_3 := -(A + B), so the implied gamma_3 is exactly zero.
Ground truth is exact wavefunction propagation.

This is NOT the same dynamics as the gamma_3 = 0 closure propagating gamma_2 alone, which
was isospectral and preserved D. Here gamma_1 evolves by its own exact equation and gamma_2
is reconstructed as gamma_1 wedge gamma_1 + lambda_2, so nothing guarantees isospectrality
and D may well move.

Monitored at each output time:
    D, Q, G   smallest eigenvalue, built complex-safely here rather than through
              dqg.antisym_matrix, which allocates a real array and would discard the
              imaginary part of a complex 2-RDM
    Pauli     smallest natural occupation of the propagated gamma_1
    consistency  |contract(gamma_2) - gamma_1|, which the exact theory keeps at zero and
              which nothing in this closure enforces
"""
import sys
import numpy as np
sys.path.insert(0, ".")
from qdiss.eom import fullci, gwt
from qdiss.eom.evaluate import assemble
from qdiss.eom.hamiltonian import random_V, random_h
from qdiss.representability import dqg
from qdiss.symbolic import mkno

UP3, LO3 = ("p", "v", "x"), ("q", "w", "y")
ORB = {"A": [], "B": []}
for s, f in mkno.gamma_to_cumulants(UP3, LO3):
    if any(n == "L" and r == 3 for n, r, u, l in f):
        continue
    k = "B" if any(n == "L" and r == 2 for n, r, u, l in f) else "A"
    ORB[k].append((s, ",".join("".join(u) + "".join(l) for n, r, u, l in f) + "->pvxqwy",
                   [(n, r) for n, r, u, l in f]))
CUM1, CUM2 = gwt.terms_cumulant(1), gwt.terms_cumulant(2)

print("=== probe: does dqg.antisym_matrix survive a complex tensor? ===")
Tc = np.zeros((4, 4, 4, 4), dtype=complex)
Tc[0, 1, 0, 1] = 1.0 + 2.0j
Tc[1, 0, 1, 0] = 1.0 + 2.0j
Tc[0, 1, 1, 0] = -(1.0 + 2.0j)
Tc[1, 0, 0, 1] = -(1.0 + 2.0j)
try:
    Mx, _ = dqg.antisym_matrix(Tc)
    print("   returned dtype %s, element [0,0] = %s   (input had imaginary part 2.0)"
          % (Mx.dtype, Mx[0, 0]))
    print("   VERDICT: imaginary part is %s"
          % ("DISCARDED" if not np.iscomplexobj(Mx) else "preserved"))
except Exception as exc:
    print("   raised %s: %s" % (type(exc).__name__, exc))
print()


def wedge(a, b):
    return np.einsum("pq,vw->pvqw", a, b) - np.einsum("pw,vq->pvqw", a, b)


def orbit(key, g1, l2, M):
    t = {("G", 1): g1, ("L", 2): l2}
    out = np.zeros((M,) * 6, dtype=complex)
    for s, sub, keys in ORB[key]:
        out = out + s * np.einsum(sub, *[t[k] for k in keys], optimize=True)
    return out


def pair_matrix(T, M):
    P = [(p, v) for p in range(M) for v in range(p + 1, M)]
    A = np.array([[2.0 * T[p, v, q, w] for (q, w) in P] for (p, v) in P])
    return 0.5 * (A + A.conj().T)


def mineig_dqg(g2, g1, M):
    D = pair_matrix(g2, M)
    Q = pair_matrix(dqg.build_Q(g2, g1), M)
    Gt = dqg.build_G(g2, g1).reshape(M * M, M * M)
    Gt = 0.5 * (Gt + Gt.conj().T)
    return (np.linalg.eigvalsh(D).min(), np.linalg.eigvalsh(Q).min(),
            np.linalg.eigvalsh(Gt).min())


def rk4(f, y, dt):
    def add(a, b, c):
        return tuple(x + c * k for x, k in zip(a, b))
    k1 = f(y); k2 = f(add(y, k1, dt / 2)); k3 = f(add(y, k2, dt / 2)); k4 = f(add(y, k3, dt))
    return tuple(x + dt / 6 * (a + 2 * b + 2 * c + d)
                 for x, a, b, c, d in zip(y, k1, k2, k3, k4))


def run(M, N, dt=0.002, nsteps=500, every=100, seed=7):
    EYE = np.eye(M, dtype=complex)
    rng = np.random.default_rng(seed)
    h, V = random_h(M, rng), random_V(M, rng)
    psi0 = fullci.random_state(M, N, rng)
    H = fullci.build_H(h, V, M)
    w, U = np.linalg.eigh(H)

    def rhs(state):
        g1, l2 = state
        l3 = -(orbit("A", g1, l2, M) + orbit("B", g1, l2, M))
        t = {("D", 1): EYE, ("G", 1): g1, ("L", 2): l2, ("L", 3): l3}
        d1 = assemble(CUM1[0], CUM1[1], h, V, t, "pq", M)
        d2 = assemble(CUM2[0], CUM2[1], h, V, t, "pvqw", M)
        return (d1, d2 - (wedge(d1, g1) + wedge(g1, d1)))

    g1_0 = fullci.rdm(psi0, M, 1)
    l2_0 = fullci.rdm(psi0, M, 2) - wedge(g1_0, g1_0)
    st = (g1_0.copy(), l2_0.copy())
    print("=== MK_gamma30, M=%d N=%d  (gamma_3 = 0 is %s here) ==="
          % (M, N, "EXACT" if N == 2 else "FALSE"))
    print("%6s | %9s %9s %9s | %9s %9s | %10s"
          % ("t", "min D", "min Q", "min G", "n_min", "n_min FCI", "consist"))
    for step in range(nsteps + 1):
        if step % every == 0:
            t = step * dt
            g1, l2 = st
            g2 = wedge(g1, g1) + l2
            dd, qq, gg = mineig_dqg(g2, g1, M)
            psit = U @ (np.exp(-1j * w * t) * (U.conj().T @ psi0))
            e1 = fullci.rdm(psit, M, 1)
            cons = np.abs(np.einsum("pxqx->pq", g2) / (N - 1) - g1).max()
            print("%6.2f | %9.5f %9.5f %9.5f | %9.5f %9.5f | %10.2e"
                  % (t, dd, qq, gg,
                     np.linalg.eigvalsh(0.5 * (g1 + g1.conj().T)).min(),
                     np.linalg.eigvalsh(0.5 * (e1 + e1.conj().T)).min(), cons))
        if step < nsteps:
            st = rk4(rhs, st, dt)
    print()


run(4, 2)
run(6, 3)
print("Exact dynamics keeps min D, min Q, min G >= 0 and consist = 0 at all times.")
