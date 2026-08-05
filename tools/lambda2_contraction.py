#!/usr/bin/env python3
"""The lambda_2 contraction identity is the particle-hole blocking factor.

The MKNO 1-RDM EOM is a function of delta, gamma_1 and lambda_2 only -- factors
D1 x14, G1 x18, L2 x4 -- with eta_1 = delta - gamma_1 a rewriting rather than a separate
object. So in MK_gamma30 both arguments of gamma_1-dot are propagated state variables and
there is no external gamma_2 input. gamma_1-dot is exact in its arguments; the trajectory
degrades only because lambda_2 is wrong, and the whole approximation sits in lambda_2-dot.

Contracting gamma_2 = gamma_1 wedge gamma_1 + lambda_2 gives an identity every physical
lambda_2 must satisfy:

    sum_x gamma_2[px,qx]            = (N-1) gamma_1
    sum_x (gamma_1 wedge gamma_1)   = N gamma_1 - gamma_1^2
    ==> sum_x lambda_2[px,qx]       = -(gamma_1 - gamma_1^2) = -gamma_1 eta_1

so the contraction of lambda_2 IS the particle-hole product, and since gamma_1 eta_1 >= 0
exactly when 0 <= gamma_1 <= 1, the identity ties lambda_2 directly to the Pauli bound.

Measured here: the identity on exact states, then its violation along MK_gamma30 at N=3,
alongside the smallest eigenvalue of -sum_x lambda_2, which the identity says should equal
the smallest eigenvalue of gamma_1 eta_1 and be non-negative.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
from qdiss.eom import fullci, gwt
from qdiss.eom.evaluate import assemble
from qdiss.eom.hamiltonian import random_V, random_h
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


def wedge(a, b):
    return np.einsum("pq,vw->pvqw", a, b) - np.einsum("pw,vq->pvqw", a, b)


def orbit(key, g1, l2, M):
    t = {("G", 1): g1, ("L", 2): l2}
    out = np.zeros((M,) * 6, dtype=complex)
    for s, sub, keys in ORB[key]:
        out = out + s * np.einsum(sub, *[t[k] for k in keys], optimize=True)
    return out


print("=== the identity on exact states ===")
print("%8s | %22s %14s" % ("(M,N)", "|sum_x lam2 + g1 eta1|", "|g1 eta1|"))
for M, N in ((4, 2), (6, 3), (8, 4)):
    rng = np.random.default_rng(5)
    psi = fullci.random_state(M, N, rng)
    g1 = fullci.rdm(psi, M, 1)
    l2 = fullci.rdm(psi, M, 2) - wedge(g1, g1)
    ge = g1 - g1 @ g1
    print("%8s | %22.3e %14.4f"
          % ("(%d,%d)" % (M, N), np.abs(np.einsum("pxqx->pq", l2) + ge).max(),
             np.abs(ge).max()))
print()

M, N = 6, 3
EYE = np.eye(M, dtype=complex)
rng = np.random.default_rng(7)
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


def rk4(f, y, dt):
    def add(a, b, c):
        return tuple(x + c * k for x, k in zip(a, b))
    k1 = f(y); k2 = f(add(y, k1, dt / 2)); k3 = f(add(y, k2, dt / 2)); k4 = f(add(y, k3, dt))
    return tuple(x + dt / 6 * (a + 2 * b + 2 * c + d)
                 for x, a, b, c, d in zip(y, k1, k2, k3, k4))


g1_0 = fullci.rdm(psi0, M, 1)
l2_0 = fullci.rdm(psi0, M, 2) - wedge(g1_0, g1_0)
st = (g1_0.copy(), l2_0.copy())
print("=== MK_gamma30 at M=6 N=3: does lambda_2 keep the identity? ===")
print("%6s | %12s | %11s %11s | %10s"
      % ("t", "identity err", "min -sum l2", "min g1 eta1", "n_min g1"))
dt, nsteps, every = 0.002, 500, 50
for step in range(nsteps + 1):
    if step % every == 0:
        g1, l2 = st
        ge = g1 - g1 @ g1
        ge = 0.5 * (ge + ge.conj().T)
        S = -np.einsum("pxqx->pq", l2)
        S = 0.5 * (S + S.conj().T)
        print("%6.2f | %12.3e | %11.5f %11.5f | %10.5f"
              % (step * dt, np.abs(S - ge).max(),
                 np.linalg.eigvalsh(S).min(), np.linalg.eigvalsh(ge).min(),
                 np.linalg.eigvalsh(0.5 * (g1 + g1.conj().T)).min()))
    if step < nsteps:
        st = rk4(rhs, st, dt)
print()
print("min g1 eta1 < 0 is exactly a Pauli violation. If min -sum_x lambda_2 goes negative")
print("BEFORE min g1 eta1 does, the cumulant contraction is the leading indicator.")
