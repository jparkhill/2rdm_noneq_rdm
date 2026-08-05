#!/usr/bin/env python3
"""wick_lambda_30 versus MK_lambda_30: are they the same theory?

Derivation methods, stated explicitly:

    both closures use the SAME gamma_2-dot expression. qdiss.eom.wick (plain vacuum Wick,
    bare moments) and qdiss.eom.gwt (MKNO generalized normal ordering) produce identical
    term lists, and each is validated against qdiss.eom.fullci at 1e-15 by
    qdiss.eom.validate. The gamma_1-dot equation contains no gamma_3 at all, so lambda_3 = 0
    does not enter it, and it too is identical between the derivations.

What differs is the closed system, i.e. which quantities are integrated:

    wick_lambda_30   state = gamma_2 alone.
                     gamma_2-dot from the Wick EOM with gamma_3 -> separable.
                     gamma_1(t) is whatever the contraction of gamma_2(t) gives.

    MK_lambda_30     state = (gamma_1, lambda_2), the MKNO-native variables.
                     gamma_1-dot is the EXACT rank-1 equation, needing no closure.
                     lambda_2-dot = gamma_2-dot|_{lambda_3=0} - d/dt(gamma_1 wedge gamma_1),
                     so the whole truncation sits in lambda_2-dot.

Ground truth is exact wavefunction propagation, psi(t) = exp(-iHt) psi(0), with H built by
qdiss.eom.fullci from the single convention in qdiss.eom.hamiltonian.

These are algebraically equivalent only if the truncated gamma_2-dot contracts to the exact
gamma_1-dot, which it does not: the separable gamma_3 violates
sum_x gamma_3[pvx,qwx] = (N-2) gamma_2. So they should diverge, and this measures by how much.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
from qdiss.eom import fullci, gwt, wick
from qdiss.eom.hamiltonian import random_V, random_h
from qdiss.symbolic import mkno

M, N = 6, 3
UP3, LO3 = ("p", "v", "x"), ("q", "w", "y")
SEP3 = [(s, ",".join("".join(u) + "".join(l) for n, r, u, l in f) + "->pvxqwy",
         [(n, r) for n, r, u, l in f])
        for s, f in mkno.gamma_to_cumulants(UP3, LO3)
        if not any(n == "L" and r == 3 for n, r, u, l in f)]
EYE = np.eye(M, dtype=complex)


def wedge(a, b):
    return np.einsum("pq,vw->pvqw", a, b) - np.einsum("pw,vq->pvqw", a, b)


def sep_gamma3(g1, l2):
    t = {("G", 1): g1, ("L", 2): l2}
    out = np.zeros((M,) * 6, dtype=complex)
    for s, sub, keys in SEP3:
        out = out + s * np.einsum(sub, *[t[k] for k in keys], optimize=True)
    return out


def contract1(g2):
    return np.einsum("pxqx->pq", g2) / (N - 1)


def tens(g1, g2, g3):
    return {("D", 1): EYE, ("G", 1): g1, ("G", 2): g2, ("G", 3): g3}


def rhs_wick(g2, h, V):
    """state gamma_2. gamma_1 and lambda_2 are derived from it at each step."""
    g1 = contract1(g2)
    l2 = g2 - wedge(g1, g1)
    return wick.gamma_dot(2, h, V, tens(g1, g2, sep_gamma3(g1, l2)), M)


def rhs_mk(state, h, V):
    """state (gamma_1, lambda_2). gamma_1-dot is exact; only lambda_2-dot is truncated."""
    g1, l2 = state
    g2 = wedge(g1, g1) + l2
    d1 = wick.gamma_dot(1, h, V, tens(g1, g2, None), M)
    d2 = wick.gamma_dot(2, h, V, tens(g1, g2, sep_gamma3(g1, l2)), M)
    dl2 = d2 - (wedge(d1, g1) + wedge(g1, d1))
    return (d1, dl2)


def rk4(f, y, dt):
    if isinstance(y, tuple):
        def add(a, b, c):
            return tuple(x + c * k for x, k in zip(a, b))
        k1 = f(y)
        k2 = f(add(y, k1, dt / 2))
        k3 = f(add(y, k2, dt / 2))
        k4 = f(add(y, k3, dt))
        return tuple(x + dt / 6 * (a + 2 * b + 2 * c + d)
                     for x, a, b, c, d in zip(y, k1, k2, k3, k4))
    k1 = f(y)
    k2 = f(y + dt / 2 * k1)
    k3 = f(y + dt / 2 * k2)
    k4 = f(y + dt * k3)
    return y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


rng = np.random.default_rng(2024)
h, V = random_h(M, rng), random_V(M, rng)
psi0 = fullci.random_state(M, N, rng)
H = fullci.build_H(h, V, M)
w, U = np.linalg.eigh(H)

g2_0 = fullci.rdm(psi0, M, 2)
g1_0 = fullci.rdm(psi0, M, 1)
l2_0 = g2_0 - wedge(g1_0, g1_0)
print("consistency of the initial condition:")
print("   |contract(gamma_2) - gamma_1| = %.3e" % np.abs(contract1(g2_0) - g1_0).max())
print("   |sep(gamma_3) - exact gamma_3| = %.4f"
      % np.abs(sep_gamma3(g1_0, l2_0) - fullci.rdm(psi0, M, 3)).max())
print()
print("term lists: wick and gwt identical? rank1 %s  rank2 %s"
      % (wick.terms(1) == gwt.terms(1), wick.terms(2) == gwt.terms(2)))
print()

dt, nsteps, every = 0.002, 500, 50
sw, sm = g2_0.copy(), (g1_0.copy(), l2_0.copy())
print("%6s | %11s %11s | %11s | %9s %9s %9s"
      % ("t", "wick vs FCI", "MK vs FCI", "wick vs MK", "n_min FCI", "n_min W", "n_min MK"))
for step in range(nsteps + 1):
    if step % every == 0:
        t = step * dt
        psit = U @ (np.exp(-1j * w * t) * (U.conj().T @ psi0))
        g1_ex = fullci.rdm(psit, M, 1)
        g1_w = contract1(sw)
        g1_m = sm[0]
        print("%6.2f | %11.3e %11.3e | %11.3e | %9.5f %9.5f %9.5f"
              % (t, np.abs(g1_w - g1_ex).max(), np.abs(g1_m - g1_ex).max(),
                 np.abs(g1_w - g1_m).max(),
                 np.linalg.eigvalsh(0.5 * (g1_ex + g1_ex.conj().T)).min(),
                 np.linalg.eigvalsh(0.5 * (g1_w + g1_w.conj().T)).min(),
                 np.linalg.eigvalsh(0.5 * (g1_m + g1_m.conj().T)).min()))
    if step < nsteps:
        sw = rk4(lambda y: rhs_wick(y, h, V), sw, dt)
        sm = rk4(lambda y: rhs_mk(y, h, V), sm, dt)
print()
print("n_min is the smallest natural occupation; below zero violates the Pauli bound.")
