#!/usr/bin/env python3
"""Where does the Pauli violation come from, given the blocking factors are present?

The MKNO 1-RDM EOM splits exactly into a lambda_2-free part and a lambda_2-carrying part:
factors D1 x14, G1 x18, L2 x4, so 4 of the 12 two-body terms carry lambda_2.

The lambda_2-free part is the same EOM with lambda_2 := 0, which is TDHF, already identified
as i[(h+W)^T, gamma_1] at 3.9e-16 -- a commutator, hence isospectral, hence it moves NO
occupation and can never violate Pauli. That is the blocking factors doing their job.

So every occupation change, and therefore the entire Pauli violation, must come from the 4
lambda_2 terms. lambda_2 is an independent propagated tensor carrying no gamma-eta structure
of its own, so nothing in the terms blocks it. Its only protection is the contraction
identity sum_x lambda_2[px,qx] = -gamma_1 eta_1, which MK_gamma30 does not enforce.

Measured along the MK_gamma30 trajectory: the occupation velocity split into the two parts.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
from qdiss.eom import fullci, gwt
from qdiss.eom.evaluate import assemble
from qdiss.eom.hamiltonian import random_V, random_h
from qdiss.symbolic import mkno

M, N = 6, 3
EYE = np.eye(M, dtype=complex)
UP3, LO3 = ("p", "v", "x"), ("q", "w", "y")
ORB = {"A": [], "B": []}
for s, f in mkno.gamma_to_cumulants(UP3, LO3):
    if any(n == "L" and r == 3 for n, r, u, l in f):
        continue
    k = "B" if any(n == "L" and r == 2 for n, r, u, l in f) else "A"
    ORB[k].append((s, ",".join("".join(u) + "".join(l) for n, r, u, l in f) + "->pvxqwy",
                   [(n, r) for n, r, u, l in f]))
CUM1, CUM2 = gwt.terms_cumulant(1), gwt.terms_cumulant(2)
HAS_L2 = [t for t in CUM1[1] if any(n == "L" and r == 2 for n, r, u, l in t[1])]
NO_L2 = [t for t in CUM1[1] if not any(n == "L" and r == 2 for n, r, u, l in t[1])]
print("rank-1 two-body terms: %d carrying lambda_2, %d not" % (len(HAS_L2), len(NO_L2)))


def wedge(a, b):
    return np.einsum("pq,vw->pvqw", a, b) - np.einsum("pw,vq->pvqw", a, b)


def orbit(key, g1, l2):
    t = {("G", 1): g1, ("L", 2): l2}
    out = np.zeros((M,) * 6, dtype=complex)
    for s, sub, keys in ORB[key]:
        out = out + s * np.einsum(sub, *[t[k] for k in keys], optimize=True)
    return out


rng = np.random.default_rng(7)
h, V = random_h(M, rng), random_V(M, rng)
psi0 = fullci.random_state(M, N, rng)


def parts(g1, l2):
    t = {("D", 1): EYE, ("G", 1): g1, ("L", 2): l2}
    free = assemble(CUM1[0], NO_L2, h, V, t, "pq", M)
    lam = assemble([], HAS_L2, h, V, t, "pq", M)
    return free, lam


def rhs(state):
    g1, l2 = state
    l3 = -(orbit("A", g1, l2) + orbit("B", g1, l2))
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
print()
print("occupation velocity ndot_p = <p|gamma_1-dot|p> in the natural-orbital basis")
print("%6s | %10s | %14s %14s | %12s"
      % ("t", "n_min", "max|ndot| free", "max|ndot| lam2", "identity err"))
dt, nsteps, every = 0.002, 400, 50
for step in range(nsteps + 1):
    if step % every == 0:
        g1, l2 = st
        lam, U = np.linalg.eigh(0.5 * (g1 + g1.conj().T))
        free, l2p = parts(g1, l2)
        nf = np.abs(np.real(np.diag(U.conj().T @ free @ U))).max()
        nl = np.abs(np.real(np.diag(U.conj().T @ l2p @ U))).max()
        ident = np.abs(-np.einsum("pxqx->pq", l2) - (g1 - g1 @ g1)).max()
        print("%6.2f | %10.5f | %14.3e %14.3e | %12.3e"
              % (step * dt, lam.min(), nf, nl, ident))
    if step < nsteps:
        st = rk4(rhs, st, dt)
print()
print("The lambda_2-free part moves no occupation at all: that IS the blocking working.")
print("Everything that moves n comes from the 4 lambda_2 terms, which carry no blocking")
print("factor because lambda_2 is an independent tensor, not a product of gamma and eta.")
