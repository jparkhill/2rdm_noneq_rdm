#!/usr/bin/env python3
"""Is MK_lambda_30 exact for two electrons, and if not, is the cause lambda_3 or a bug?

Derivation: gamma_2-dot and gamma_1-dot from qdiss.eom.wick (plain vacuum Wick, bare
moments), identical to qdiss.eom.gwt, both validated against qdiss.eom.fullci at 1e-15.
Ground truth here is exact wavefunction propagation.

For N = 2 the exact gamma_3 vanishes identically, since three particles cannot be
annihilated from a two-particle state. The cumulant decomposition then reads

    0 = gamma_3 = A + B + lambda_3,        A = the gamma_1^3 orbit, B = the gamma_1*lambda_2 orbit

so lambda_3 = -(A + B). The lambda_3 = 0 closure is exact at N = 2 if and only if A + B = 0,
which is a nontrivial identity relating the two orbits.

This matters twice over. It decides whether MK_lambda_30 is exact for two electrons, and it
is the first test of the RELATIVE COEFFICIENT of the two orbits: the earlier validation used
a Slater determinant, where lambda_2 = 0, so it pinned A and said nothing about B. If
A + B is nonzero but A + c*B vanishes for some c other than 1, the expansion carries a wrong
coefficient and everything built on lambda_3 = 0 is affected.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
from qdiss.eom import fullci, wick
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


def orbit(key, g1, l2, M):
    t = {("G", 1): g1, ("L", 2): l2}
    out = np.zeros((M,) * 6, dtype=complex)
    for s, sub, keys in ORB[key]:
        out = out + s * np.einsum(sub, *[t[k] for k in keys], optimize=True)
    return out


def wedge(a, b):
    return np.einsum("pq,vw->pvqw", a, b) - np.einsum("pw,vq->pvqw", a, b)


print("=== is lambda_3 zero for two electrons? ===")
print("%8s | %10s %10s %10s %10s | %12s %10s"
      % ("(M,N)", "|g3 exact|", "|A|", "|B|", "|A+B|", "best c", "|A+cB|"))
for M, N in ((4, 2), (6, 2), (8, 2), (6, 3)):
    rng = np.random.default_rng(99)
    psi = fullci.random_state(M, N, rng)
    g1 = fullci.rdm(psi, M, 1)
    g2 = fullci.rdm(psi, M, 2)
    g3 = fullci.rdm(psi, M, 3)
    l2 = g2 - wedge(g1, g1)
    A = orbit("A", g1, l2, M)
    B = orbit("B", g1, l2, M)
    av, bv = A.ravel(), B.ravel()
    c = -np.vdot(bv, av) / np.vdot(bv, bv)
    print("%8s | %10.3e %10.4f %10.4f %10.3e | %12.6f %10.3e"
          % ("(%d,%d)" % (M, N), np.abs(g3).max(), np.abs(A).max(), np.abs(B).max(),
             np.abs(A + B).max(), c.real, np.abs(A + c * B).max()))
print()
print("For N=2 the exact gamma_3 is zero, so |A+B| IS |lambda_3|. A nonzero best c other")
print("than 1.0 would indicate a wrong relative coefficient rather than genuine lambda_3.")
print()

print("=== propagation at N = 2: does MK_lambda_30 reproduce exact dynamics? ===")
M, N = 4, 2
EYE = np.eye(M, dtype=complex)
rng = np.random.default_rng(7)
h, V = random_h(M, rng), random_V(M, rng)
psi0 = fullci.random_state(M, N, rng)
H = fullci.build_H(h, V, M)
w, U = np.linalg.eigh(H)


def sep3(g1, l2):
    return orbit("A", g1, l2, M) + orbit("B", g1, l2, M)


def tens(g1, g2, g3):
    return {("D", 1): EYE, ("G", 1): g1, ("G", 2): g2, ("G", 3): g3}


def rhs_mk(state, closure):
    g1, l2 = state
    g2 = wedge(g1, g1) + l2
    g3 = sep3(g1, l2) if closure == "l3" else np.zeros((M,) * 6, dtype=complex)
    d1 = wick.gamma_dot(1, h, V, tens(g1, g2, None), M)
    d2 = wick.gamma_dot(2, h, V, tens(g1, g2, g3), M)
    return (d1, d2 - (wedge(d1, g1) + wedge(g1, d1)))


def rk4(f, y, dt):
    def add(a, b, c):
        return tuple(x + c * k for x, k in zip(a, b))
    k1 = f(y); k2 = f(add(y, k1, dt / 2)); k3 = f(add(y, k2, dt / 2)); k4 = f(add(y, k3, dt))
    return tuple(x + dt / 6 * (a + 2 * b + 2 * c + d)
                 for x, a, b, c, d in zip(y, k1, k2, k3, k4))


g1_0 = fullci.rdm(psi0, M, 1)
l2_0 = fullci.rdm(psi0, M, 2) - wedge(g1_0, g1_0)
st = {k: (g1_0.copy(), l2_0.copy()) for k in ("l3", "g3")}
dt, nsteps, every = 0.002, 500, 50
print("%6s | %14s %14s | %14s %14s"
      % ("t", "MK_l3=0 g1 err", "MK_g3=0 g1 err", "MK_l3=0 g2 err", "MK_g3=0 g2 err"))
for step in range(nsteps + 1):
    if step % every == 0:
        t = step * dt
        psit = U @ (np.exp(-1j * w * t) * (U.conj().T @ psi0))
        e1, e2 = fullci.rdm(psit, M, 1), fullci.rdm(psit, M, 2)
        row = []
        for k in ("l3", "g3"):
            a, b = st[k]
            row.append((np.abs(a - e1).max(), np.abs(wedge(a, a) + b - e2).max()))
        print("%6.2f | %14.3e %14.3e | %14.3e %14.3e"
              % (t, row[0][0], row[1][0], row[0][1], row[1][1]))
    if step < nsteps:
        for k in ("l3", "g3"):
            st[k] = rk4(lambda y, kk=k: rhs_mk(y, kk), st[k], dt)
