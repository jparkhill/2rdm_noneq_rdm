#!/usr/bin/env python3
"""Is the MK_gamma30 residual at N=2 integrator error or model error?

RK4 has global error O(dt^4), so halving dt must reduce the residual by about 16 if the
closure is exact. A plateau would mean genuine model error. Method: qdiss.eom.gwt cumulant
form, state (gamma_1, lambda_2), lambda_3 solved from gamma_3 = 0. Truth: exact propagation.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
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


M, N, T = 6, 2, 0.5
EYE = np.eye(M, dtype=complex)
rng = np.random.default_rng(7)
h, V = random_h(M, rng), random_V(M, rng)
psi0 = fullci.random_state(M, N, rng)
H = fullci.build_H(h, V, M)
w, U = np.linalg.eigh(H)
psiT = U @ (np.exp(-1j * w * T) * (U.conj().T @ psi0))
e1 = fullci.rdm(psiT, M, 1)


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
print("MK_gamma30, M=%d N=%d, error in gamma_1 at t=%.2f" % (M, N, T))
print("%10s %8s | %13s | %8s" % ("dt", "steps", "err", "ratio"))
prev = None
for dt in (0.008, 0.004, 0.002, 0.001):
    st = (g1_0.copy(), l2_0.copy())
    for _ in range(int(round(T / dt))):
        st = rk4(rhs, st, dt)
    err = float(np.abs(st[0] - e1).max())
    print("%10.4f %8d | %13.4e | %8s"
          % (dt, int(round(T / dt)), err, "-" if prev is None else "%.1f" % (prev / err)))
    prev = err
print()
print("A ratio near 16 per halving is RK4 convergence, i.e. the closure itself is exact.")
