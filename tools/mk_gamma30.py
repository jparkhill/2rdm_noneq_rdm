#!/usr/bin/env python3
"""MK_gamma30: the MKNO cumulant EOM closed by the identity gamma_3 = 0.

Derivation method: qdiss.eom.gwt, MKNO generalized normal ordering and the generalized
Wick theorem, in CUMULANT form -- 12 one-body and 140 two-body terms over gamma_1,
lambda_2 and lambda_3. Validated against qdiss.eom.fullci at 1e-15 by qdiss.eom.validate.
Ground truth here is exact wavefunction propagation.

Propagated state: (gamma_1, lambda_2), the MKNO-native variables.

    gamma_1-dot   the exact rank-1 cumulant EOM. Contains no lambda_3, so no closure.
    lambda_2-dot  gamma_2-dot - d/dt(gamma_1 wedge gamma_1)

Closure. For two electrons gamma_3 vanishes identically, and the cumulant decomposition
reads gamma_3 = A + B + lambda_3 with A the gamma_1^3 orbit and B the gamma_1*lambda_2
orbit. Rather than setting lambda_3 = 0, which is false at N = 2 and injects a spurious
gamma_3 of about 4.6e-2, SOLVE the identity:

    lambda_3 := -(A + B)     so that the implied gamma_3 is exactly zero

This is a genuinely different code path from zeroing gamma_3 in the bare-moment EOM: it
drives the 140-term cumulant expression with a nonzero lambda_3 tensor, which no previous
propagation has done.

At N = 2 the resulting system is closed and exact, so it must reproduce FullCI to
integrator error. At N = 3 gamma_3 is not zero, so the same closure must fail -- included
as a control so the test is known to have teeth.
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
CUM1 = gwt.terms_cumulant(1)
CUM2 = gwt.terms_cumulant(2)
print("MKNO cumulant term lists: rank1 %d + %d, rank2 %d + %d"
      % (len(CUM1[0]), len(CUM1[1]), len(CUM2[0]), len(CUM2[1])))
print("lambda_3 factors present in the rank-2 list: %d"
      % sum(1 for s, f in CUM2[1] for n, r, u, l in f if n == "L" and r == 3))
print()


def wedge(a, b):
    return np.einsum("pq,vw->pvqw", a, b) - np.einsum("pw,vq->pvqw", a, b)


def orbit(key, g1, l2, M):
    t = {("G", 1): g1, ("L", 2): l2}
    out = np.zeros((M,) * 6, dtype=complex)
    for s, sub, keys in ORB[key]:
        out = out + s * np.einsum(sub, *[t[k] for k in keys], optimize=True)
    return out


def rk4(f, y, dt):
    def add(a, b, c):
        return tuple(x + c * k for x, k in zip(a, b))
    k1 = f(y); k2 = f(add(y, k1, dt / 2)); k3 = f(add(y, k2, dt / 2)); k4 = f(add(y, k3, dt))
    return tuple(x + dt / 6 * (a + 2 * b + 2 * c + d)
                 for x, a, b, c, d in zip(y, k1, k2, k3, k4))


def run(M, N, label, dt=0.002, nsteps=500, every=100, seed=7):
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
    l3_0 = -(orbit("A", g1_0, l2_0, M) + orbit("B", g1_0, l2_0, M))
    implied = orbit("A", g1_0, l2_0, M) + orbit("B", g1_0, l2_0, M) + l3_0
    print("=== %s : M=%d N=%d ===" % (label, M, N))
    print("   implied gamma_3 from the closure : %.3e" % np.abs(implied).max())
    print("   exact gamma_3                    : %.3e"
          % np.abs(fullci.rdm(psi0, M, 3)).max())
    print("   |lambda_3 used|                  : %.4f" % np.abs(l3_0).max())
    st = (g1_0.copy(), l2_0.copy())
    print("   %6s | %13s %13s" % ("t", "gamma_1 err", "gamma_2 err"))
    for step in range(nsteps + 1):
        if step % every == 0:
            t = step * dt
            psit = U @ (np.exp(-1j * w * t) * (U.conj().T @ psi0))
            e1, e2 = fullci.rdm(psit, M, 1), fullci.rdm(psit, M, 2)
            g1, l2 = st
            print("   %6.2f | %13.3e %13.3e"
                  % (t, np.abs(g1 - e1).max(), np.abs(wedge(g1, g1) + l2 - e2).max()))
        if step < nsteps:
            st = rk4(rhs, st, dt)
    print()


run(4, 2, "MK_gamma30, two electrons")
run(6, 2, "MK_gamma30, two electrons, larger basis")
run(6, 3, "MK_gamma30, THREE electrons (control: must fail)")
print("At N=2 the closure is an identity, so the system is closed and exact and the residual")
print("is integrator error only. At N=3 gamma_3 is not zero and the same closure is wrong.")
