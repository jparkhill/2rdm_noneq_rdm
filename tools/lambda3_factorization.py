#!/usr/bin/env python3
"""Factor the rank-3 trace identity into particle/hole form.

Rank 2 factors into a particle-hole product: sum_x lambda_2[px,qx] = -(g1 - g1^2) = -g1 eta.
Substituting the verified closed forms for sum_y A and sum_y B into

    sum_y lambda_3 = (N-2) gamma_2 - sum_y A - sum_y B

and using gamma_2 = g1 wedge g1 + lambda_2, the N-dependence cancels and, by hand, so do the
two groups involving G = -(g1 - g1^2):

    from (N-2)gamma_2 and sum_y A :  -2(g1^g1) + (g1^g1^2) + (g1^2^g1) = -(g1 ^ g1 eta) - (g1 eta ^ g1)
    from -(g1 ^ G) - (G ^ g1)      :  +(g1 ^ g1 eta) + (g1 eta ^ g1)

leaving only the lambda_2 pieces. With Gam = g1 (+) g1 the one-body operator on pair space,
that predicts

    sum_y lambda_3 = E lambda_2 + lambda_2 E ,      E = Gam - 1

and since the pair-space hole operator is Eta = eta (+) eta = 2 - Gam, E = (Gam - Eta)/2, so

    sum_y lambda_3 = (1/2) [ (Gam - Eta) lambda_2 + lambda_2 (Gam - Eta) ]

a particle-minus-hole anticommutator rather than a particle-times-hole product. Hand algebra,
tested here against the exact contraction.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
from qdiss.eom import fullci
from qdiss.symbolic import mkno

UP3, LO3 = ("p", "v", "x"), ("q", "w", "y")
ORB = {"A": [], "B": []}
for s, f in mkno.gamma_to_cumulants(UP3, LO3):
    if any(n == "L" and r == 3 for n, r, u, l in f):
        continue
    k = "B" if any(n == "L" and r == 2 for n, r, u, l in f) else "A"
    ORB[k].append((s, ",".join("".join(u) + "".join(l) for n, r, u, l in f) + "->pvxqwy",
                   [(n, r) for n, r, u, l in f]))


def wedge(a, b):
    return np.einsum("pq,vw->pvqw", a, b) - np.einsum("pw,vq->pvqw", a, b)


def orbit(key, g1, l2, M):
    t = {("G", 1): g1, ("L", 2): l2}
    out = np.zeros((M,) * 6, dtype=complex)
    for s, sub, keys in ORB[key]:
        out = out + s * np.einsum(sub, *[t[k] for k in keys], optimize=True)
    return out


def left(X, l2):
    """(X (+) X) acting on the upper index pair of lambda_2."""
    return (np.einsum("py,yvqw->pvqw", X, l2) + np.einsum("vy,pyqw->pvqw", X, l2))


def right(l2, X):
    """(X (+) X) acting on the lower index pair of lambda_2."""
    return (np.einsum("pvyw,yq->pvqw", l2, X) + np.einsum("pvqy,yw->pvqw", l2, X))


print("%8s | %14s %14s %14s | %11s"
      % ("(M,N)", "E l2 + l2 E", "(Gam-Eta)/2 form", "-2 l2 only", "|sum_y l3|"))
for M, N in ((4, 2), (6, 3), (8, 4), (10, 4), (10, 5)):
    rng = np.random.default_rng(5)
    psi = fullci.random_state(M, N, rng)
    g1 = fullci.rdm(psi, M, 1)
    g2 = fullci.rdm(psi, M, 2)
    g3 = fullci.rdm(psi, M, 3)
    l2 = g2 - wedge(g1, g1)
    l3 = g3 - orbit("A", g1, l2, M) - orbit("B", g1, l2, M)
    tr3 = np.einsum("pvyqwy->pvqw", l3)
    I = np.eye(M, dtype=complex)
    eta = I - g1
    E = g1 - 0.5 * I
    predE = left(E, l2) + right(l2, E)
    predPH = 0.5 * (left(g1 - eta, l2) + right(l2, g1 - eta))
    print("%8s | %14.3e %14.3e %14.3e | %11.4f"
          % ("(%d,%d)" % (M, N), np.abs(tr3 - predE).max(),
             np.abs(tr3 - predPH).max(), np.abs(tr3 + 2 * l2).max(),
             np.abs(tr3).max()))
print()
print("E = g1 - I/2 per index, so (E (+) E) = Gam - 1 on the pair. The -2 l2 column is the")
print("same expression with the gamma_1 contractions dropped, shown for scale.")
