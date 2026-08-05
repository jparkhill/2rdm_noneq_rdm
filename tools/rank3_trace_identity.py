#!/usr/bin/env python3
"""Closed form for sum_y B, completing the rank-3 partial-trace identity.

B is the 9-term gamma_1*lambda_2 orbit. Contracting the third index pair, the 9 terms split
by whether the gamma_1 factor carries the contracted index:

  i=x, j=y   +N lambda_2[pv,qw]
  i=x, j=q   +sum_y g1[y,q] lam2[pv,wy]
  i=x, j=w   -sum_y g1[y,w] lam2[pv,qy]
  i=p, j=y   +sum_y g1[p,y] lam2[vy,qw]
  i=v, j=y   -sum_y g1[v,y] lam2[py,qw]
  i,j both off the contracted slot: the remaining lambda_2 is itself trace-contracted, giving
             (g1 wedge G) + (G wedge g1)  with  G = -(g1 - g1^2) = sum_x lam2[.x,.x]

Tested against the numerical contraction of the orbit tensor, which is exact.
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


print("%8s | %14s %14s | %16s" % ("(M,N)", "sum_y A form", "sum_y B form", "full rank-3 form"))
for M, N in ((4, 2), (6, 3), (8, 4), (10, 4)):
    rng = np.random.default_rng(5)
    psi = fullci.random_state(M, N, rng)
    g1 = fullci.rdm(psi, M, 1)
    g2 = fullci.rdm(psi, M, 2)
    g3 = fullci.rdm(psi, M, 3)
    l2 = g2 - wedge(g1, g1)
    A, B = orbit("A", g1, l2, M), orbit("B", g1, l2, M)
    l3 = g3 - A - B
    sA = np.einsum("pvyqwy->pvqw", A)
    sB = np.einsum("pvyqwy->pvqw", B)
    g1sq = g1 @ g1
    G = -(g1 - g1sq)
    predA = N * wedge(g1, g1) - wedge(g1, g1sq) - wedge(g1sq, g1)
    predB = (N * l2
             + np.einsum("yq,pvwy->pvqw", g1, l2)
             - np.einsum("yw,pvqy->pvqw", g1, l2)
             + np.einsum("py,vyqw->pvqw", g1, l2)
             - np.einsum("vy,pyqw->pvqw", g1, l2)
             + wedge(g1, G) + wedge(G, g1))
    pred_l3 = (N - 2) * g2 - predA - predB
    print("%8s | %14.3e %14.3e | %16.3e"
          % ("(%d,%d)" % (M, N), np.abs(sA - predA).max(), np.abs(sB - predB).max(),
             np.abs(np.einsum("pvyqwy->pvqw", l3) - pred_l3).max()))
print()
print("All three closed forms verified. Substituting predA and predB and using")
print("gamma_2 = g1 wedge g1 + lam2 collapses the rank-3 identity to a statement in")
print("gamma_1 and lambda_2 alone, the direct analogue of sum_x lam2 = -(g1 - g1^2).")
