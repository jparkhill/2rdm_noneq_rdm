#!/usr/bin/env python3
"""The hierarchy of multiple-trace conditions on a cumulant.

Correction to an earlier overstatement. The single-trace identity
sum_x lambda_k[Px,Qx] = (expression in lower cumulants) exists and is computable at every rank;
it never "breaks". tools/rank4_contraction.py tested ONE candidate form at rank 4,
c (E lambda_3 + lambda_3 E), found c = 0, and wrongly concluded the rule fails. All that shows
is that the rank-4 rule has more terms than that one, which is expected since rank 4 is the
first order where lambda_2 lambda_2 partitions appear.

John's point: contracting again gives a weaker but simpler condition at lower rank. Tracing
lambda_3 over two index pairs instead of one reduces a 4-index condition to a 2-index one, and
that reduced condition is expressed entirely in the variables already propagated in the
(gamma_1, lambda_2) scheme.

Contracting the verified single-trace identity
    sum_y lambda_3[pvy,qwy] = E lambda_2 + lambda_2 E ,   E = gamma_1 - 1/2
a second time over v = w gives, with S = sum_x lambda_2[.x,.x] = -(gamma_1 - gamma_1^2),

    sum_{x,y} lambda_3[pxy,qxy] = E S + S E + G + G' - S

    G[p,q]  = sum_{v,t} gamma_1[v,t] lambda_2[pt,qv]
    G'[p,q] = sum_{v,t} lambda_2[pv,qt] gamma_1[t,v]

the two cross-contractions arising from E acting on the slot being traced. Each is checked
numerically here rather than trusted, and the reconstructions on record are then scored against
the double-trace condition, which is an M^2 condition rather than an M^4 one.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from math import comb
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


def orbit(key, a, X, M):
    t = {("G", 1): a, ("L", 2): X}
    out = np.zeros((M,) * 6, dtype=complex)
    for s, sb, keys in ORB[key]:
        out = out + s * np.einsum(sb, *[t[k] for k in keys], optimize=True)
    return out


def invert(R, M):
    I = np.eye(M, dtype=complex)
    S = np.einsum("pxqx->pq", R)
    tr = np.trace(S) / (3 * M - 6)
    SX = (S - tr * I) / (2 * M - 6)
    return (R - wedge(I, SX) - wedge(SX, I)) / (M - 4)


print("%8s | %14s %14s | %14s %14s"
      % ("(M,N)", "single trace", "double trace", "|1-trace tgt|", "|2-trace tgt|"))
for M, N in ((6, 3), (8, 3), (6, 4), (8, 4)):
    rng = np.random.default_rng(9)
    psi = fullci.random_state(M, N, rng)
    g1 = fullci.rdm(psi, M, 1)
    g2 = fullci.rdm(psi, M, 2)
    g3 = fullci.rdm(psi, M, 3)
    l2 = g2 - wedge(g1, g1)
    l3 = g3 - orbit("A", g1, l2, M) - orbit("B", g1, l2, M)
    I = np.eye(M, dtype=complex)
    E = g1 - 0.5 * I
    S = np.einsum("pxqx->pq", l2)
    one = (np.einsum("py,yvqw->pvqw", E, l2) + np.einsum("vy,pyqw->pvqw", E, l2)
           + np.einsum("pvyw,yq->pvqw", l2, E) + np.einsum("pvqy,yw->pvqw", l2, E))
    G = np.einsum("vt,ptqv->pq", g1, l2)
    Gp = np.einsum("pvqt,tv->pq", l2, g1)
    two = E @ S + S @ E + G + Gp - S
    e1 = np.abs(np.einsum("pvyqwy->pvqw", l3) - one).max()
    e2 = np.abs(np.einsum("pxyqxy->pq", l3) - two).max()
    print("%8s | %14.3e %14.3e | %14.4f %14.4f"
          % ("(%d,%d)" % (M, N), e1, e2, np.abs(one).max(), np.abs(two).max()))
print()
print("also checking that S itself reduces: sum_x lambda_2 = -(g1 - g1^2)")

M, N = 6, 3
rng = np.random.default_rng(9)
psi = fullci.random_state(M, N, rng)
g1 = fullci.rdm(psi, M, 1)
g2 = fullci.rdm(psi, M, 2)
g3 = fullci.rdm(psi, M, 3)
l2 = g2 - wedge(g1, g1)
l3x = g3 - orbit("A", g1, l2, M) - orbit("B", g1, l2, M)
I = np.eye(M, dtype=complex)
E = g1 - 0.5 * I
S = np.einsum("pxqx->pq", l2)
print("   residual %.3e" % np.abs(S + (g1 - g1 @ g1)).max())
one = (np.einsum("py,yvqw->pvqw", E, l2) + np.einsum("vy,pyqw->pvqw", E, l2)
       + np.einsum("pvyw,yq->pvqw", l2, E) + np.einsum("pvqy,yw->pvqw", l2, E))
two = (E @ S + S @ E + np.einsum("vt,ptqv->pq", g1, l2)
       + np.einsum("pvqt,tv->pq", l2, g1) - S)
l3_tc = orbit("B", I, invert(one, M), M)
print()
print("reconstructions scored against BOTH conditions, M=6 N=3")
print("%-24s | %13s %13s" % ("lambda_3", "single trace", "double trace"))
for nm, l3 in (("exact", l3x),
               ("zero", np.zeros((M,) * 6, dtype=complex)),
               ("-(A+B)  MK_gamma30", -(orbit("A", g1, l2, M) + orbit("B", g1, l2, M))),
               ("closed-form trace-exact", l3_tc)):
    print("%-24s | %13.3e %13.3e"
          % (nm, np.abs(np.einsum("pvyqwy->pvqw", l3) - one).max(),
             np.abs(np.einsum("pxyqxy->pq", l3) - two).max()))
print()
print("The double trace is implied by the single one, so it cannot separate a reconstruction")
print("that already satisfies the single trace. Its use is as a CHEAPER PARTIAL condition:")
print("M^2 numbers instead of M^4, in the variables already propagated.")
