#!/usr/bin/env python3
"""Independent validation of gamma_to_cumulants at rank 3.

Every closure result rests on this one function, and it has only ever been used with the
lambda_3 block deleted, never checked against anything independent.

Two non-circular tests.

1. Single determinant. All cumulants vanish exactly for a Slater determinant, so the
   gamma_1^3 orbit ALONE must reproduce gamma_3. This pins the 6 terms and their signs
   with no reference to the cumulant definition.

2. Antisymmetry. gamma_3 is antisymmetric under exchange of any two upper indices and of
   any two lower indices. Each orbit tensor must be separately antisymmetric, since the
   orbits are closed under the permutation group. A wrong sign anywhere in the set-partition
   expansion breaks this, and it is independent of any numerical reference.

Also reports the rank-2 expansion, which should be the familiar gamma_2 = gamma_1 wedge
gamma_1 + lambda_2, as a sanity anchor.
"""
import itertools
import sys
import numpy as np
sys.path.insert(0, ".")
from qdiss.symbolic import mkno

M = 6
UP2, LO2 = ("p", "v"), ("q", "w")
UP3, LO3 = ("p", "v", "x"), ("q", "w", "y")

print("=== rank-2 expansion (anchor) ===")
for s, f in mkno.gamma_to_cumulants(UP2, LO2):
    print("   %+d  %s" % (s, " ".join("%s%d[%s;%s]"
          % ({"G": "g", "L": "lam"}[n], r, "".join(u), "".join(l)) for n, r, u, l in f)))
print()

ORB = {"g1g1g1": [], "g1l2": [], "l3": []}
for s, f in mkno.gamma_to_cumulants(UP3, LO3):
    if any(n == "L" and r == 3 for n, r, u, l in f):
        k = "l3"
    elif any(n == "L" and r == 2 for n, r, u, l in f):
        k = "g1l2"
    else:
        k = "g1g1g1"
    ORB[k].append((s, ",".join("".join(u) + "".join(l) for n, r, u, l in f) + "->pvxqwy",
                   [(n, r) for n, r, u, l in f]))
print("rank-3 orbit sizes: gamma_1^3 %d, gamma_1*lambda_2 %d, lambda_3 %d   (total %d)"
      % (len(ORB["g1g1g1"]), len(ORB["g1l2"]), len(ORB["l3"]),
         sum(len(v) for v in ORB.values())))
print()


def ad(p, v):
    out = np.zeros_like(v); bit = 1 << p
    for s in np.nonzero(v)[0]:
        if s & bit: continue
        out[s | bit] += ((-1.0) ** bin(s & (bit - 1)).count("1")) * v[s]
    return out


def aa(p, v):
    out = np.zeros_like(v); bit = 1 << p
    for s in np.nonzero(v)[0]:
        if not (s & bit): continue
        out[s ^ bit] += ((-1.0) ** bin(s & (bit - 1)).count("1")) * v[s]
    return out


def build(key, tens):
    out = np.zeros((M,) * 6, dtype=complex)
    for s, sub, keys in ORB[key]:
        out = out + s * np.einsum(sub, *[tens[k] for k in keys], optimize=True)
    return out


print("=== 1. single determinant: cumulants vanish, gamma_1^3 orbit must give gamma_3 ===")
dim = 1 << M
det = np.zeros(dim, dtype=complex)
det[0b000111] = 1.0
Q1 = np.array([aa(p, det) for p in range(M)])
g1 = Q1.conj() @ Q1.T
Q2 = np.array([aa(v, aa(p, det)) for p in range(M) for v in range(M)])
g2 = (Q2.conj() @ Q2.T).reshape(M, M, M, M)
Q3 = np.array([aa(x, aa(v, aa(p, det)))
               for p in range(M) for v in range(M) for x in range(M)])
g3 = (Q3.conj() @ Q3.T).reshape((M,) * 6)
l2 = g2 - (np.einsum("pq,vw->pvqw", g1, g1) - np.einsum("pw,vq->pvqw", g1, g1))
print("   |lambda_2| for the determinant                 = %.3e   (must be 0)" % np.abs(l2).max())
A = build("g1g1g1", {("G", 1): g1, ("L", 2): l2})
print("   |gamma_1^3 orbit - exact gamma_3|              = %.3e" % np.abs(A - g3).max())
print("   scale |gamma_3| = %.4f" % np.abs(g3).max())
print()

print("=== 2. antisymmetry of each orbit, on a correlated state ===")
r = np.random.default_rng(5)
popc = np.array([bin(s).count("1") for s in range(dim)])
psi = np.where(popc == 3, r.normal(size=dim) + 1j * r.normal(size=dim), 0.0)
psi = psi / np.linalg.norm(psi)
Q1 = np.array([aa(p, psi) for p in range(M)])
g1c = Q1.conj() @ Q1.T
Q2 = np.array([aa(v, aa(p, psi)) for p in range(M) for v in range(M)])
g2c = (Q2.conj() @ Q2.T).reshape(M, M, M, M)
Q3 = np.array([aa(x, aa(v, aa(p, psi)))
               for p in range(M) for v in range(M) for x in range(M)])
g3c = (Q3.conj() @ Q3.T).reshape((M,) * 6)
l2c = g2c - (np.einsum("pq,vw->pvqw", g1c, g1c) - np.einsum("pw,vq->pvqw", g1c, g1c))
tens = {("G", 1): g1c, ("L", 2): l2c}
print("%12s | %14s %14s | %s" % ("orbit", "swap upper 0,1", "swap lower 0,1", "scale"))
for key in ("g1g1g1", "g1l2"):
    T = build(key, tens)
    su = np.abs(T + np.transpose(T, (1, 0, 2, 3, 4, 5))).max()
    sl = np.abs(T + np.transpose(T, (0, 1, 2, 4, 3, 5))).max()
    print("%12s | %14.3e %14.3e | %.4f" % (key, su, sl, np.abs(T).max()))
T = g3c
print("%12s | %14.3e %14.3e | %.4f"
      % ("exact g3", np.abs(T + np.transpose(T, (1, 0, 2, 3, 4, 5))).max(),
         np.abs(T + np.transpose(T, (0, 1, 2, 4, 3, 5))).max(), np.abs(T).max()))
print()
l3 = g3c - build("g1g1g1", tens) - build("g1l2", tens)
print("=== implied lambda_3 = gamma_3 - orbits, and its contraction ===")
print("   |lambda_3| = %.4f" % np.abs(l3).max())
print("   |sum_x lambda_3[pvx,qwx]| = %.4f     |lambda_2| = %.4f"
      % (np.abs(np.einsum("pvxqwx->pvqw", l3)).max(), np.abs(l2c).max()))
print("   ratio of those two = %.6f" % (np.abs(np.einsum("pvxqwx->pvqw", l3)).max()
                                        / np.abs(l2c).max()))
