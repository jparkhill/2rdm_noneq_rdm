#!/usr/bin/env python3
"""The tangency algebra at ranks 2 AND 3, so nothing is generalised from lambda_2 alone.

SETUP. At rank k the propagated object is gamma_k, whose independent content is a matrix indexed
by unordered k-tuples of spin-orbitals. Write COMB_k for the C(M,k) sorted tuples, and define

    Mat_k(T)[I,J] = k! * T[I + J]                       for I, J in COMB_k
    Ten_k(Z)[sigma(I) + tau(J)] = sign(sigma) sign(tau) Z[I,J] / k!

so Mat_k(Ten_k(Z)) = Z. For k = 2 this is the pair matrix D[I,J] = 2 gamma_2[p,v,q,w] used
earlier. Mat_k(gamma_k) is a Gram matrix of the k-fold annihilated state, hence positive
semidefinite, and the relevant cone at every rank is the PSD cone of that matrix.

THE CONDITION. At a boundary point with kernel projector P_K the tangent cone to the PSD cone is
{Y : P_K Y P_K >= 0}, so Nagumo's no-outward-flow requirement is P_K (d/dt Mat_k) P_K >= 0. Two
forms give it unconditionally, for ANY A and L_i:

    A X + X A^dag        compresses to exactly 0 (X P_K = P_K X = 0)   -- tangent
    sum_i L_i X L_i^dag  compresses to sum_i (P_K L_i) X (P_K L_i)^dag -- inward, since X >= 0

THE TEST. A Lyapunov map A X + X A^dag sends a rank-1 input |Vec><Vec| to
(A Vec)<Vec| + |Vec>(A Vec)^dag, of rank at most 2. So feed rank-1 matrices and measure the
output rank; anything above 2 falsifies the Lyapunov form. This is pure algebra on the term list
and does not depend on the state being N-representable, so a random Hermitian gamma_1 and a random
rank-1 Mat_k are legitimate probes.

The rank-k EOM couples gamma_k to gamma_k and to gamma_{k+1} and to nothing lower, so the term
list partitions cleanly: the gamma_k terms are the drift under test, and the gamma_{k+1} terms are
exactly what a reconstruction closure has to supply. Term lists come from qdiss.eom.wick.terms(k),
bare moments, validated against qdiss.eom.fullci at 8.7e-16 at rank 3 and 1e-15 at ranks 1-2.
"""
import itertools
import sys
import numpy as np
sys.path.insert(0, ".")
from qdiss.eom import wick
from qdiss.eom.hamiltonian import random_h, random_V

M = 6
rng = np.random.default_rng(11)
h, V = random_h(M, rng), random_V(M, rng)
EYE = np.eye(M, dtype=complex)
g1r = rng.normal(size=(M, M)) + 1j * rng.normal(size=(M, M))
g1r = 0.5 * (g1r + g1r.conj().T)


def perm_sign(p):
    s, p = 1, list(p)
    for i in range(len(p)):
        for j in range(len(p) - 1 - i):
            if p[j] > p[j + 1]:
                p[j], p[j + 1] = p[j + 1], p[j]; s = -s
    return s


def combs(k):
    return list(itertools.combinations(range(M), k))


def ten_from_mat(Z, k):
    C = combs(k)
    fac = float(np.math.factorial(k))
    T = np.zeros((M,) * (2 * k), dtype=complex)
    for i, I in enumerate(C):
        for j, J in enumerate(C):
            z = Z[i, j] / fac
            for su in itertools.permutations(range(k)):
                pu = perm_sign(su); Iu = tuple(I[a] for a in su)
                for sl in itertools.permutations(range(k)):
                    T[Iu + tuple(J[b] for b in sl)] = pu * perm_sign(sl) * z
    return T


def mat_from_ten(T, k):
    C = combs(k)
    fac = float(np.math.factorial(k))
    return np.array([[fac * T[I + J] for J in C] for I in C])


def partition(k):
    one, two = wick.terms(k)
    lo, hi = ([], []), ([], [])
    for which, lst in ((0, one), (1, two)):
        for s, f in lst:
            keys = [(n, r) for n, r, u, l in f]
            subs = ["".join(u) + "".join(l) for n, r, u, l in f]
            (hi if ("G", k + 1) in keys else lo)[which].append((s, subs, keys))
    return lo, hi


def apply_terms(terms, k, tab):
    out = "".join(wick.OPERATOR[k][0]) + "".join(wick.OPERATOR[k][1])
    tot = np.zeros((M,) * (2 * k), dtype=complex)
    for s, subs, keys in terms[0]:
        tot += s * np.einsum(",".join(["rs"] + subs) + "->" + out, h,
                             *[tab[key] for key in keys], optimize=True)
    for s, subs, keys in terms[1]:
        tot += s * 0.25 * np.einsum(",".join(["rtsu"] + subs) + "->" + out, V,
                                    *[tab[key] for key in keys], optimize=True)
    return 1j * tot


print("tangency algebra at ranks 2 and 3;  M = %d spin-orbitals" % M)
print("EOM term lists from qdiss.eom.wick, validated against fullci (8.7e-16 at rank 3)")
print()
print("   %-6s %8s | %-28s | %-24s" % ("rank", "dim", "terms carrying gamma_k", "carrying gamma_{k+1}"))
for k in (2, 3):
    lo, hi = partition(k)
    print("   %-6d %8d | %2d one-body + %2d two-body      | %2d one-body + %2d two-body"
          % (k, len(combs(k)), len(lo[0]), len(lo[1]), len(hi[0]), len(hi[1])))
print()
print("rank-1 probe of the gamma_k-dependent drift: Lyapunov form requires output rank <= 2")
print("   %-6s %8s %14s %20s" % ("rank", "trial", "output rank", "sv[2:5]/sv[0]"))
for k in (2, 3):
    n = len(combs(k))
    lo, _ = partition(k)
    for t in range(3):
        vv = rng.normal(size=n) + 1j * rng.normal(size=n)
        vv /= np.linalg.norm(vv)
        Z = np.outer(vv, vv.conj())
        tab = {("D", 1): EYE, ("G", 1): g1r, ("G", k): ten_from_mat(Z, k)}
        out = mat_from_ten(apply_terms(lo, k, tab), k)
        sv = np.linalg.svd(out, compute_uv=False)
        print("   %-6d %8d %14d %20s"
              % (k, t, int((sv > 1e-10 * sv.max()).sum()),
                 np.array2string(sv[2:5] / sv.max(), precision=2)))
