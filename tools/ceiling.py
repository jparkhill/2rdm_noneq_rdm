#!/usr/bin/env python3
"""What is the ceiling on a geminal occupation?

Two candidate answers. Either 2, the Pauli intuition that a geminal is a
two-electron function; or O(N), the Yang/Sasaki bound saturated by a condensate,
which in the tr gamma_2 = N(N-1) normalization used here would allow up to about N.

A random state gave 1.853 at N=3, consistent with both. A deliberately
condensate-like state -- a geminal power |G^dag)^{N/2}|vac> with all pairs in phase
-- maximizes the top eigenvalue, so if the ceiling exceeds 2 that is where it shows.

Also reports the largest occupation over random pure and mixed samples, for
comparison, and the trace so the normalization is never in doubt.
"""
import itertools
import sys
import numpy as np
sys.path.insert(0, ".")

M_LIST = [(4, 2), (6, 2), (6, 3), (8, 4), (10, 4), (10, 5)]


def build(M):
    dim = 1 << M

    def ad(p, v):
        out = np.zeros_like(v); bit = 1 << p
        src = np.nonzero(v)[0]
        for s in src:
            if s & bit: continue
            out[s | bit] += ((-1.0) ** bin(s & (bit - 1)).count("1")) * v[s]
        return out

    def aa(p, v):
        out = np.zeros_like(v); bit = 1 << p
        for s in np.nonzero(v)[0]:
            if not (s & bit): continue
            out[s ^ bit] += ((-1.0) ** bin(s & (bit - 1)).count("1")) * v[s]
        return out
    return dim, ad, aa


def geminal_occ(M, psi, ad, aa):
    """eigenvalues of gamma_2 on the antisymmetric pair space, tr = N(N-1)."""
    PAIRS = [(p, v) for p in range(M) for v in range(p + 1, M)]
    NP = len(PAIRS)
    D = np.zeros((NP, NP), dtype=complex)
    kets = {}
    for b, (q, w) in enumerate(PAIRS):
        kets[b] = aa(w, aa(q, psi))
    for a, (p, v) in enumerate(PAIRS):
        left = aa(v, aa(p, psi))
        for b in range(NP):
            D[a, b] = 2.0 * np.vdot(left, kets[b])
    D = 0.5 * (D + D.conj().T)
    return np.linalg.eigvalsh(D)


print("%8s %10s %12s %12s %12s" % ("(M,N)", "state", "tr gamma_2", "max n_g", "N(N-1)"))
for M, N in M_LIST:
    dim, ad, aa = build(M)
    popc = np.array([bin(s).count("1") for s in range(dim)])

    # random pure
    r = np.random.default_rng(11)
    best_rand = 0.0
    for _ in range(4):
        psi = np.where(popc == N, r.normal(size=dim) + 1j * r.normal(size=dim), 0.0)
        psi = psi / np.linalg.norm(psi)
        w = geminal_occ(M, psi, ad, aa)
        best_rand = max(best_rand, w.max())
        tr = w.sum()
    print("%8s %10s %12.4f %12.6f %12d" % ("(%d,%d)" % (M, N), "random", tr,
                                           best_rand, N * (N - 1)))

    # condensate: geminal power with all pairs in phase, (sum_k a+_{2k} a+_{2k+1})^{N/2}
    if N % 2 == 0:
        vac = np.zeros(dim, dtype=complex); vac[0] = 1.0
        npair = M // 2

        def Gdag(v):
            out = np.zeros_like(v)
            for k in range(npair):
                out += ad(2 * k, ad(2 * k + 1, v))
            return out
        psi = vac
        for _ in range(N // 2):
            psi = Gdag(psi)
        nrm = np.linalg.norm(psi)
        if nrm > 1e-12:
            psi = psi / nrm
            w = geminal_occ(M, psi, ad, aa)
            print("%8s %10s %12.4f %12.6f %12d"
                  % ("", "condensate", w.sum(), w.max(), N * (N - 1)))
