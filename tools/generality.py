#!/usr/bin/env python3
"""Does the Q-tangency structure generalize beyond M=6, N=3?

Everything measured so far sits at one filling and one state, and the (n, 2-n)
geminal pairing already showed how that goes wrong. Before building an SDP on the
tangency condition, check at other (M, N):

  * dim ker Q, and how many flux directions are structurally zero
  * what fraction of random GKSL generators are Q-inward
  * whether the Q-inward ones relax (non-isospectral)

If the structure is stable across fillings, the tangency condition is load-bearing
and the feasibility solve is worth setting up.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from qdiss.representability.dqg import contract_to_1rdm, build_Q


def run(M, N, nrun=24, seed=2026):
    dim = 1 << M
    eye = np.eye(M, dtype=complex)
    PAIRS = [(p, v) for p in range(M) for v in range(p + 1, M)]
    NP = len(PAIRS)

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

    rng = np.random.default_rng(seed)
    popc = np.array([bin(s).count("1") for s in range(dim)])
    psi = np.where(popc == N, rng.normal(size=dim) + 1j * rng.normal(size=dim), 0.0)
    psi = psi / np.linalg.norm(psi)

    g2 = np.zeros((M,) * 4, dtype=complex)
    for p in range(M):
        for v in range(M):
            for q in range(M):
                for w in range(M):
                    g2[p, v, q, w] = np.vdot(psi, ad(p, ad(v, aa(w, aa(q, psi)))))
    g1 = contract_to_1rdm(g2, N)

    def to_pair(T):
        A = np.zeros((NP, NP), dtype=complex)
        for a, (p, v) in enumerate(PAIRS):
            for b, (q, w) in enumerate(PAIRS):
                A[a, b] = 2.0 * T[p, v, q, w]
        return A

    def from_pair(A):
        T = np.zeros((M,) * 4, dtype=complex)
        for a, (p, v) in enumerate(PAIRS):
            for b, (q, w) in enumerate(PAIRS):
                x = 0.5 * A[a, b]
                T[p, v, q, w] = x; T[v, p, q, w] = -x
                T[p, v, w, q] = -x; T[v, p, w, q] = x
        return T

    def mat(T):
        A = T.reshape(M * M, M * M)
        return 0.5 * (A + A.conj().T)

    def qlin(dg2):
        dg1 = contract_to_1rdm(dg2, N)
        return (- np.einsum("ki,lj->ijkl", eye, dg1)
                + np.einsum("kj,li->ijkl", eye, dg1)
                + np.einsum("li,kj->ijkl", eye, dg1)
                - np.einsum("lj,ki->ijkl", eye, dg1)
                + np.einsum("klij->ijkl", dg2))

    rho = to_pair(g2)
    wq, Uq = np.linalg.eigh(mat(build_Q(g2, g1)))
    nullQ = Uq[:, wq < 1e-9]
    wd = np.linalg.eigvalsh(mat(g2))

    nq = 0
    nz_counts = []
    relax = []
    for r_ in range(nrun):
        r = np.random.default_rng(9000 + r_)
        A = r.normal(size=(NP, NP)) + 1j * r.normal(size=(NP, NP))
        K = 0.5 * (A + A.conj().T)
        Ls = [0.3 * (r.normal(size=(NP, NP)) + 1j * r.normal(size=(NP, NP)))
              for _ in range(1 + r_ % 3)]
        d = -1j * (K @ rho - rho @ K)
        for L in Ls:
            d = d + L @ rho @ L.conj().T - 0.5 * (L.conj().T @ L @ rho
                                                  + rho @ L.conj().T @ L)
        dg2 = from_pair(d)
        vals = np.real(np.einsum("ai,ab,bi->i", nullQ.conj(), mat(qlin(dg2)), nullQ))
        nz_counts.append(int((np.abs(vals) < 1e-12).sum()))
        if vals.min() > -1e-9:
            nq += 1
            relax.append(abs(2 * np.trace(rho @ d).real))

    print("M=%d N=%d  half-filled=%-5s | pair dim %2d | rank g2 %2d | dim ker Q %2d"
          % (M, N, N * 2 == M, NP, int((wd > 1e-9).sum()), nullQ.shape[1]))
    print("    structurally-zero flux directions: %s"
          % (("%d (constant)" % nz_counts[0]) if len(set(nz_counts)) == 1
             else "varies %s" % sorted(set(nz_counts))))
    print("    Q inward: %d/%d" % (nq, nrun))
    if relax:
        print("    of those, max relaxation |d tr rho^2| = %.3e ; isospectral: %d"
              % (max(relax), sum(1 for x in relax if x < 1e-9)))
    print()


for M, N in ((6, 3), (6, 2), (8, 3), (8, 4)):
    run(M, N)
