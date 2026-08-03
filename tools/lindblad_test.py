#!/usr/bin/env python3
"""Is complete positivity on the pair space enough to keep Q and G?

Construction rather than search. Build the generator in manifest Lindblad form on
the antisymmetric two-particle space,

    rho-dot = -i[K, rho] + sum_k ( L_k rho L_k^dag - 1/2 {L_k^dag L_k, rho} )

which preserves rho >= 0 by construction and, unlike gamma_3 = 0, is NOT
isospectral -- so it can relax the geminal occupations. Then measure the boundary
flux in all three sectors on an exactly representable gamma_2.

D must come out inward: that is what CP means. The question is Q and G, since
they are DIFFERENT affine functions of gamma_2 and gamma_1, and nothing about a
CP map on gamma_2 obviously controls them.

If Q and G come out inward too, the CP route is the answer and the search for a
closure is over. If they do not, complete positivity on the pair space is
necessary but not sufficient, and the extra constraint has to be identified.
"""
import itertools
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from qdiss.representability.dqg import contract_to_1rdm, build_Q, build_G

M, N = 6, 3
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


def rdm2(psi):
    g = np.zeros((M,) * 4, dtype=complex)
    for p in range(M):
        for v in range(M):
            for q in range(M):
                for w in range(M):
                    g[p, v, q, w] = np.vdot(psi, ad(p, ad(v, aa(w, aa(q, psi)))))
    return g


def to_pair(T):
    """tensor -> matrix on the antisymmetric pair basis (orthonormal)."""
    A = np.zeros((NP, NP), dtype=complex)
    for a, (p, v) in enumerate(PAIRS):
        for b, (q, w) in enumerate(PAIRS):
            A[a, b] = 2.0 * T[p, v, q, w]
    return A


def from_pair(A):
    """matrix on the antisymmetric pair basis -> antisymmetric tensor."""
    T = np.zeros((M,) * 4, dtype=complex)
    for a, (p, v) in enumerate(PAIRS):
        for b, (q, w) in enumerate(PAIRS):
            x = 0.5 * A[a, b]
            T[p, v, q, w] = x
            T[v, p, q, w] = -x
            T[p, v, w, q] = -x
            T[v, p, w, q] = x
    return T


# ---------- an exactly representable gamma_2 ----------
rng = np.random.default_rng(2026)
psi = np.where(np.array([bin(s).count("1") for s in range(dim)]) == N,
               rng.normal(size=dim) + 1j * rng.normal(size=dim), 0.0)
psi = psi / np.linalg.norm(psi)
g2 = rdm2(psi)
g1 = contract_to_1rdm(g2, N)

rt = from_pair(to_pair(g2))
print("round trip tensor -> pair -> tensor : %.3e" % np.abs(rt - g2).max())
rho = to_pair(g2)
print("pair-space rho: %dx%d, tr = %.6f (should be N(N-1) = %d)"
      % (NP, NP, np.trace(rho).real, N * (N - 1)))
w = np.linalg.eigvalsh(0.5 * (rho + rho.conj().T))
print("rho spectrum min %.3e  max %.6f  rank %d"
      % (w.min(), w.max(), int((w > 1e-9).sum())))
print()


def qdot_from(dg2):
    dg1 = contract_to_1rdm(dg2, N)
    return (- np.einsum("ki,lj->ijkl", eye, dg1)
            + np.einsum("kj,li->ijkl", eye, dg1)
            + np.einsum("li,kj->ijkl", eye, dg1)
            - np.einsum("lj,ki->ijkl", eye, dg1)
            + np.einsum("klij->ijkl", dg2))


def gdot_from(dg2):
    dg1 = contract_to_1rdm(dg2, N)
    return np.einsum("jl,ik->ijkl", eye, dg1) - dg2.transpose(0, 3, 2, 1)


def mat(T):
    A = T.reshape(M * M, M * M)
    return 0.5 * (A + A.conj().T)


def fluxes(dg2, label):
    out = {}
    for name, T, dot in (("D", g2, dg2),
                         ("Q", build_Q(g2, g1), qdot_from(dg2)),
                         ("G", build_G(g2, g1), gdot_from(dg2))):
        wv, U = np.linalg.eigh(mat(T))
        null = U[:, wv < 1e-9]
        vals = np.real(np.einsum("ai,ab,bi->i", null.conj(), mat(dot), null))
        out[name] = vals.min()
    print("  %-26s D %+11.4e   Q %+11.4e   G %+11.4e"
          % (label, out["D"], out["Q"], out["G"]))
    return out


def lindblad(rho, K, Ls):
    d = -1j * (K @ rho - rho @ K)
    for L in Ls:
        d = d + L @ rho @ L.conj().T - 0.5 * (L.conj().T @ L @ rho
                                              + rho @ L.conj().T @ L)
    return d


print("boundary flux, min over null directions (>= 0 is inward):")
print()

nrun = 8
tally = {"D": 0, "Q": 0, "G": 0}
for run in range(nrun):
    r = np.random.default_rng(500 + run)
    A = r.normal(size=(NP, NP)) + 1j * r.normal(size=(NP, NP))
    K = 0.5 * (A + A.conj().T)
    Ls = []
    for _ in range(3):
        B = r.normal(size=(NP, NP)) + 1j * r.normal(size=(NP, NP))
        Ls.append(0.3 * B)
    drho = lindblad(rho, K, Ls)
    dg2 = from_pair(drho)
    f = fluxes(dg2, "Lindblad run %d" % run)
    for k in tally:
        if f[k] > -1e-9:
            tally[k] += 1

print()
print("inward in %d/%d runs:  D %d   Q %d   G %d"
      % (nrun, nrun, tally["D"], tally["Q"], tally["G"]))
print()
print("is the generator isospectral? (it must not be, or it cannot relax)")
drho = lindblad(rho, np.zeros((NP, NP), dtype=complex),
                [0.3 * (np.random.default_rng(9).normal(size=(NP, NP))
                        + 1j * np.random.default_rng(10).normal(size=(NP, NP)))])
print("  d(tr rho)   = %+.3e   (Lindblad preserves trace)" % np.trace(drho).real)
print("  d(tr rho^2) = %+.3e   (nonzero => spectrum moves => relaxation)"
      % (2 * np.trace(rho @ drho).real))
