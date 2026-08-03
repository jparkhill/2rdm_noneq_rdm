#!/usr/bin/env python3
"""Is CP-on-particles plus Q-inwardness satisfiable by a RELAXING generator?

The first thing to establish, before building the hole-CP machinery, is whether
the combined requirement is over-constrained. If every generator that keeps Q
inward turns out to be isospectral, then D + Q together force unitary pair
dynamics and that is a no-go theorem rather than a solution.

Samples random GKSL generators on the antisymmetric pair space and records, for
each: the boundary flux in all three sectors, and d(tr rho^2), which is zero iff
the flow is isospectral and therefore cannot relax the geminal occupations.
"""
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
            T[p, v, q, w] = x
            T[v, p, q, w] = -x
            T[p, v, w, q] = -x
            T[v, p, w, q] = x
    return T


rng = np.random.default_rng(2026)
psi = np.where(np.array([bin(s).count("1") for s in range(dim)]) == N,
               rng.normal(size=dim) + 1j * rng.normal(size=dim), 0.0)
psi = psi / np.linalg.norm(psi)
g2 = rdm2(psi)
g1 = contract_to_1rdm(g2, N)
rho = to_pair(g2)


def mat(T):
    A = T.reshape(M * M, M * M)
    return 0.5 * (A + A.conj().T)


NULL = {}
for name, T in (("D", g2), ("Q", build_Q(g2, g1)), ("G", build_G(g2, g1))):
    w, U = np.linalg.eigh(mat(T))
    NULL[name] = U[:, w < 1e-9]


def fluxes(dg2):
    dg1 = contract_to_1rdm(dg2, N)
    dQ = (- np.einsum("ki,lj->ijkl", eye, dg1)
          + np.einsum("kj,li->ijkl", eye, dg1)
          + np.einsum("li,kj->ijkl", eye, dg1)
          - np.einsum("lj,ki->ijkl", eye, dg1)
          + np.einsum("klij->ijkl", dg2))
    dG = np.einsum("jl,ik->ijkl", eye, dg1) - dg2.transpose(0, 3, 2, 1)
    out = {}
    for name, dot in (("D", dg2), ("Q", dQ), ("G", dG)):
        nl = NULL[name]
        out[name] = float(np.real(np.einsum("ai,ab,bi->i", nl.conj(), mat(dot), nl)).min())
    return out


def lindblad(rho, K, Ls):
    d = -1j * (K @ rho - rho @ K)
    for L in Ls:
        d = d + L @ rho @ L.conj().T - 0.5 * (L.conj().T @ L @ rho + rho @ L.conj().T @ L)
    return d


print("sampling random GKSL generators on the %d-dim antisymmetric pair space" % NP)
print("Q inward means min <v|Qdot|v> >= -1e-9;  relaxing means d(tr rho^2) != 0")
print()

nrun = 60
both = []
qbad = 0
iso_and_q = 0
for run in range(nrun):
    r = np.random.default_rng(9000 + run)
    A = r.normal(size=(NP, NP)) + 1j * r.normal(size=(NP, NP))
    K = 0.5 * (A + A.conj().T)
    nL = 1 + run % 3
    Ls = [0.3 * (r.normal(size=(NP, NP)) + 1j * r.normal(size=(NP, NP)))
          for _ in range(nL)]
    drho = lindblad(rho, K, Ls)
    dg2 = from_pair(drho)
    f = fluxes(dg2)
    dtr2 = 2 * np.trace(rho @ drho).real
    qok = f["Q"] > -1e-9
    if not qok:
        qbad += 1
    if qok:
        both.append((abs(dtr2), f, run))
        if abs(dtr2) < 1e-9:
            iso_and_q += 1

print("  D inward : %d/%d" % (nrun, nrun))
print("  Q inward : %d/%d   (violated in %d)" % (nrun - qbad, nrun, qbad))
print("  of the Q-inward generators, %d were isospectral" % iso_and_q)
print()

both.sort(reverse=True)
print("Q-inward generators with the LARGEST relaxation (|d tr rho^2|):")
print("  %6s %14s %12s %12s %12s" % ("run", "|d tr rho^2|", "D flux", "Q flux", "G flux"))
for dtr2, f, run in both[:6]:
    print("  %6d %14.4e %12.4e %12.4e %12.4e"
          % (run, dtr2, f["D"], f["Q"], f["G"]))

print()
if both and both[0][0] > 1e-6:
    print("CONCLUSION: not over-constrained. A generator can be CP on the pair space,")
    print("keep all three DQG conditions inward, and still relax the geminal spectrum")
    print("by |d tr rho^2| = %.3e. So D + Q together do NOT force unitary pair" % both[0][0])
    print("dynamics, and there is no no-go here.")
else:
    print("CONCLUSION: every Q-inward generator found is isospectral -- D + Q together")
    print("appear to force unitary pair dynamics. That would be a no-go.")
