#!/usr/bin/env python3
"""Does the Q violation survive at M=8, N=3, or was it an M=6 N=3 artifact?

The generality check showed that at M=8 N=3 the Q null space is exactly the
symmetric complement, so Q is full rank on the physical antisymmetric pair space
and its boundary is not active at t=0. But the earlier propagation at M=6 N=3
drove Q to -0.58 over a trajectory. Those are different questions: a state can
start interior and reach the boundary.

So propagate at M=8 N=3 under both closures and watch min D and min Q, exactly as
tests/test_n3_closures.py did at M=6.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from qdiss.symbolic import mkno
from qdiss.dynamics.coherent2 import compile_terms, evaluate, rk4
from qdiss.representability.dqg import contract_to_1rdm, build_Q
from tools.eom2_coherent import ONE, TWO

M, N = 8, 3
dim = 1 << M
ONE_C = compile_terms(ONE, "rs")
TWO_C = compile_terms(TWO, "rtsu")

UP3, LO3 = ("p", "v", "x"), ("q", "w", "y")
sep3 = [(s, f) for s, f in mkno.gamma_to_cumulants(UP3, LO3)
        if not any(n == "L" and r == 3 for n, r, u, l in f)]
SEP3 = [(s, ",".join("".join(u) + "".join(l) for n, r, u, l in f) + "->pvxqwy",
         [(n, r) for n, r, u, l in f]) for s, f in sep3]


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


def sep_gamma3(g1, l2):
    t = {("G", 1): g1, ("L", 2): l2}
    out = np.zeros((M,) * 6, dtype=complex)
    for sign, sub, keys in SEP3:
        out = out + sign * np.einsum(sub, *[t[k] for k in keys], optimize=True)
    return out


def lam2_of(g2, g1):
    return g2 - (np.einsum("pq,vw->pvqw", g1, g1) - np.einsum("pw,vq->pvqw", g1, g1))


rng = np.random.default_rng(2026)
A = rng.normal(size=(M, M)) + 1j * rng.normal(size=(M, M))
h = 0.25 * (A + A.conj().T)
B = rng.normal(size=(M,) * 4) + 1j * rng.normal(size=(M,) * 4)
V = B - np.transpose(B, (1, 0, 2, 3))
V = V - np.transpose(V, (0, 1, 3, 2))
V = 0.25 * (V + np.transpose(V.conj(), (2, 3, 0, 1)))
psi = np.where(np.array([bin(s).count("1") for s in range(dim)]) == N,
               rng.normal(size=dim) + 1j * rng.normal(size=dim), 0.0)
psi = psi / np.linalg.norm(psi)
g2_0 = rdm2(psi)

PAIRD = M * (M - 1) // 2


def mat(T):
    Am = T.reshape(M * M, M * M)
    return 0.5 * (Am + Am.conj().T)


def make(alpha):
    def rhs(g2):
        g1 = contract_to_1rdm(g2, N)
        g3 = (alpha * sep_gamma3(g1, lam2_of(g2, g1)) if alpha
              else np.zeros((M,) * 6, dtype=complex))
        gam = {2: g2, 3: g3}
        return 1j * (evaluate(ONE_C, h, gam, M, out_rank=4)
                     + 0.25 * evaluate(TWO_C, V, gam, M, out_rank=4))
    return rhs


def Hpsi(v):
    out = np.zeros_like(v)
    for r in range(M):
        for s in range(M):
            out += h[r, s] * ad(r, aa(s, v))
    for r in range(M):
        for t in range(M):
            for s in range(M):
                for u in range(M):
                    c = V[r, t, s, u]
                    if c != 0:
                        out += 0.25 * c * ad(r, ad(t, aa(u, aa(s, v))))
    return out


print("M=%d N=%d  pair dim %d  rank g2(0) = %d"
      % (M, N, PAIRD, int((np.linalg.eigvalsh(mat(g2_0)) > 1e-9).sum())))
print()
dt, T = 0.004, 0.8
nstep = int(T / dt)
states = {0.0: g2_0.copy(), 1.0: g2_0.copy()}
pk = psi.copy()
print("%6s | %-24s | %-24s" % ("t", "alpha=0 (gamma_3=0)", "alpha=1 (lambda_3=0)"))
print("%6s | %10s %12s | %10s %12s" % ("", "min D", "min Q", "min D", "min Q"))
for step in range(nstep + 1):
    if step % 50 == 0:
        row = ["%6.2f |" % (step * dt)]
        for a in (0.0, 1.0):
            g = states[a]
            g1 = contract_to_1rdm(g, N)
            mD = float(np.linalg.eigvalsh(mat(g)).min())
            mQ = float(np.linalg.eigvalsh(mat(build_Q(g, g1))).min())
            row.append(" %10.3e %12.3e |" % (mD, mQ))
        print("".join(row))
    for a in (0.0, 1.0):
        states[a] = rk4(make(a), states[a], dt)
    pk = rk4(lambda v: -1j * Hpsi(v), pk, dt)
