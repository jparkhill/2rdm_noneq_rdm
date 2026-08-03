#!/usr/bin/env python3
"""Convex version: is the exact rho-dot in the image of the GKSL cone at rho?

The jump-operator fit reached 8.87% and stopped, but that parameterization is
non-convex, so a floor there could be a stuck optimizer rather than an
obstruction. In the Kossakowski form the problem IS convex:

    L(rho) = -i[K,rho] + sum_ab c_ab ( F_a rho F_b^dag - 1/2 {F_b^dag F_a, rho} )

with K Hermitian and c >= 0. Both terms are linear in the unknowns, so
minimizing || L(rho) - target ||_F over (K, c >= 0) is an SDP. Its optimum is
global: if it is nonzero, the exact instantaneous dynamics genuinely is not
CP-generated at this rho, and that is a real obstruction rather than a numerical
one.

To keep the SDP tractable the jump basis is restricted to a random subspace of
dimension NB, and NB is increased until the residual stops improving.
"""
import sys
import numpy as np
import cvxpy as cp
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from qdiss.symbolic import mkno
from qdiss.dynamics.coherent2 import compile_terms, evaluate
from qdiss.representability.dqg import contract_to_1rdm
from tools.eom2_coherent import ONE, TWO

M, N = 6, 3
dim = 1 << M
PAIRS = [(p, v) for p in range(M) for v in range(p + 1, M)]
NP = len(PAIRS)
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


def rdm(psi, k):
    import itertools
    g = np.zeros((M,) * (2 * k), dtype=complex)
    for t in itertools.product(*([range(M)] * (2 * k))):
        up, lo = t[:k], t[k:]
        v = psi
        for q in lo:
            v = aa(q, v)
        for p in reversed(up):
            v = ad(p, v)
        g[t] = np.vdot(psi, v)
    return g


def to_pair(T):
    A = np.zeros((NP, NP), dtype=complex)
    for a, (p, v) in enumerate(PAIRS):
        for b, (q, w) in enumerate(PAIRS):
            A[a, b] = 2.0 * T[p, v, q, w]
    return A


rng = np.random.default_rng(2026)
A_ = rng.normal(size=(M, M)) + 1j * rng.normal(size=(M, M))
h = 0.25 * (A_ + A_.conj().T)
B = rng.normal(size=(M,) * 4) + 1j * rng.normal(size=(M,) * 4)
V = B - np.transpose(B, (1, 0, 2, 3))
V = V - np.transpose(V, (0, 1, 3, 2))
V = 0.25 * (V + np.transpose(V.conj(), (2, 3, 0, 1)))
psi = np.where(np.array([bin(s).count("1") for s in range(dim)]) == N,
               rng.normal(size=dim) + 1j * rng.normal(size=dim), 0.0)
psi = psi / np.linalg.norm(psi)

g2 = rdm(psi, 2)
g3 = rdm(psi, 3)
g1 = contract_to_1rdm(g2, N)
rho = to_pair(g2)


def gdot(g3t):
    gam = {1: g1, 2: g2, 3: g3t}
    return 1j * (evaluate(ONE_C, h, gam, M, out_rank=4)
                 + 0.25 * evaluate(TWO_C, V, gam, M, out_rank=4))


target = to_pair(gdot(g3))
nt = np.linalg.norm(target)
print("|exact rho-dot| = %.4f" % nt)
print()

for NB in (48, 96, 160):
    r = np.random.default_rng(7)
    basis = [r.normal(size=(NP, NP)) + 1j * r.normal(size=(NP, NP)) for _ in range(NB)]
    # precompute the linear maps  D_ab(rho) = F_a rho F_b^dag - 1/2 {F_b^dag F_a, rho}
    Dm = np.zeros((NB, NB, NP, NP), dtype=complex)
    for a in range(NB):
        for b in range(NB):
            Fa, Fb = basis[a], basis[b]
            Dm[a, b] = (Fa @ rho @ Fb.conj().T
                        - 0.5 * (Fb.conj().T @ Fa @ rho + rho @ Fb.conj().T @ Fa))
    Kv = cp.Variable((NP, NP), hermitian=True)
    C = cp.Variable((NB, NB), hermitian=True)
    diss = sum(C[a, b] * Dm[a, b] for a in range(NB) for b in range(NB))
    L = -1j * (Kv @ rho - rho @ Kv) + diss
    prob = cp.Problem(cp.Minimize(cp.norm(L - target, "fro")), [C >> 0])
    try:
        prob.solve(solver=cp.SCS, verbose=False, max_iters=20000)
        res = prob.value
        print("  NB=%2d  SDP optimum residual %.4f   (%.2f%% of |exact|)   status %s"
              % (NB, res, 100 * res / nt, prob.status))
    except Exception as e:
        print("  NB=%2d  failed: %s" % (NB, str(e)[:60]))
