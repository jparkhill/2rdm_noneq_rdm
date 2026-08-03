#!/usr/bin/env python3
"""How well can a CP generator on pair space reproduce the EXACT gamma_2-dot?

gamma_3 = 0 is positivity-safe but isospectral, so it cannot relax; lambda_3 = 0
relaxes but leaves the cone. The question is whether a completely positive
generator can do both -- match the true dynamics while being CP by construction.

Fit  L(rho) = -i[K,rho] + sum_k ( L_k rho L_k^dag - 1/2 {L_k^dag L_k, rho} )
to the exact rho-dot at a fixed rho, over K Hermitian and free jump operators.
Parameterizing by jump operators keeps the Kossakowski matrix PSD automatically,
so every candidate is genuinely CP.

Baselines to beat: the residual from gamma_3 = 0 and from lambda_3 = 0.

A necessary condition first: for ANY GKSL generator the commutator and
anticommutator parts vanish on ker rho, leaving sum_k <v|L_k rho L_k^dag|v> >= 0.
So D flux >= 0 is necessary for CP-representability, and the exact dynamics must
satisfy it -- it does, at ~1e-16, i.e. marginally. That marginality forces
L_k^dag to map ker rho into ker rho, which is a real structural constraint on the
jump operators rather than a free fit.
"""
import sys
import numpy as np
from scipy.optimize import minimize
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
NJ = 12
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


def sep_gamma3():
    l2 = g2 - (np.einsum("pq,vw->pvqw", g1, g1) - np.einsum("pw,vq->pvqw", g1, g1))
    t = {("G", 1): g1, ("L", 2): l2}
    out = np.zeros((M,) * 6, dtype=complex)
    for sign, sub, keys in SEP3:
        out = out + sign * np.einsum(sub, *[t[k] for k in keys], optimize=True)
    return out


target = to_pair(gdot(g3))
base0 = to_pair(gdot(np.zeros_like(g3)))
base1 = to_pair(gdot(sep_gamma3()))
nt = np.linalg.norm(target)
print("M=%d N=%d  pair dim %d  |exact rho-dot| = %.4f" % (M, N, NP, nt))
print("  gamma_3 = 0   residual %.4f   (%.1f%% of |exact|)"
      % (np.linalg.norm(base0 - target), 100 * np.linalg.norm(base0 - target) / nt))
print("  lambda_3 = 0  residual %.4f   (%.1f%%)"
      % (np.linalg.norm(base1 - target), 100 * np.linalg.norm(base1 - target) / nt))

w, U = np.linalg.eigh(0.5 * (rho + rho.conj().T))
nullr = U[:, w < 1e-9]
dflux = np.real(np.einsum("ai,ab,bi->i", nullr.conj(),
                          0.5 * (target + target.conj().T), nullr))
print("  D flux of the exact dynamics on ker rho: min %.3e max %.3e"
      % (dflux.min(), dflux.max()))
print()


def unpack(x):
    i = 0
    kr = x[i:i + NP * NP].reshape(NP, NP); i += NP * NP
    ki = x[i:i + NP * NP].reshape(NP, NP); i += NP * NP
    K = kr + kr.T + 1j * (ki - ki.T)
    Ls = []
    for _ in range(NJ):
        lr = x[i:i + NP * NP].reshape(NP, NP); i += NP * NP
        li = x[i:i + NP * NP].reshape(NP, NP); i += NP * NP
        Ls.append(lr + 1j * li)
    return K, Ls


def gen(K, Ls):
    d = -1j * (K @ rho - rho @ K)
    for L in Ls:
        LdL = L.conj().T @ L
        d = d + L @ rho @ L.conj().T - 0.5 * (LdL @ rho + rho @ LdL)
    return d


def obj(x):
    K, Ls = unpack(x)
    r = gen(K, Ls) - target
    return float(np.real(np.vdot(r, r)))


npar = 2 * NP * NP * (1 + NJ)
print("fitting %d real parameters (K + %d jump operators) to %d complex entries"
      % (npar, NJ, NP * NP))
best = None
for trial in range(3):
    x0 = np.random.default_rng(100 + trial).normal(size=npar) * 0.1
    res = minimize(obj, x0, method="L-BFGS-B",
                   options={"maxiter": 4000, "maxfun": 60000})
    if best is None or res.fun < best.fun:
        best = res
    print("  trial %d: residual %.6e" % (trial, np.sqrt(max(res.fun, 0))))

K, Ls = unpack(best.x)
fit = gen(K, Ls)
rel = np.linalg.norm(fit - target) / nt
print()
print("best CP fit residual %.4f  (%.2f%% of |exact|)"
      % (np.linalg.norm(fit - target), 100 * rel))
print("  d(tr rho^2) exact %.4e   fit %.4e"
      % (2 * np.trace(rho @ target).real, 2 * np.trace(rho @ fit).real))
print("  trace preservation: exact %.2e  fit %.2e"
      % (abs(np.trace(target).real), abs(np.trace(fit).real)))
