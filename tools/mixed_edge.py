#!/usr/bin/env python3
"""Mixed states, and sampling the PSD edge exactly.

Two concerns, both John's.

1. Every measurement so far used a PURE state. They were multi-determinant
   (random full-CI vectors), not single determinants, but pure. That matters
   because X_KK = 0 -- the rank conservation that forced the CP construction --
   may be a pure-state artifact. Here gamma_2 and gamma_3 are built from a mixed
   ensemble rho = sum p_i |psi_i><psi_i|.

2. The question is not whether a non-PSD density gets mapped back, it is whether
   a density ON the boundary escapes. So sample the edge exactly: for any rho >= 0
   and unit v, rho' = P rho P with P = I - vv^dag satisfies rho' v = 0 exactly and
   stays PSD. Rescale the trace and evaluate the flux <v|gamma_2-dot|v> there
   under each closure. No trajectory, no integrator.
"""
import itertools
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from qdiss.symbolic import mkno
from qdiss.dynamics.coherent2 import compile_terms, evaluate
from qdiss.representability.dqg import contract_to_1rdm, build_Q, build_G
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


def apply_str(up, lo, v):
    for q in lo:
        v = aa(q, v)
    for p in reversed(up):
        v = ad(p, v)
    return v


def rdm_mixed(states, probs, k):
    g = np.zeros((M,) * (2 * k), dtype=complex)
    for t in itertools.product(*([range(M)] * (2 * k))):
        up, lo = t[:k], t[k:]
        acc = 0.0 + 0j
        for p_, psi in zip(probs, states):
            acc += p_ * np.vdot(psi, apply_str(up, lo, psi))
        g[t] = acc
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
            T[p, v, q, w] = x; T[v, p, q, w] = -x
            T[p, v, w, q] = -x; T[v, p, w, q] = x
    return T


rng = np.random.default_rng(2026)
A_ = rng.normal(size=(M, M)) + 1j * rng.normal(size=(M, M))
h = 0.25 * (A_ + A_.conj().T)
B = rng.normal(size=(M,) * 4) + 1j * rng.normal(size=(M,) * 4)
V = B - np.transpose(B, (1, 0, 2, 3))
V = V - np.transpose(V, (0, 1, 3, 2))
V = 0.25 * (V + np.transpose(V.conj(), (2, 3, 0, 1)))
popc = np.array([bin(s).count("1") for s in range(dim)])


def randstate(r):
    v = np.where(popc == N, r.normal(size=dim) + 1j * r.normal(size=dim), 0.0)
    return v / np.linalg.norm(v)


def gdot(g1, g2, g3):
    gam = {1: g1, 2: g2, 3: g3}
    return 1j * (evaluate(ONE_C, h, gam, M, out_rank=4)
                 + 0.25 * evaluate(TWO_C, V, gam, M, out_rank=4))


def sep_g3(g1, g2):
    l2 = g2 - (np.einsum("pq,vw->pvqw", g1, g1) - np.einsum("pw,vq->pvqw", g1, g1))
    t = {("G", 1): g1, ("L", 2): l2}
    out = np.zeros((M,) * 6, dtype=complex)
    for sign, sub, keys in SEP3:
        out = out + sign * np.einsum(sub, *[t[k] for k in keys], optimize=True)
    return out


print("=== part 1: does X_KK = 0 survive for MIXED states? ===")
print("%14s %8s %10s %14s %14s" % ("ensemble", "rank g2", "dim ker", "|X_KK|", "|X_KR|"))
for nc in (1, 2, 3, 5, 10):
    r = np.random.default_rng(4000 + nc)
    states = [randstate(r) for _ in range(nc)]
    p = r.random(nc); p /= p.sum()
    g2 = rdm_mixed(states, p, 2)
    g3 = rdm_mixed(states, p, 3)
    g1 = contract_to_1rdm(g2, N)
    X = to_pair(gdot(g1, g2, g3))
    X = 0.5 * (X + X.conj().T)
    rho = to_pair(g2); rho = 0.5 * (rho + rho.conj().T)
    w, U = np.linalg.eigh(rho)
    Kb = U[:, w < 1e-9]; Rb = U[:, w >= 1e-9]
    xkk = abs(Kb.conj().T @ X @ Kb).max() if Kb.shape[1] else 0.0
    xkr = abs(Kb.conj().T @ X @ Rb).max() if Kb.shape[1] else 0.0
    print("%14s %8d %10d %14.3e %14.3e"
          % ("pure" if nc == 1 else "mixed(%d)" % nc,
             int((w > 1e-9).sum()), Kb.shape[1], xkk, xkr))

print()
print("=== part 2: exact edge sampling -- does the flow escape? ===")
print("rho' = P rho P with P = I - vv^dag, so rho'v = 0 exactly; trace rescaled.")
print("flux is <v|gamma_2-dot|v> at that boundary point. negative = ESCAPES.")
print()
r = np.random.default_rng(77)
states = [randstate(r) for _ in range(4)]
p = r.random(4); p /= p.sum()
g2_m = rdm_mixed(states, p, 2)
g1_m = contract_to_1rdm(g2_m, N)
rho_m = 0.5 * (to_pair(g2_m) + to_pair(g2_m).conj().T)
print("base mixed state: rank gamma_2 = %d of %d"
      % (int((np.linalg.eigvalsh(rho_m) > 1e-9).sum()), NP))
print()
print("%5s | %12s %12s | %10s %10s" % ("trial", "gamma_3=0", "lambda_3=0", "minQ", "minG"))
bad0 = bad1 = 0
ntr = 20
for t in range(ntr):
    rr = np.random.default_rng(9100 + t)
    v = rr.normal(size=NP) + 1j * rr.normal(size=NP)
    v /= np.linalg.norm(v)
    P = np.eye(NP, dtype=complex) - np.outer(v, v.conj())
    rp = P @ rho_m @ P
    rp = rp * (N * (N - 1) / np.trace(rp).real)
    g2p = from_pair(rp)
    g1p = contract_to_1rdm(g2p, N)
    ev = np.linalg.eigvalsh(0.5 * (rp + rp.conj().T))
    fl = {}
    for lbl, g3p in (("g0", np.zeros((M,) * 6, dtype=complex)),
                     ("g1", sep_g3(g1p, g2p))):
        X = to_pair(gdot(g1p, g2p, g3p))
        X = 0.5 * (X + X.conj().T)
        fl[lbl] = float(np.real(v.conj() @ X @ v))
    def mn(T):
        A = T.reshape(M * M, M * M)
        return float(np.linalg.eigvalsh(0.5 * (A + A.conj().T)).min())
    mq, mg = mn(build_Q(g2p, g1p)), mn(build_G(g2p, g1p))
    if fl["g0"] < -1e-9: bad0 += 1
    if fl["g1"] < -1e-9: bad1 += 1
    print("%5d | %12.4e %12.4e | %10.2e %10.2e" % (t, fl["g0"], fl["g1"], mq, mg))

print()
print("escapes at the edge:  gamma_3=0 %d/%d    lambda_3=0 %d/%d"
      % (bad0, ntr, bad1, ntr))
