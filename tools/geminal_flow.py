#!/usr/bin/env python3
"""Sources and drains for a single geminal: are they properly blocked?

At M=6 N=3 the exact gamma_2 has six fractionally occupied geminals. For each,
n_g-dot = <g|gamma_2-dot|g>. The question is whether the terms feeding a geminal
carry a factor that dies as it fills, and the terms draining it a factor that dies
as it empties -- the rank-2 analogue of n_q (1-n_p) versus n_p (1-n_q).

Two probes.

A. Decompose n_g-dot at the true state into positive contributions (sources) and
   negative ones (drains), by term, so the balance is visible.

B. Sweep the occupation. Replace n_g by s while holding the other geminals fixed,
   gamma_2(s) = gamma_2 + (s - n_g) |g><g|, and watch n_g-dot(s). Blocking means
   n_g-dot >= 0 as s -> 0 and n_g-dot <= 0 as s approaches its ceiling. If instead
   n_g-dot is roughly constant in s, there is no blocking and the flow will push
   straight through both ends.
"""
import itertools
import sys
import numpy as np
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


def apply_str(up, lo, v):
    for q in lo:
        v = aa(q, v)
    for p in reversed(up):
        v = ad(p, v)
    return v


def rdm(psi, k):
    g = np.zeros((M,) * (2 * k), dtype=complex)
    for t in itertools.product(*([range(M)] * (2 * k))):
        g[t] = np.vdot(psi, apply_str(t[:k], t[k:], psi))
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
psi = np.where(np.array([bin(s).count("1") for s in range(dim)]) == N,
               rng.normal(size=dim) + 1j * rng.normal(size=dim), 0.0)
psi = psi / np.linalg.norm(psi)
g2 = rdm(psi, 2)
g3 = rdm(psi, 3)


def sep_g3(g1, g2t):
    l2 = g2t - (np.einsum("pq,vw->pvqw", g1, g1) - np.einsum("pw,vq->pvqw", g1, g1))
    t = {("G", 1): g1, ("L", 2): l2}
    out = np.zeros((M,) * 6, dtype=complex)
    for sign, sub, keys in SEP3:
        out = out + sign * np.einsum(sub, *[t[k] for k in keys], optimize=True)
    return out


def term_contrib(g1, g2t, g3t, gvec):
    """Per-term contribution to <g|gamma_2-dot|g>, split by which RDM it carries."""
    eye = np.eye(M, dtype=complex)
    gam = {1: g1, 2: g2t, 3: g3t}
    out = []
    for tag, terms, Hop, pref in (("h", ONE, h, 1.0), ("V", TWO, V, 0.25)):
        for sign, facs in terms:
            subs = ["rs" if Hop.ndim == 2 else "rtsu"]
            ops = [Hop]
            rank = 0
            for n, r, u, l in facs:
                subs.append("".join(u) + "".join(l))
                ops.append(eye if n == "D" else gam[r])
                if n == "G":
                    rank = r
            T = sign * pref * np.einsum(",".join(subs) + "->pvqw", *ops, optimize=True)
            Tp = to_pair(1j * T)
            Tp = 0.5 * (Tp + Tp.conj().T)
            out.append((tag, rank, float(np.real(gvec.conj() @ Tp @ gvec))))
    return out


rho = 0.5 * (to_pair(g2) + to_pair(g2).conj().T)
w, U = np.linalg.eigh(rho)
occ = w[::-1]
Uo = U[:, ::-1]
g1 = contract_to_1rdm(g2, N)
nocc = int((occ > 1e-9).sum())
print("occupied geminals: %s" % np.array2string(occ[:nocc], precision=4))
print()

print("=== A. sources and drains at the true state (exact gamma_3) ===")
print("%4s %9s | %10s %10s %10s | %10s" % ("g", "n_g", "sources", "drains", "net",
                                           "ndot check"))
for j in range(nocc):
    gv = Uo[:, j]
    contrib = term_contrib(g1, g2, g3, gv)
    src = sum(c for _, _, c in contrib if c > 0)
    drn = sum(c for _, _, c in contrib if c < 0)
    Xp = to_pair(1j * (evaluate(ONE_C, h, {1: g1, 2: g2, 3: g3}, M, out_rank=4)
                       + 0.25 * evaluate(TWO_C, V, {1: g1, 2: g2, 3: g3}, M, out_rank=4)))
    Xp = 0.5 * (Xp + Xp.conj().T)
    print("%4d %9.4f | %10.4f %10.4f %10.4f | %10.4f"
          % (j, occ[j], src, drn, src + drn, float(np.real(gv.conj() @ Xp @ gv))))

print()
print("=== B. sweep the occupation of one geminal, closure gamma_3 = 0 ===")
print("blocking would mean ndot >= 0 near s = 0 and ndot <= 0 near the ceiling")
print()
for j in (0, 3, 5):
    gv = Uo[:, j]
    print("geminal %d (true n_g = %.4f):" % (j, occ[j]))
    print("   %8s %12s %12s %12s" % ("s", "ndot", "sources", "drains"))
    for s in (0.0, 0.25, 0.5, 1.0, 1.5, 2.0):
        rho_s = rho + (s - occ[j]) * np.outer(gv, gv.conj())
        g2s = from_pair(0.5 * (rho_s + rho_s.conj().T))
        g1s = contract_to_1rdm(g2s, N)
        z3 = np.zeros((M,) * 6, dtype=complex)
        contrib = term_contrib(g1s, g2s, z3, gv)
        src = sum(c for _, _, c in contrib if c > 0)
        drn = sum(c for _, _, c in contrib if c < 0)
        print("   %8.3f %12.5f %12.5f %12.5f" % (s, src + drn, src, drn))
    print()
