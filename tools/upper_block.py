#!/usr/bin/env python3
"""Search for an upper blocking factor in the geminal source terms.

Lower blocking is established: every source carries gamma, so it dies as a geminal
empties. The question here is the other end -- does any term carry a factor that
vanishes as n_g approaches its ceiling, Coleman's N(M-N+2)/M?

A degree argument says where to look. Deform gamma_2(s) = gamma_2 + (s-n_g)|g><g|
and every RDM becomes a polynomial in s: gamma_1 via contraction is degree 1, and a
term's contribution has degree equal to its number of gamma factors. A term linear
in s that vanishes at s = 0 is c*s, which has no second root -- so it can block at
the bottom or at a ceiling, never both. Upper blocking therefore cannot come from
the h terms or from the gamma_3 = 0 truncation, only from the gamma_3 terms, whose
separable reconstruction is quadratic and cubic in s.

So: fit every term's contribution to <g|gamma_2-dot|g> as a polynomial in s,
confirm the degrees, and look for roots sitting at the ceiling.
"""
import itertools
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from qdiss.symbolic import mkno
from qdiss.representability.dqg import contract_to_1rdm
from tools.eom2_coherent import ONE, TWO

M, N = 6, 3
CEIL = N * (M - N + 2) / M
dim = 1 << M
PAIRS = [(p, v) for p in range(M) for v in range(p + 1, M)]
NP = len(PAIRS)
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
    g = np.zeros((M,) * (2 * k), dtype=complex)
    for t in itertools.product(*([range(M)] * (2 * k))):
        v = psi
        for q in t[k:][::-1]:
            v = aa(q, v)
        for p in t[:k]:
            v = ad(p, v)
        g[t] = np.vdot(psi, v)
    return g


def to_pair(T):
    return np.array([[2.0 * T[p, v, q, w] for (q, w) in PAIRS] for (p, v) in PAIRS])


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

TERMS = [("h", ONE, h, 1.0, "rs"), ("V", TWO, V, 0.25, "rtsu")]

print("term inventory for gamma_2-dot, by the RDM rank each term carries")
for tag, terms, _, _, _ in TERMS:
    cnt = {}
    for sign, facs in terms:
        r = max([r for n, r, u, l in facs if n == "G"] + [0])
        cnt[r] = cnt.get(r, 0) + 1
    print("   %s : %d terms, %s" % (tag, len(terms),
          ", ".join("%d carry gamma_%d" % (v, k) for k, v in sorted(cnt.items()))))
print()


def sep_g3(g1, g2t):
    l2 = g2t - (np.einsum("pq,vw->pvqw", g1, g1) - np.einsum("pw,vq->pvqw", g1, g1))
    t = {("G", 1): g1, ("L", 2): l2}
    out = np.zeros((M,) * 6, dtype=complex)
    for sign, sub, keys in SEP3:
        out = out + sign * np.einsum(sub, *[t[k] for k in keys], optimize=True)
    return out


def contribs(g1, g2t, g3t, gvec):
    eye = np.eye(M, dtype=complex)
    gam = {1: g1, 2: g2t, 3: g3t}
    out = []
    for tag, terms, Hop, pref, hsub in TERMS:
        for ti, (sign, facs) in enumerate(terms):
            subs = [hsub]; ops = [Hop]; rank = 0
            for n, r, u, l in facs:
                subs.append("".join(u) + "".join(l))
                ops.append(eye if n == "D" else gam[r])
                if n == "G":
                    rank = r
            T = sign * pref * np.einsum(",".join(subs) + "->pvqw", *ops, optimize=True)
            Tp = to_pair(1j * T)
            Tp = 0.5 * (Tp + Tp.conj().T)
            out.append(("%s%d" % (tag, ti), rank,
                        float(np.real(gvec.conj() @ Tp @ gvec))))
    return out


rho = 0.5 * (to_pair(g2) + to_pair(g2).conj().T)
w, U = np.linalg.eigh(rho)
occ, Uo = w[::-1], U[:, ::-1]
nocc = int((occ > 1e-9).sum())
print("ceiling N(M-N+2)/M = %.4f     occupied geminals %s"
      % (CEIL, np.array2string(occ[:nocc], precision=4)))
print()

SAMP = np.linspace(0.0, 1.25 * CEIL, 13)

for closure in ("gamma_3 = 0", "lambda_3 = 0"):
    print("=== closure %s ===" % closure)
    print("%5s %7s | %6s %10s %10s %10s | %s"
          % ("gem", "n_g", "deg", "ndot(0)", "ndot(ceil)", "ndot(1.25c)", "term roots at ceiling"))
    for j in range(min(nocc, 3)):
        gv = Uo[:, j]
        rows = []
        for s in SAMP:
            rho_s = rho + (s - occ[j]) * np.outer(gv, gv.conj())
            g2s = from_pair(0.5 * (rho_s + rho_s.conj().T))
            g1s = contract_to_1rdm(g2s, N)
            g3s = (np.zeros((M,) * 6, dtype=complex) if closure == "gamma_3 = 0"
                   else sep_g3(g1s, g2s))
            rows.append(contribs(g1s, g2s, g3s, gv))
        names = [r[0] for r in rows[0]]
        ranks = [r[1] for r in rows[0]]
        Y = np.array([[r[2] for r in row] for row in rows])
        tot = Y.sum(axis=1)
        deg = 0
        for d in range(1, 5):
            c = np.polyfit(SAMP, tot, d)
            if np.abs(np.polyval(c, SAMP) - tot).max() < 1e-8 * max(1, np.abs(tot).max()):
                deg = d; break
        hits = []
        for k, nm in enumerate(names):
            if np.abs(Y[:, k]).max() < 1e-10:
                continue
            dk = 0
            for d in range(1, 5):
                ck = np.polyfit(SAMP, Y[:, k], d)
                if np.abs(np.polyval(ck, SAMP) - Y[:, k]).max() < 1e-8 * max(1, np.abs(Y[:, k]).max()):
                    dk = d; break
            if dk >= 2:
                rt = np.roots(np.polyfit(SAMP, Y[:, k], dk))
                for z in rt:
                    if abs(z.imag) < 1e-6 and abs(z.real - CEIL) < 0.15:
                        hits.append("%s@%.3f" % (nm, z.real))
        print("%5d %7.4f | %6d %10.4f %10.4f %10.4f | %s"
              % (j, occ[j], deg, np.interp(0.0, SAMP, tot),
                 np.interp(CEIL, SAMP, tot), tot[-1],
                 ",".join(hits) if hits else "none"))
    print()
