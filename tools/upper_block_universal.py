#!/usr/bin/env python3
"""Is there a universal root in n_g-dot(s) at the geminal ceiling?

tools/upper_block.py established the degrees: with gamma_3 = 0 the geminal flow is
linear in the occupation s, and with lambda_3 = 0 it is cubic, so only the latter can
host a root at a nonzero ceiling. This asks whether it does.

A blocking factor is a structural property of the equation, not of a particular
Hamiltonian. If n_g-dot carries a factor (c - n_g), then n_g-dot(s) vanishes at s = c
for EVERY h and V. So sample many random Hamiltonians, solve for the real roots of the
cubic, and see whether any root sits at the same place each time. A root that wanders
with H is an accident of that H; a root pinned at 2.5 = N(M-N+2)/M, or at N = 3, is a
blocking factor.

Also reports max |n_g-dot| over the whole sweep for the gamma_3 = 0 closure, to pin
down whether it is merely small or identically zero.
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
    for t in itertools.product(*([range(M)] * 4)):
        g[t] = np.vdot(psi, ad(t[0], ad(t[1], aa(t[3], aa(t[2], psi)))))
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


def sep_g3(g1, g2t):
    l2 = g2t - (np.einsum("pq,vw->pvqw", g1, g1) - np.einsum("pw,vq->pvqw", g1, g1))
    t = {("G", 1): g1, ("L", 2): l2}
    out = np.zeros((M,) * 6, dtype=complex)
    for sign, sub, keys in SEP3:
        out = out + sign * np.einsum(sub, *[t[k] for k in keys], optimize=True)
    return out


def ndot(h, V, g1, g2t, g3t, gvec):
    eye = np.eye(M, dtype=complex)
    gam = {1: g1, 2: g2t, 3: g3t}
    tot = 0.0
    for terms, Hop, pref, hsub in ((ONE, h, 1.0, "rs"), (TWO, V, 0.25, "rtsu")):
        for sign, facs in terms:
            subs = [hsub]; ops = [Hop]
            for n, r, u, l in facs:
                subs.append("".join(u) + "".join(l))
                ops.append(eye if n == "D" else gam[r])
            T = sign * pref * np.einsum(",".join(subs) + "->pvqw", *ops, optimize=True)
            Tp = to_pair(1j * T)
            Tp = 0.5 * (Tp + Tp.conj().T)
            tot += float(np.real(gvec.conj() @ Tp @ gvec))
    return tot


SAMP = np.linspace(0.0, 3.6, 13)
popc = np.array([bin(s).count("1") for s in range(dim)])
print("ceiling N(M-N+2)/M = %.4f      asymptotic N = %d" % (CEIL, N))
print()
print("%4s %6s | %11s | %s" % ("trial", "n_g", "max|ndot| g3=0", "real roots of ndot(s), lambda_3 = 0"))
allroots = []
frozen = 0.0
for trial in range(14):
    rng = np.random.default_rng(1000 + trial)
    A_ = rng.normal(size=(M, M)) + 1j * rng.normal(size=(M, M))
    h = 0.25 * (A_ + A_.conj().T)
    B = rng.normal(size=(M,) * 4) + 1j * rng.normal(size=(M,) * 4)
    V = B - np.transpose(B, (1, 0, 2, 3))
    V = V - np.transpose(V, (0, 1, 3, 2))
    V = 0.25 * (V + np.transpose(V.conj(), (2, 3, 0, 1)))
    psi = np.where(popc == N, rng.normal(size=dim) + 1j * rng.normal(size=dim), 0.0)
    psi = psi / np.linalg.norm(psi)
    rho = 0.5 * (to_pair(rdm2(psi)) + to_pair(rdm2(psi)).conj().T)
    w, U = np.linalg.eigh(rho)
    occ, gv = w[::-1][1], U[:, ::-1][:, 1]

    z3 = np.zeros((M,) * 6, dtype=complex)
    ys, yl = [], []
    for s in SAMP:
        rho_s = rho + (s - occ) * np.outer(gv, gv.conj())
        g2s = from_pair(0.5 * (rho_s + rho_s.conj().T))
        g1s = contract_to_1rdm(g2s, N)
        ys.append(ndot(h, V, g1s, g2s, z3, gv))
        yl.append(ndot(h, V, g1s, g2s, sep_g3(g1s, g2s), gv))
    frozen = max(frozen, np.abs(ys).max())
    c = np.polyfit(SAMP, yl, 3)
    rt = sorted(z.real for z in np.roots(c) if abs(z.imag) < 1e-7)
    allroots.extend(r for r in rt if 0.0 < r < 5.0)
    print("%4d %6.3f | %11.2e | %s" % (trial, occ, np.abs(ys).max(),
          "  ".join("%7.3f" % r for r in rt) if rt else "no real root"))

print()
ar = np.array(allroots)
print("real roots in (0,5) across trials: n=%d  mean %.3f  std %.3f  min %.3f  max %.3f"
      % (len(ar), ar.mean(), ar.std(), ar.min(), ar.max()) if len(ar) else "none")
for target, name in ((CEIL, "ceiling 2.500"), (float(N), "N = 3.000")):
    near = int((np.abs(ar - target) < 0.05).sum()) if len(ar) else 0
    print("   roots within 0.05 of %s : %d of %d trials" % (name, near, 14))
