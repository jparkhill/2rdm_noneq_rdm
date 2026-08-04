#!/usr/bin/env python3
"""Why does every term vanish at saturation -- a blocking factor, or stationarity?

A multiplicative blocking factor (c - n_g) makes each term LINEAR in the gap. A
critical point of n_g on the representable set makes ndot linear in the displacement
t while the gap falls as t^2, so each term goes as sqrt(gap). Fitting the exponent
per term separates the two.

The aggregate already scales as sqrt. This asks whether every individual term does,
which would say they vanish for one common reason, or whether the exponents differ,
which would say several mechanisms are at work.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from tools.eom2_coherent import ONE, TWO

M, N = 8, 4
CEIL = N * (M - N + 2) / M
dim = 1 << M
PAIRS = [(p, v) for p in range(M) for v in range(p + 1, M)]


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


def rdms(psi):
    P1 = np.array([aa(p, psi) for p in range(M)])
    P2 = np.array([aa(v, aa(p, psi)) for p in range(M) for v in range(M)])
    P3 = np.array([aa(x, aa(v, aa(p, psi)))
                   for p in range(M) for v in range(M) for x in range(M)])
    return (P1.conj() @ P1.T, (P2.conj() @ P2.T).reshape(M, M, M, M),
            (P3.conj() @ P3.T).reshape((M,) * 6))


def to_pair(T):
    return np.array([[2.0 * T[p, v, q, w] for (q, w) in PAIRS] for (p, v) in PAIRS])


LBL = ["h%d" % i for i in range(len(ONE))] + ["V%d" % i for i in range(len(TWO))]


def terms(h, V, g1, g2, g3, gv):
    eye = np.eye(M, dtype=complex)
    gam = {1: g1, 2: g2, 3: g3}
    out = []
    for tl, Hop, pref, hsub in ((ONE, h, 1.0, "rs"), (TWO, V, 0.25, "rtsu")):
        for sign, facs in tl:
            subs = [hsub]; ops = [Hop]
            for n, r, u, l in facs:
                subs.append("".join(u) + "".join(l))
                ops.append(eye if n == "D" else gam[r])
            T = sign * pref * np.einsum(",".join(subs) + "->pvqw", *ops, optimize=True)
            Tp = to_pair(1j * T); Tp = 0.5 * (Tp + Tp.conj().T)
            out.append(float(np.real(gv.conj() @ Tp @ gv)))
    return np.array(out)


vac = np.zeros(dim, dtype=complex); vac[0] = 1.0
agp = vac
for _ in range(N // 2):
    acc = np.zeros_like(agp)
    for k in range(M // 2):
        acc = acc + ad(2 * k, ad(2 * k + 1, agp))
    agp = acc
agp = agp / np.linalg.norm(agp)
popc = np.array([bin(s).count("1") for s in range(dim)])
r = np.random.default_rng(4242)
A_ = r.normal(size=(M, M)) + 1j * r.normal(size=(M, M))
h = 0.25 * (A_ + A_.conj().T)
Bt = r.normal(size=(M,) * 4) + 1j * r.normal(size=(M,) * 4)
V = Bt - np.transpose(Bt, (1, 0, 2, 3))
V = V - np.transpose(V, (0, 1, 3, 2))
V = 0.25 * (V + np.transpose(V.conj(), (2, 3, 0, 1)))
mix = np.where(popc == N, r.normal(size=dim) + 1j * r.normal(size=dim), 0.0)
mix = mix / np.linalg.norm(mix)

TS = [0.01, 0.02, 0.04, 0.08]
gaps, tvals, ndots = [], [], []
for t in TS:
    psi = np.cos(t) * agp + np.sin(t) * mix
    psi = psi / np.linalg.norm(psi)
    g1, g2, g3 = rdms(psi)
    rho = 0.5 * (to_pair(g2) + to_pair(g2).conj().T)
    w, U = np.linalg.eigh(rho)
    gaps.append(CEIL - w[-1])
    tvals.append(t)
    ndots.append(terms(h, V, g1, g2, g3, U[:, -1]))
gaps = np.array(gaps); Y = np.array(ndots)

print("displacement t      : %s" % "  ".join("%9.4f" % x for x in tvals))
print("gap c - n_g         : %s" % "  ".join("%9.2e" % x for x in gaps))
print("gap / t^2           : %s" % "  ".join("%9.4f" % (g / t ** 2) for g, t in zip(gaps, tvals)))
print("  -> gap is quadratic in t, i.e. n_g is at a critical point")
print()
print("%6s | %s | %8s %8s" % ("term", "  ".join("%9s" % ("t=%.2f" % t) for t in TS),
                              "exp in t", "exp in gap"))
for i, lb in enumerate(LBL):
    y = Y[:, i]
    if np.abs(y).max() < 1e-12:
        print("%6s | %s | %8s %8s" % (lb, "  ".join("%9.1e" % v for v in y), "zero", "zero"))
        continue
    et = np.polyfit(np.log(tvals), np.log(np.abs(y)), 1)[0]
    eg = np.polyfit(np.log(gaps), np.log(np.abs(y)), 1)[0]
    print("%6s | %s | %8.3f %8.3f" % (lb, "  ".join("%+9.2e" % v for v in y), et, eg))
tot = Y.sum(axis=1)
print("%6s | %s | %8.3f %8.3f" % ("SUM", "  ".join("%+9.2e" % v for v in tot),
      np.polyfit(np.log(tvals), np.log(np.abs(tot)), 1)[0],
      np.polyfit(np.log(gaps), np.log(np.abs(tot)), 1)[0]))
print()
print("exponent 1 in t and 0.5 in the gap = stationarity of n_g, no blocking factor.")
print("exponent 1 in the GAP would be a genuine multiplicative (c - n_g) factor.")
