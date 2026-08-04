#!/usr/bin/env python3
"""Calibration: does rank 1 have the first-order blocking that rank 2 lacks?

At rank 2 every term of gamma_2-dot vanishes at a saturated geminal as gap^0.5, which
is stationarity of n_g rather than a multiplicative (c - n_g) factor. Rank 1 is the
control, because there the hole factor eta = delta - gamma is known to appear
explicitly and the ceiling n_p <= 1 is exact, universal, and attained on a whole face
of the feasible set rather than at one extremal point.

Build a state with an exactly saturated orbital, psi = a+_0 |chi> with chi an (N-1)
electron state on the remaining orbitals, so n_0 = 1. Displace toward a random
N-electron state and fit how each term of gamma_1-dot scales in the gap 1 - n_0.

    exponent 1.0 in the gap : a genuine multiplicative blocking factor
    exponent 0.5 in the gap : stationarity only, same as rank 2

The two-body part uses the EXACT gamma_2, so this is the exact gamma_1-dot, not a
closure. If first-order blocking exists anywhere it must show up here.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from tools.eom1 import one as ONE1, two as TWO1

M, N = 8, 4
dim = 1 << M


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
    return P1.conj() @ P1.T, (P2.conj() @ P2.T).reshape(M, M, M, M)


LBL = ["h%d" % i for i in range(len(ONE1))] + ["V%d" % i for i in range(len(TWO1))]


def terms(h, V, g1, g2, u):
    eye = np.eye(M, dtype=complex)
    gam = {1: g1, 2: g2}
    out = []
    for tl, Hop, pref, hsub in ((ONE1, h, 1.0, "rs"), (TWO1, V, 0.25, "rtsu")):
        for sign, facs in tl:
            subs = [hsub]; ops = [Hop]
            for n, r, up, lo in facs:
                subs.append("".join(up) + "".join(lo))
                ops.append(eye if n == "D" else gam[r])
            T = sign * pref * np.einsum(",".join(subs) + "->pq", *ops, optimize=True)
            T = 1j * T
            T = 0.5 * (T + T.conj().T)
            out.append(float(np.real(u.conj() @ T @ u)))
    return np.array(out)


r = np.random.default_rng(808)
A_ = r.normal(size=(M, M)) + 1j * r.normal(size=(M, M))
h = 0.25 * (A_ + A_.conj().T)
Bt = r.normal(size=(M,) * 4) + 1j * r.normal(size=(M,) * 4)
V = Bt - np.transpose(Bt, (1, 0, 2, 3))
V = V - np.transpose(V, (0, 1, 3, 2))
V = 0.25 * (V + np.transpose(V.conj(), (2, 3, 0, 1)))

popc = np.array([bin(s).count("1") for s in range(dim)])
chi = np.where((popc == N - 1) & (np.arange(dim) & 1 == 0),
               r.normal(size=dim) + 1j * r.normal(size=dim), 0.0)
chi = chi / np.linalg.norm(chi)
sat = ad(0, chi)
sat = sat / np.linalg.norm(sat)
mix = np.where(popc == N, r.normal(size=dim) + 1j * r.normal(size=dim), 0.0)
mix = mix / np.linalg.norm(mix)

g1s, _ = rdms(sat)
print("saturated state: natural occupations %s"
      % np.array2string(np.linalg.eigvalsh(g1s)[::-1], precision=6))
print()

TS = [0.01, 0.02, 0.04, 0.08]
gaps, rows = [], []
for t in TS:
    psi = np.cos(t) * sat + np.sin(t) * mix
    psi = psi / np.linalg.norm(psi)
    g1, g2 = rdms(psi)
    w, U = np.linalg.eigh(g1)
    gaps.append(1.0 - w[-1])
    rows.append(terms(h, V, g1, g2, U[:, -1]))
gaps = np.array(gaps); Y = np.array(rows)

print("displacement t   : %s" % "  ".join("%9.4f" % t for t in TS))
print("gap 1 - n_p      : %s" % "  ".join("%9.2e" % g for g in gaps))
print("gap / t          : %s" % "  ".join("%9.4f" % (g / t) for g, t in zip(gaps, TS)))
print("gap / t^2        : %s" % "  ".join("%9.4f" % (g / t ** 2) for g, t in zip(gaps, TS)))
print()
print("%6s | %s | %8s %8s" % ("term", "  ".join("%9s" % ("t=%.2f" % t) for t in TS),
                              "exp in t", "exp in gap"))
for i, lb in enumerate(LBL):
    y = Y[:, i]
    if np.abs(y).max() < 1e-12:
        print("%6s | %s | %8s %8s" % (lb, "  ".join("%9.1e" % v for v in y), "zero", "zero"))
        continue
    et = np.polyfit(np.log(TS), np.log(np.abs(y)), 1)[0]
    eg = np.polyfit(np.log(gaps), np.log(np.abs(y)), 1)[0]
    print("%6s | %s | %8.3f %8.3f" % (lb, "  ".join("%+9.2e" % v for v in y), et, eg))
tot = Y.sum(axis=1)
print("%6s | %s | %8.3f %8.3f" % ("SUM", "  ".join("%+9.2e" % v for v in tot),
      np.polyfit(np.log(TS), np.log(np.abs(tot)), 1)[0],
      np.polyfit(np.log(gaps), np.log(np.abs(tot)), 1)[0]))
print()
print("rank 2 reference (tools/why_vanishing.py): every term exp 1.0 in t, 0.5 in gap.")
