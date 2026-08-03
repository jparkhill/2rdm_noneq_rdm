#!/usr/bin/env python3
"""Term-by-term boundary flux, and the geminal-basis structure of the exact EOM.

Part A: decompose <v|Ddot|v> and <v|Qdot|v> term by term over the 4 one-body
gamma_2 terms, the 4 two-body gamma_2 terms and the 8 two-body gamma_3 terms, so
we can see which individual terms carry the outward flux rather than only the
total.

Part B: the question John actually poses. The D-null directions ARE the empty
geminals, so <v|Ddot|v> is n_g-dot at n_g = 0. Under EXACT dynamics that must be
>= 0 (the true 2-RDM never leaves the cone). Look at the exact ndot in the
geminal basis and ask whether it has gain/loss structure with blocking factors --
the rank-2 analogue of ndot_p = sum_q [W n_q (1-n_p) - W n_p (1-n_q)].
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from qdiss.symbolic import mkno
from qdiss.dynamics.coherent2 import compile_terms
from qdiss.representability.dqg import contract_to_1rdm, build_Q
from tools.eom2_coherent import ONE, TWO

M, N = 6, 3
dim = 1 << M
eye = np.eye(M, dtype=complex)


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
    idx = [range(M)] * (2 * k)
    g = np.zeros((M,) * (2 * k), dtype=complex)
    import itertools
    for t in itertools.product(*idx):
        up, lo = t[:k], t[k:]
        v = psi
        for q in lo:
            v = aa(q, v)
        for p in reversed(up):
            v = ad(p, v)
        g[t] = np.vdot(psi, v)
    return g


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

g2 = rdm(psi, 2)
g3 = rdm(psi, 3)
g1 = contract_to_1rdm(g2, N)


def one_term(sign, facs, H, gam):
    subs, ops = [], []
    subs.append("rs" if H.ndim == 2 else "rtsu")
    ops.append(H)
    for n, r, u, l in facs:
        subs.append("".join(u) + "".join(l))
        ops.append(eye if n == "D" else gam[r])
    return sign * np.einsum(",".join(subs) + "->pvqw", *ops, optimize=True)


def qdot_from(dg2):
    dg1 = contract_to_1rdm(dg2, N)
    return (- np.einsum("ki,lj->ijkl", eye, dg1)
            + np.einsum("kj,li->ijkl", eye, dg1)
            + np.einsum("li,kj->ijkl", eye, dg1)
            - np.einsum("lj,ki->ijkl", eye, dg1)
            + np.einsum("klij->ijkl", dg2))


def mat(T):
    A_ = T.reshape(M * M, M * M)
    return 0.5 * (A_ + A_.conj().T)


gam = {1: g1, 2: g2, 3: g3}
wD, UD = np.linalg.eigh(mat(g2))
nullD = UD[:, wD < 1e-9]
wQ, UQ = np.linalg.eigh(mat(build_Q(g2, g1)))
nullQ = UQ[:, wQ < 1e-9]

print("=== Part A: term-by-term boundary flux (exact gamma_3) ===")
print("gamma_2 rank %d, so %d empty geminals; Q has %d null directions"
      % (len(wD) - nullD.shape[1], nullD.shape[1], nullQ.shape[1]))
print()
print("%-28s %14s %14s" % ("term", "min <v|Ddot|v>", "min <v|Qdot|v>"))

rows = []
for i, (s, f) in enumerate(ONE):
    dg2 = 1j * one_term(s, f, h, gam)
    rows.append(("1body g2 #%d" % i, dg2))
for i, (s, f) in enumerate(TWO):
    rank = max([r for n, r, u, l in f if n == "G"] or [0])
    dg2 = 1j * 0.25 * one_term(s, f, V, gam)
    rows.append(("2body g%d #%d" % (rank, i), dg2))

totD = np.zeros((M,) * 4, dtype=complex)
for label, dg2 in rows:
    fD = np.real(np.einsum("ai,ab,bi->i", nullD.conj(), mat(dg2), nullD))
    fQ = np.real(np.einsum("ai,ab,bi->i", nullQ.conj(), mat(qdot_from(dg2)), nullQ))
    print("%-28s %14.4e %14.4e" % (label, fD.min(), fQ.min()))
    totD = totD + dg2

fD = np.real(np.einsum("ai,ab,bi->i", nullD.conj(), mat(totD), nullD))
fQ = np.real(np.einsum("ai,ab,bi->i", nullQ.conj(), mat(qdot_from(totD)), nullQ))
print("%-28s %14.4e %14.4e" % ("TOTAL (exact dynamics)", fD.min(), fQ.min()))
print()

print("=== Part B: geminal occupations under exact dynamics ===")
occ = wD[::-1]
print("geminal occupations:", np.array2string(occ[occ > 1e-10], precision=6))
print("empty geminals: %d" % int((wD < 1e-9).sum()))
ndot = np.real(np.einsum("ai,ab,bi->i", UD.conj(), mat(totD), UD))
print()
print("%8s %12s %12s" % ("n_g", "ndot_g", "sign at n=0"))
order = np.argsort(-wD)
for j in order:
    tag = ""
    if wD[j] < 1e-9:
        tag = "OK" if ndot[j] > -1e-9 else "NEGATIVE"
    print("%8.6f %12.4e %12s" % (max(wD[j], 0.0), ndot[j], tag))
