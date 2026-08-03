#!/usr/bin/env python3
"""Derive Q^{ij}_{kl} = <a_i a_j a+_l a+_k> with the verified normal-ordering
engine, check it against explicit Fock space, and emit a correct build_Q.

dqg.py's build_Q disagrees with the exact two-hole matrix by 3.9e-01 and reports
a min eigenvalue of -1.55 on an exactly N-representable 2-RDM. Rather than
hand-correct the index order, generate the expression.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from qdiss.symbolic import mkno
from qdiss.representability.dqg import contract_to_1rdm

M, N = 6, 3

# a_i a_j a+_l a+_k  as an operator word
word = [("a", "i"), ("a", "j"), ("c", "l"), ("c", "k")]
terms = []
for s, deltas, cre, ann in mkno.normal_order(word):
    if len(cre) != len(ann):
        continue
    up = tuple(cre)
    lo = tuple(reversed(ann))
    facs = deltas + ((("G", len(up), up, lo),) if up else ())
    terms.append((s, facs))
terms = mkno.collect(terms)

print("Q^{ij}_{kl} = <a_i a_j a+_l a+_k>  :  %d terms" % len(terms))
for s, f in terms:
    parts = []
    for n, r, u, l in f:
        parts.append(("d" if n == "D" else "g") + str(r) + "[" + "".join(u) + ";" + "".join(l) + "]")
    print("   %+d  %s" % (s, " ".join(parts)))

# ---- numerical validation ----
rng = np.random.default_rng(2026)
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


psi = np.where(np.array([bin(s).count("1") for s in range(dim)]) == N,
               rng.normal(size=dim) + 1j * rng.normal(size=dim), 0.0)
psi = psi / np.linalg.norm(psi)

g2 = np.zeros((M,) * 4, dtype=complex)
Qx = np.zeros((M,) * 4, dtype=complex)
for i in range(M):
    for j in range(M):
        for k in range(M):
            for l in range(M):
                g2[i, j, k, l] = np.vdot(psi, ad(i, ad(j, aa(l, aa(k, psi)))))
                Qx[i, j, k, l] = np.vdot(psi, aa(i, aa(j, ad(l, ad(k, psi)))))
g1 = contract_to_1rdm(g2, N)

# build Q from the generated terms
eye = np.eye(M, dtype=complex)
OUT = "ijkl"
Qb = np.zeros((M,) * 4, dtype=complex)
for sign, facs in terms:
    subs, ops = [], []
    for n, r, u, l in facs:
        subs.append("".join(u) + "".join(l))
        ops.append(eye if n == "D" else (g1 if r == 1 else g2))
    Qb = Qb + sign * np.einsum(",".join(subs) + "->" + OUT, *ops, optimize=True)

print()
print("generated Q vs exact : %.3e" % np.abs(Qb - Qx).max())


def mineig(T):
    Mx = T.reshape(M * M, M * M)
    Mx = 0.5 * (Mx + Mx.conj().T)
    return float(np.linalg.eigvalsh(Mx).min())


print("min eig, exact Q     : %+.3e" % mineig(Qx))
print("min eig, generated Q : %+.3e" % mineig(Qb))
