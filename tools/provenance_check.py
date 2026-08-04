#!/usr/bin/env python3
"""Provenance and consistency checks raised by John's question.

1. Does the lambda_3 = 0 reconstruction satisfy the contraction sum rule
   sum_x gamma_3[pvx,qwx] = (N-2) gamma_2[pv,qw]?  Drafted as a claim, never tested.
2. tools/eom2_coherent.py's docstring states H = sum h + (1/2) sum V^rt_su a^rt_su, but
   every consumer applies 1/4. Which prefactor reproduces the exact gamma_2-dot?
"""
import itertools
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from qdiss.symbolic import mkno
from tools.eom2_coherent import ONE, TWO

M, N = 6, 3
dim = 1 << M
UP3, LO3 = ("p", "v", "x"), ("q", "w", "y")
SEP3 = [(s, ",".join("".join(u) + "".join(l) for n, r, u, l in f) + "->pvxqwy",
         [(n, r) for n, r, u, l in f])
        for s, f in mkno.gamma_to_cumulants(UP3, LO3)
        if not any(n == "L" and r == 3 for n, r, u, l in f)]


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


r = np.random.default_rng(11)
popc = np.array([bin(s).count("1") for s in range(dim)])
psi = np.where(popc == N, r.normal(size=dim) + 1j * r.normal(size=dim), 0.0)
psi = psi / np.linalg.norm(psi)
Q1 = np.array([aa(p, psi) for p in range(M)])
g1 = Q1.conj() @ Q1.T
Q2 = np.array([aa(v, aa(p, psi)) for p in range(M) for v in range(M)])
g2 = (Q2.conj() @ Q2.T).reshape(M, M, M, M)
Q3 = np.array([aa(x, aa(v, aa(p, psi)))
               for p in range(M) for v in range(M) for x in range(M)])
g3 = (Q3.conj() @ Q3.T).reshape((M,) * 6)
l2 = g2 - (np.einsum("pq,vw->pvqw", g1, g1) - np.einsum("pw,vq->pvqw", g1, g1))
t = {("G", 1): g1, ("L", 2): l2}
g3s = np.zeros((M,) * 6, dtype=complex)
for s, sub, keys in SEP3:
    g3s = g3s + s * np.einsum(sub, *[t[k] for k in keys], optimize=True)

print("=== 1. contraction sum rule  sum_x gamma_3[pvx,qwx] = (N-2) gamma_2 ===")
tgt = (N - 2) * g2
for nm, G in (("exact gamma_3", g3), ("lambda_3 = 0", g3s),
              ("gamma_3 = 0", np.zeros((M,) * 6, dtype=complex))):
    err = np.abs(np.einsum("pvxqwx->pvqw", G) - tgt).max()
    print("   %-14s  violation %.4e   (scale |(N-2) gamma_2| = %.4f)"
          % (nm, err, np.abs(tgt).max()))
print()
print("   also: sum_x gamma_2[px,qx] = (N-1) gamma_1 violation %.3e   tr gamma_1 = %.4f"
      % (np.abs(np.einsum("pxqx->pq", g2) - (N - 1) * g1).max(), np.trace(g1).real))
print()

print("=== 2. which two-body prefactor reproduces the exact gamma_2-dot? ===")
h = 0.25 * (lambda A: A + A.conj().T)(r.normal(size=(M, M)) + 1j * r.normal(size=(M, M)))
Bt = r.normal(size=(M,) * 4) + 1j * r.normal(size=(M,) * 4)
V = Bt - np.transpose(Bt, (1, 0, 2, 3))
V = V - np.transpose(V, (0, 1, 3, 2))
V = 0.25 * (V + np.transpose(V.conj(), (2, 3, 0, 1)))
Hm = np.zeros((dim, dim), dtype=complex)
basis = np.eye(dim, dtype=complex)
for p in range(M):
    for q in range(M):
        if abs(h[p, q]) < 1e-14: continue
        for k in range(dim):
            Hm[:, k] += h[p, q] * ad(p, aa(q, basis[:, k]))
for p, v, q, w in itertools.product(*([range(M)] * 4)):
    if abs(V[p, v, q, w]) < 1e-14: continue
    for k in range(dim):
        Hm[:, k] += 0.25 * V[p, v, q, w] * ad(p, ad(v, aa(w, aa(q, basis[:, k]))))
Hm = 0.5 * (Hm + Hm.conj().T)


def exact_dot():
    out = np.zeros((M,) * 4, dtype=complex)
    Hp = Hm @ psi
    for p, v, q, w in itertools.product(*([range(M)] * 4)):
        O = ad(p, ad(v, aa(w, aa(q, psi))))
        OH = ad(p, ad(v, aa(w, aa(q, Hp))))
        out[p, v, q, w] = 1j * (np.vdot(Hp, O) - np.vdot(psi, OH))
    return out


ex = exact_dot()
print("   anti-Hermitian part of H: %.2e   scale of exact gamma_2-dot: %.4f" % (0.0, np.abs(ex).max()))
eye = np.eye(M, dtype=complex)
for pf in (0.5, 0.25):
    gam = {1: g1, 2: g2, 3: g3}
    tot = np.zeros((M,) * 4, dtype=complex)
    for tl, Hop, p_, hs in ((ONE, h, 1.0, "rs"), (TWO, V, pf, "rtsu")):
        for s, facs in tl:
            subs = [hs]; ops = [Hop]
            for n, rk, u, l in facs:
                subs.append("".join(u) + "".join(l))
                ops.append(eye if n == "D" else gam[rk])
            tot = tot + s * p_ * np.einsum(",".join(subs) + "->pvqw", *ops, optimize=True)
    print("   prefactor %.2f : max error vs exact Fock space = %.4e" % (pf, np.abs(1j * tot - ex).max()))
