#!/usr/bin/env python3
"""Certify the derived full-cumulant 1-2-2 closed EOM against an EXACT
   Fock-space evaluation.  No Wolfram, no symbolic trust.

   Claim under test, index by index:

     < [ [ {a^p_q}, {a^rt_su} ], {a^vx_wy} ] >  ==  sum over the 368 derived
     terms, each a product of delta / gamma / lambda2 / lambda3.

   The 368 terms are obtained from the 1312-term commutator by taking the MK
   expectation value:  <{a_1}> = 0,  <{a_2}> = lambda2,  <{a_3}> = lambda3.

   Conventions are inherited from tools/verify_gwt.py, which verifies them.
"""
import sys, itertools, numpy as np
sys.path.insert(0, 'tools')
from term_audit import parse, canon

M, N = 6, 3
rng = np.random.default_rng(20260802)

# ---------------- Fock space ----------------
dim = 1 << M
def cre(p):
    A = np.zeros((dim, dim))
    for s in range(dim):
        if (s >> p) & 1: continue
        A[s | (1 << p), s] = (-1) ** bin(s & ((1 << p) - 1)).count("1")
    return A
Ad = [cre(p) for p in range(M)]
Aa = [a.T.copy() for a in Ad]

idx = [s for s in range(dim) if bin(s).count("1") == N]
psi = np.zeros(dim)
for k, s in enumerate(idx): psi[s] = rng.normal()
psi /= np.linalg.norm(psi)
ev = lambda O: float(psi @ (O @ psi))

def astr(up, lo):
    O = np.eye(dim)
    for p in up: O = O @ Ad[p]
    for q in reversed(lo): O = O @ Aa[q]
    return O

R = range(M)
g1 = np.array([[ev(astr([p], [q])) for q in R] for p in R])
g2 = np.zeros((M,) * 4)
for t in itertools.product(R, repeat=4): g2[t] = ev(astr(list(t[:2]), list(t[2:])))
g3 = np.zeros((M,) * 6)
for t in itertools.product(R, repeat=6): g3[t] = ev(astr(list(t[:3]), list(t[3:])))

I = np.eye(M)
l2 = g2 - np.einsum('pr,qs->pqrs', g1, g1) + np.einsum('ps,qr->pqrs', g1, g1)

def wedge3(g, l2):
    A = np.zeros((M,) * 6); B = np.zeros((M,) * 6)
    perms = [(p, int(round(np.linalg.det(np.eye(3)[list(p)])))) for p in itertools.permutations(range(3))]
    for up in itertools.product(R, repeat=3):
        for lo in itertools.product(R, repeat=3):
            a = sum(sg * g[up[0], lo[pi[0]]] * g[up[1], lo[pi[1]]] * g[up[2], lo[pi[2]]]
                    for pi, sg in perms)
            b = 0.0
            for i in range(3):
                for j in range(3):
                    ru = [up[k] for k in range(3) if k != i]
                    rl = [lo[k] for k in range(3) if k != j]
                    b += (-1) ** (i + j) * g[up[i], lo[j]] * l2[ru[0], ru[1], rl[0], rl[1]]
            A[up + lo] = a; B[up + lo] = b
    return A, B
A3, B3 = wedge3(g1, l2)
l3 = g3 - A3 - B3
print(f"reference: M={M} spin orbitals, N={N} electrons, random correlated state")
print(f"  |lambda2| = {np.linalg.norm(l2):.5f}   |lambda3| = {np.linalg.norm(l3):.5f}")

# ---------------- MK normal-ordered operators ----------------
def gno1(p, q): return astr([p], [q]) - g1[p, q] * np.eye(dim)
def gno2(p, q, r, s):
    O = astr([p, q], [r, s]) - g2[p, q, r, s] * np.eye(dim)
    O -= g1[p, r] * gno1(q, s); O += g1[p, s] * gno1(q, r)
    O += g1[q, r] * gno1(p, s); O -= g1[q, s] * gno1(p, r)
    return O

# ---------------- the derived closed EOM ----------------
def unpack(ts):
    out = []
    for n, r, up, lo in ts:
        if n == 'D' and r > 1: out += [('D', 1, (a,), (b,)) for a, b in zip(up, lo)]
        else: out.append((n, r, up, lo))
    return out

base = {}
for sg, ts in parse('results/stage_joined_raw.m'):
    s2, k = canon(sg, unpack(ts)); base[k] = base.get(k, 0) + s2
closed = {}
for k, v in base.items():
    if v == 0: continue
    # MK generalized normal ordering: <{a_k}> = 0 for EVERY k >= 1.  This is the
    # defining property of MK GNO and is verified numerically in
    # tools/verify_gwt.py -- it is NOT ordinary normal ordering, where <{a_2}>
    # would be lambda2.  So every residual-operator term drops, and the closed
    # equation is exactly the scalar slice.
    if any(t[0] == 'ta' for t in k):
        continue
    closed[k] = closed.get(k, 0) + v
closed = {k: v for k, v in closed.items() if v != 0}
print(f"  derived closed EOM: {len(closed)} terms")

TENSOR = {'D': lambda *i: I[i], 'G': lambda *i: g1[i],
          'L': lambda *i: (l2[i] if len(i) == 4 else l3[i])}

def evaluate(term_dict, sub):
    tot = 0.0
    for k, c in term_dict.items():
        val = float(c)
        for n, r, up, lo in k:
            i = tuple(sub[a] for a in up) + tuple(sub[a] for a in lo)
            val *= TENSOR[n](*i)
            if val == 0.0: break
        tot += val
    return tot

EXT = ['p', 'q', 'r', 't', 's', 'u', 'v', 'x', 'w', 'y']
def exact(sub):
    A = gno1(sub['p'], sub['q'])
    B = gno2(sub['r'], sub['t'], sub['s'], sub['u'])
    C = gno2(sub['v'], sub['x'], sub['w'], sub['y'])
    AB = A @ B - B @ A
    return ev(AB @ C - C @ AB)

ntrial = int(sys.argv[1]) if len(sys.argv) > 1 else 60
worst = 0.0; scale = 0.0; fails = 0
for n in range(ntrial):
    sub = {a: int(rng.integers(M)) for a in EXT}
    e = exact(sub); d = evaluate(closed, sub)
    worst = max(worst, abs(e - d)); scale = max(scale, abs(e))
    if abs(e - d) > 1e-9:
        fails += 1
        if fails <= 3:
            print(f"  MISMATCH {sub}: exact {e:+.9f}  derived {d:+.9f}")
print(f"\n{ntrial} random index tuples: max |exact - derived| = {worst:.3e}"
      f"   (typical magnitude {scale:.3f})")
print("VERDICT:", "EXACT MATCH" if worst < 1e-9 else f"MISMATCH in {fails} tuples")

# and what the old generator would have given
old = {k: v for k, v in closed.items() if not any(t[0] == 'L' for t in k)}
errs = []
for n in range(ntrial):
    sub = {a: int(rng.integers(M)) for a in EXT}
    errs.append(abs(exact(sub) - evaluate(old, sub)))
errs = np.array(errs)
print(f"\nNoCumulants ({len(old)} terms) against the same exact values:")
print(f"  max error {errs.max():.4f}, mean error {errs.mean():.4f}"
      f"  -- nonzero: the deleted terms are not a small correction")
