#!/usr/bin/env python3
"""The coherent 2-RDM EOM: i d<a^pv_qw>/dt = <[H, a^pv_qw]>, no dissipator.

H = sum h^r_s a^r_s  +  (1/4) sum V^rt_su a^rt_su, V antisymmetrized

Both commutators are derived with the machinery verified in tools/eom1.py and
each is checked against explicit Fock-space matrices before use.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from qdiss.symbolic import mkno
from gno import Fock


def bare_product(u1, l1, u2, l2):
    word = mkno.ops_of(u1, l1) + mkno.ops_of(u2, l2)
    out = []
    for s, deltas, cre, ann in mkno.normal_order(word):
        if len(cre) == len(ann):
            out.append((s, deltas, tuple(cre), tuple(reversed(ann))))
    return out


def commutator_expectation(u1, l1, u2, l2):
    terms = []
    for s, d, up, lo in bare_product(u1, l1, u2, l2):
        terms.append((s, d + ((("G", len(up), up, lo),) if up else ())))
    for s, d, up, lo in bare_product(u2, l2, u1, l1):
        terms.append((-s, d + ((("G", len(up), up, lo),) if up else ())))
    return mkno.collect(terms)


D2 = ("p", "v"), ("q", "w")
ONE = commutator_expectation(("r",), ("s",), D2[0], D2[1])
TWO = commutator_expectation(("r", "t"), ("s", "u"), D2[0], D2[1])


def show(label, terms):
    print("%s : %d terms" % (label, len(terms)))
    for s, f in terms:
        parts = []
        for n, r, u, l in f:
            sym = {"G": "g", "D": "d"}.get(n, n)
            parts.append("%s%d[%s;%s]" % (sym, r, "".join(u), "".join(l)))
        print("   %+d  %s" % (s, " ".join(parts)))


show("<[a^r_s, a^pv_qw]>", ONE)
print()
by = {}
for s, f in TWO:
    g = [r for n, r, u, l in f if n == "G"]
    by[g[0] if g else 0] = by.get(g[0] if g else 0, 0) + 1
print("<[a^rt_su, a^pv_qw]> : %d terms  (%s)" %
      (len(TWO), ", ".join("gamma_%d x%d" % (k, v) for k, v in sorted(by.items()))))
print()

# ---------- numerical verification of both ----------
M = 6
rng = np.random.default_rng(31)
dim = 1 << M
popc = np.array([bin(s).count("1") for s in range(dim)])
psi = np.where(popc == 3, rng.normal(size=dim), 0.0)
F = Fock(M, psi / np.linalg.norm(psi))


def evaluate(terms, sub):
    tot = 0.0
    for c, facs in terms:
        val = float(c)
        for n, r, u, l in facs:
            iu = [sub[i] for i in u]
            il = [sub[i] for i in l]
            if n == "D":
                val *= 1.0 if iu[0] == il[0] else 0.0
            elif n == "G":
                val *= F.rdm(iu, il)
            if val == 0.0:
                break
        tot += val
    return tot


def brute(u1, l1, u2, l2, sub):
    def ap(word, v):
        for kind, idx in reversed(word):
            v = F.ad(idx, v) if kind == "c" else F.aa(idx, v)
        return v
    a = mkno.ops_of([sub[i] for i in u1], [sub[i] for i in l1])
    b = mkno.ops_of([sub[i] for i in u2], [sub[i] for i in l2])
    return float(F.psi @ (ap(a, ap(b, F.psi)) - ap(b, ap(a, F.psi))))


w1 = w2 = 0.0
for _ in range(200):
    sub = {k: int(rng.integers(M)) for k in "pqrstuvw"}
    w1 = max(w1, abs(evaluate(ONE, sub) - brute(("r",), ("s",), D2[0], D2[1], sub)))
    w2 = max(w2, abs(evaluate(TWO, sub) - brute(("r", "t"), ("s", "u"), D2[0], D2[1], sub)))
print("numerical check, M=%d N=3, 200 tuples:" % M)
print("  one-body : %.3e" % w1)
print("  two-body : %.3e" % w2)
print("  VERDICT:", "BOTH EXACT" if max(w1, w2) < 1e-12 else "MISMATCH")
