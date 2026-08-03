#!/usr/bin/env python3
"""The gamma_2 EOM's three-body term, and what lambda_3 = 0 injects into it.

Hypothesis under test (John): closing on the CUMULANT rather than on the DENSITY
destroys an exact cancellation. At N = 2 the BBGKY term carrying gamma_3
vanishes identically, so the 2-RDM EOM closes exactly. Rewriting gamma_3 as
separable(gamma_1, lambda_2) + lambda_3 and then setting lambda_3 = 0 leaves the
separable piece behind, injecting a three-body contribution where the truth is
zero. If so, instability observed at small N is an artifact of the
parameterization, not a failure of the physics.

Derived with the machinery already verified in tools/eom1.py, then checked
numerically before anything is concluded.
"""
import itertools
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from qdiss.symbolic import mkno
from gno import Fock
import gno


def bare_product(u1, l1, u2, l2):
    word = mkno.ops_of(u1, l1) + mkno.ops_of(u2, l2)
    out = []
    for s, deltas, cre, ann in mkno.normal_order(word):
        if len(cre) != len(ann):
            continue
        out.append((s, deltas, tuple(cre), tuple(reversed(ann))))
    return out


def commutator_expectation(u1, l1, u2, l2):
    terms = []
    for s, d, up, lo in bare_product(u1, l1, u2, l2):
        terms.append((s, d + ((("G", len(up), up, lo),) if up else ())))
    for s, d, up, lo in bare_product(u2, l2, u1, l1):
        terms.append((-s, d + ((("G", len(up), up, lo),) if up else ())))
    return mkno.collect(terms)


# <[ a^rt_su , a^pv_qw ]> : two-body H piece against the 2-RDM operator
H2 = ("r", "t"), ("s", "u")
D2 = ("p", "v"), ("q", "w")
terms = commutator_expectation(H2[0], H2[1], D2[0], D2[1])

by_rank = {}
for s, f in terms:
    g = [r for n, r, u, l in f if n == "G"]
    by_rank.setdefault(g[0] if g else 0, []).append((s, f))
print("<[a^rt_su, a^pv_qw]> : %d terms" % len(terms))
for r in sorted(by_rank):
    print("   gamma_%d terms: %d" % (r, len(by_rank[r])))
print()

# --- numerical check of the whole commutator ---
M = 6
rng = np.random.default_rng(23)
dim = 1 << M
popc = np.array([bin(s).count("1") for s in range(dim)])


def make(N):
    psi = np.where(popc == N, rng.normal(size=dim), 0.0)
    return Fock(M, psi / np.linalg.norm(psi))


def evaluate(terms, F, sub, gamma3=None):
    tot = 0.0
    for c, facs in terms:
        val = float(c)
        for n, r, u, l in facs:
            iu = [sub[i] for i in u]
            il = [sub[i] for i in l]
            if n == "D":
                val *= 1.0 if iu[0] == il[0] else 0.0
            elif n == "G":
                val *= (gamma3(iu, il) if (r == 3 and gamma3) else F.rdm(iu, il))
            if val == 0.0:
                break
        tot += val
    return tot


def brute(F, u1, l1, u2, l2, sub):
    def ap(word, v):
        for kind, idx in reversed(word):
            v = F.ad(idx, v) if kind == "c" else F.aa(idx, v)
        return v
    a = mkno.ops_of([sub[i] for i in u1], [sub[i] for i in l1])
    b = mkno.ops_of([sub[i] for i in u2], [sub[i] for i in l2])
    return float(F.psi @ (ap(a, ap(b, F.psi)) - ap(b, ap(a, F.psi))))


F3 = make(3)
worst = 0.0
labels = ["p", "q", "r", "s", "t", "u", "v", "w"]
for _ in range(200):
    sub = {k: int(rng.integers(M)) for k in labels}
    worst = max(worst, abs(evaluate(terms, F3, sub)
                           - brute(F3, H2[0], H2[1], D2[0], D2[1], sub)))
print("commutator verified numerically (N=3, 200 tuples): max error %.3e" % worst)
print("VERDICT:", "EXACT" if worst < 1e-12 else "MISMATCH")
print()

# --- the injected term ---
def separable_g3(F):
    """gamma_3 as reconstructed with lambda_3 = 0."""
    def f(up, lo):
        tot = 0.0
        for part in gno._set_partitions(3):
            if len(part) == 1:
                continue
            tot += gno._partition_value(F, up, lo, part)
        return tot
    return f


for N in (2, 3):
    F = make(N)
    sep = separable_g3(F)
    exact_c = spur_c = scale = 0.0
    for _ in range(400):
        sub = {k: int(rng.integers(M)) for k in labels}
        e = evaluate(terms, F, sub)                       # true gamma_3
        s = evaluate(terms, F, sub, gamma3=sep)           # lambda_3 = 0
        exact_c = max(exact_c, abs(e))
        spur_c = max(spur_c, abs(s - e))
        scale = max(scale, abs(e))
    print("N = %d" % N)
    print("   max |exact  d<gamma_2>/dt term|        = %.4e" % exact_c)
    print("   max |error injected by lambda_3 = 0|   = %.4e" % spur_c)
    if scale > 0:
        print("   relative                               = %.1f%%" % (100 * spur_c / scale))
    print()
