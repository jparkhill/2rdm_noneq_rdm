#!/usr/bin/env python3
"""Decompose the derived commutator by operator rank and cumulant content.

Two questions this answers:
  * which terms survive taking an expectation value (only the pure scalars,
    since <{a_k}> = 0 for k >= 1)
  * what the closure has to kill for a gamma_2 theory to be closed
"""
import sys
from collections import Counter
sys.path.insert(0, "tools")
from term_audit import parse, canon
from verify_stage_operator import unpack

path = sys.argv[1] if len(sys.argv) > 1 else "results/py_comm122_lam.m"
d = {}
for sg, ts in parse(path):
    s, k = canon(sg, unpack(ts))
    d[k] = d.get(k, 0) + s
d = {k: v for k, v in d.items() if v}

by_op = Counter()
by_lam = Counter()
for k in d:
    o = [r for n, r, u, l in k if n == "ta"]
    L = [r for n, r, u, l in k if n == "L"]
    by_op[o[0] if o else 0] += 1
    by_lam[tuple(sorted(L))] += 1

print(f"{path}: {len(d)} terms")
print()
print("by operator rank (0 = pure scalar, the part that survives <.>):")
for r in sorted(by_op):
    label = "scalar" if r == 0 else "{a_%d}" % r
    print("   %-8s %5d" % (label, by_op[r]))
print()
print("by cumulant content:")
for L in sorted(by_lam, key=lambda x: (len(x), x)):
    label = "none" if not L else "lambda" + "".join("_%d" % r for r in L)
    print("   %-16s %5d" % (label, by_lam[L]))
print()
keep = sum(v for L, v in by_lam.items() if 3 not in L)
drop = sum(v for L, v in by_lam.items() if 3 in L)
print("closure lambda_3 = 0:  %d terms survive, %d dropped" % (keep, drop))

scalars = {k: v for k, v in d.items() if not any(n == "ta" for n, r, u, l in k)}
print()
print("pure-scalar part (the closed EOM contribution): %d terms" % len(scalars))
sl = Counter()
for k in scalars:
    L = [r for n, r, u, l in k if n == "L"]
    sl[tuple(sorted(L))] += 1
for L in sorted(sl, key=lambda x: (len(x), x)):
    label = "none" if not L else "lambda" + "".join("_%d" % r for r in L)
    print("   %-16s %5d" % (label, sl[L]))
