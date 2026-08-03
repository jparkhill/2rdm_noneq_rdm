#!/usr/bin/env python3
"""Compare the cumulant-free scalar part against MKNO cell 1624.

CONTINUITY records cell 1624 -- the 1-2-2 triple commutator with no projectors
-- as 96 terms: 64 of signature delta.delta.gamma.gamma.gamma and 32 of
delta.delta.delta.gamma.gamma. That is an independent, known-good fixture, so it
is a real test of the port rather than a self-consistency check.
"""
import sys
from collections import Counter
sys.path.insert(0, "tools")
from term_audit import parse, canon
from verify_stage_operator import unpack

d = {}
for sg, ts in parse("results/py_comm122_lam.m"):
    s, k = canon(sg, unpack(ts))
    d[k] = d.get(k, 0) + s
d = {k: v for k, v in d.items() if v}

scalar_nolam = [k for k in d
                if not any(n == "ta" for n, r, u, l in k)
                and not any(n == "L" for n, r, u, l in k)]

sig = Counter()
for k in scalar_nolam:
    names = sorted(n for n, r, u, l in k)
    sig[tuple(names)] += 1

print("cumulant-free scalar part: %d terms" % len(scalar_nolam))
print()
print("by signature:")
for s, c in sorted(sig.items(), key=lambda x: (-x[1], x[0])):
    pretty = ".".join("delta" if n == "D" else "gamma" if n == "G" else n for n in s)
    print("   %-40s %4d" % (pretty, c))
print()
print("CONTINUITY records for cell 1624: 96 total, 64 x delta.delta.gamma.gamma.gamma,")
print("                                            32 x delta.delta.delta.gamma.gamma")

want = {("D", "D", "G", "G", "G"): 64, ("D", "D", "D", "G", "G"): 32}
got = {k: v for k, v in sig.items()}
print()
print("MATCH:" if got == want else "DIFFERS:", "total %d" % len(scalar_nolam))
if got != want:
    for k in set(want) | set(got):
        if want.get(k) != got.get(k):
            print("   %-30s expected %s got %s" % (k, want.get(k), got.get(k)))
