#!/usr/bin/env python3
"""Derive A.B and A.B.C symbolically and write them for the numerical verifier.

The triple product is where the Wolfram pipeline goes wrong (docs/12): its
rv1v2 stage disagrees with brute force by the size of the answer. This builds
the same object from the definition instead.
"""
import sys, time
sys.path.insert(0, ".")
from qdiss.symbolic import mkno

A = (("p",), ("q",))
B = (("r", "t"), ("s", "u"))
C = (("v", "x"), ("w", "y"))


def times(terms, b):
    """Multiply a collected term list on the right by a GNO operator b."""
    out = []
    for sign, facs in terms:
        gno = [f for f in facs if f[0] == mkno.GNO]
        rest = tuple(f for f in facs if f[0] != mkno.GNO)
        if not gno:
            out.append((sign, rest + ((mkno.GNO, len(b[0]), b[0], b[1]),)))
            continue
        name, rank, up, lo = gno[0]
        for s2, f2 in mkno.gno_product((up, lo), b):
            out.append((sign * s2, rest + f2))
    return mkno.collect(out)


t0 = time.time()
ab = mkno.gno_product(A, B)
print(f"A.B      : {len(ab):6d} terms   ({time.time()-t0:.1f}s)")

t0 = time.time()
abc = times(ab, C)
print(f"A.B.C    : {len(abc):6d} terms   ({time.time()-t0:.1f}s)")

mkno.to_wolfram(ab, "results/py_AB.m")
mkno.to_wolfram(abc, "results/py_ABC.m")

ranks = {}
for s, facs in abc:
    for n, r, u, l in facs:
        ranks.setdefault(n, set()).add(r)
print("factor ranks present:", {k: sorted(v) for k, v in sorted(ranks.items())})
