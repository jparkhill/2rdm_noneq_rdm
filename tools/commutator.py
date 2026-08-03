#!/usr/bin/env python3
"""Build [[A,B],C] = ABC - BAC - CAB + CBA with the validated symbolic engine.

Each of the four products is written out and verified independently as an
operator identity; the commutator is their signed sum, so linearity carries the
verification, but tools/verify_commutator.py checks the combination directly as
well in case the signs of the combination itself are wrong.
"""
import sys, time
sys.path.insert(0, ".")
from qdiss.symbolic import mkno

A = (("p",), ("q",))
B = (("r", "t"), ("s", "u"))
C = (("v", "x"), ("w", "y"))


def times_right(terms, b):
    """terms * {a^b}"""
    out = []
    for sign, facs in terms:
        gno = [f for f in facs if f[0] == mkno.GNO]
        rest = tuple(f for f in facs if f[0] != mkno.GNO)
        if not gno:
            out.append((sign, rest + ((mkno.GNO, len(b[0]), b[0], b[1]),)))
            continue
        _, _, up, lo = gno[0]
        for s2, f2 in mkno.gno_product((up, lo), b):
            out.append((sign * s2, rest + f2))
    return mkno.collect(out)


def times_left(b, terms):
    """{a^b} * terms"""
    out = []
    for sign, facs in terms:
        gno = [f for f in facs if f[0] == mkno.GNO]
        rest = tuple(f for f in facs if f[0] != mkno.GNO)
        if not gno:
            out.append((sign, rest + ((mkno.GNO, len(b[0]), b[0], b[1]),)))
            continue
        _, _, up, lo = gno[0]
        for s2, f2 in mkno.gno_product(b, (up, lo)):
            out.append((sign * s2, rest + f2))
    return mkno.collect(out)


def scale(terms, c):
    return [(s * c, f) for s, f in terms]


t0 = time.time()
ab = mkno.gno_product(A, B)
ba = mkno.gno_product(B, A)
print(f"AB {len(ab):5d}   BA {len(ba):5d}")

prods = {
    "ABC": times_right(ab, C),
    "BAC": times_right(ba, C),
    "CAB": times_left(C, ab),
    "CBA": times_left(C, ba),
}
for k, v in prods.items():
    print(f"{k} {len(v):6d} terms")
    mkno.to_wolfram(v, f"results/py_{k}.m")

comm = mkno.collect(
    scale(prods["ABC"], 1) + scale(prods["BAC"], -1)
    + scale(prods["CAB"], -1) + scale(prods["CBA"], 1)
)
print(f"\n[[A,B],C] gamma form : {len(comm)} terms   ({time.time()-t0:.1f}s)")
mkno.to_wolfram(comm, "results/py_comm122.m")

lam = mkno.to_cumulant_form(comm)
nl = sum(1 for s, f in lam if any(x[0] == "L" for x in f))
print(f"[[A,B],C] cumulant   : {len(lam)} terms, {nl} carrying lambda")
mkno.to_wolfram(lam, "results/py_comm122_lam.m")

ranks = {}
for s, facs in lam:
    for n, r, u, l in facs:
        ranks.setdefault(n, set()).add(r)
print("factor ranks:", {k: sorted(v) for k, v in sorted(ranks.items())})
