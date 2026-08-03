#!/usr/bin/env python3
"""Verify a derived term list against the double commutator [[A,B],C].

LHS is built from explicit Fock-space operators: ABC - BAC - CAB + CBA applied
to a random vector. RHS is the term list. Same oracle as
tools/verify_stage_operator.py, different left-hand side.

Usage: python3 tools/verify_commutator.py <terms.m> [ntrial] [M]
"""
import sys
import numpy as np
sys.path.insert(0, 'tools')
from gno import Fock
from term_audit import parse, canon
from verify_stage_operator import unpack

FAC = {'A': (['p'], ['q']), 'B': (['r', 't'], ['s', 'u']),
       'C': (['v', 'x'], ['w', 'y'])}
WORDS = [("ABC", 1), ("BAC", -1), ("CAB", -1), ("CBA", 1)]


def load(path):
    d = {}
    for sg, ts in parse(path):
        s2, k = canon(sg, unpack(ts))
        d[k] = d.get(k, 0) + s2
    return {k: v for k, v in d.items() if v != 0}


def main():
    path = sys.argv[1]
    ntrial = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    M = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    N = M // 2
    rng = np.random.default_rng(20260802)
    dim = 1 << M
    popc = np.array([bin(s).count("1") for s in range(dim)])
    psi = np.where(popc == N, rng.normal(size=dim), 0.0)
    F = Fock(M, psi)
    terms = load(path)
    print(f"{path}: {len(terms)} terms;  LHS = [[A,B],C];  M={M}, N={N}")

    worst = scale = 0.0
    for trial in range(ntrial):
        sub = {a: int(rng.integers(M)) for a in
               ['p', 'q', 'r', 't', 's', 'u', 'v', 'x', 'w', 'y']}
        v = rng.normal(size=dim)

        lhs = np.zeros_like(v)
        for word, coef in WORDS:
            acc = v.copy()
            for ch in word[::-1]:
                up, lo = FAC[ch]
                acc = F.gno([sub[i] for i in up], [sub[i] for i in lo], acc)
            lhs = lhs + coef * acc

        rhs = np.zeros_like(v)
        for k, c in terms.items():
            coef = float(c)
            op = None
            for n, r, up, lo in k:
                iu = [sub[i] for i in up]
                il = [sub[i] for i in lo]
                if n == 'D':
                    coef *= 1.0 if iu[0] == il[0] else 0.0
                elif n == 'G':
                    coef *= F.rdm(iu, il)
                elif n == 'L':
                    coef *= F.lam(iu, il)
                elif n == 'ta':
                    op = (iu, il)
                if coef == 0.0:
                    break
            if coef == 0.0:
                continue
            rhs = rhs + coef * (F.gno(op[0], op[1], v, tag=trial) if op else v)

        worst = max(worst, np.abs(lhs - rhs).max())
        scale = max(scale, np.abs(lhs).max())

    print(f"  {ntrial} random index tuples + random vectors:")
    print(f"  max |LHS - RHS| = {worst:.3e}   (typical |LHS| = {scale:.3f})")
    print("  VERDICT:", "EXACT OPERATOR IDENTITY"
          if worst < 1e-9 * max(scale, 1) else "MISMATCH")


if __name__ == "__main__":
    main()
