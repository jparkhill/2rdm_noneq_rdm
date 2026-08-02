#!/usr/bin/env python3
"""Which terms of a derived stage are wrong?

The stage is an operator identity  LHS = sum_i c_i * (scalars_i) * {a_i}.
Apply both sides to random Fock vectors and look at the residual

    D(v) = LHS(v) - RHS(v).

If the derivation is wrong only in the coefficients of a few residual
operators, D lies in the span of those operators applied to v.  We fit D
against the operator basis grouped by residual rank and report which rank
carries the error and what coefficient correction it needs.
"""
import sys, numpy as np
sys.path.insert(0, 'tools')
from gno import Fock
from term_audit import parse

FAC = {'A': (['p'], ['q']), 'B': (['r', 't'], ['s', 'u']),
       'C': (['v', 'x'], ['w', 'y'])}
IDX = ['p', 'q', 'r', 't', 's', 'u', 'v', 'x', 'w', 'y']

def main(path, spec, M=8, ntrial=4, seed=20260802):
    N = M // 2
    rng = np.random.default_rng(seed)
    dim = 1 << M
    popc = np.array([bin(s).count("1") for s in range(dim)])
    psi = np.where(popc == N, rng.normal(size=dim), 0.0)
    F = Fock(M, psi)
    terms = list(parse(path))
    print(f"{path}: {len(terms)} terms;  LHS = {' . '.join(spec)};  M={M}")

    for trial in range(ntrial):
        sub = {a: int(rng.integers(M)) for a in IDX}
        v = rng.normal(size=dim)
        lhs = v.copy()
        for ch in spec[::-1]:
            up, lo = FAC[ch]
            lhs = F.gno([sub[i] for i in up], [sub[i] for i in lo], lhs)

        # RHS, accumulated separately per residual rank
        byrank = {}
        opvec = {}
        for sg, ts in terms:
            coef = float(sg); op = None
            for n, r, up, lo in ts:
                iu = [sub[i] for i in up]; il = [sub[i] for i in lo]
                if n == 'D':
                    coef *= 1.0 if all(a == b for a, b in zip(iu, il)) else 0.0
                elif n == 'G':
                    coef *= F.rdm(iu, il)
                elif n == 'L':
                    coef *= F.lam(iu, il)
                elif n == 'ta':
                    op = (r, tuple(iu), tuple(il))
                if coef == 0.0: break
            if coef == 0.0: continue
            rk = op[0] if op else 0
            if op:
                key = op
                if key not in opvec:
                    opvec[key] = F.gno(list(op[1]), list(op[2]), v, tag=trial)
                contrib = coef * opvec[key]
            else:
                contrib = coef * v
            byrank[rk] = byrank.get(rk, 0.0) + contrib

        rhs = sum(byrank.values())
        D = lhs - rhs
        line = "  ".join(f"r{k}:{np.linalg.norm(x):.3f}" for k, x in sorted(byrank.items()))
        print(f"\n trial {trial}: |LHS|={np.linalg.norm(lhs):.4f}  "
              f"|D|={np.linalg.norm(D):.4f}")
        print(f"   RHS pieces by residual rank: {line}")

        # fit D against the top-rank operator vectors
        for rk in (5, 4):
            cols = [opvec[k] for k in opvec if k[0] == rk]
            if not cols: continue
            Amat = np.array(cols).T
            c, *_ = np.linalg.lstsq(Amat, D, rcond=None)
            resid = np.linalg.norm(D - Amat @ c)
            print(f"   fit D with the {len(cols)} rank-{rk} operators: "
                  f"residual {resid:.3e}   coefficients {np.round(c, 6)[:6]}")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2],
         M=int(sys.argv[3]) if len(sys.argv) > 3 else 8,
         ntrial=int(sys.argv[4]) if len(sys.argv) > 4 else 4)
