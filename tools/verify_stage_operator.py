#!/usr/bin/env python3
"""Verify a derived MKNO stage as an OPERATOR IDENTITY, not just its expectation
   value.  Both sides are applied to a random Fock-space vector, so the residual
   normal-ordered operators {a_k} are tested too -- which is where an error would
   hide and then propagate into the next product.

   Usage: python3 tools/verify_stage_operator.py <stage.m> <lhs-spec> [ntrial]
     lhs-spec is a product of factors from  A={a^p_q}  B={a^rt_su}  C={a^vx_wy}
     e.g.  AB   BA   ABC   BAC   CBA   CAB
"""
import sys, itertools, numpy as np
sys.path.insert(0, 'tools')
from gno import Fock
from term_audit import parse, canon

def unpack(ts):
    out = []
    for n, r, up, lo in ts:
        if n == 'D' and r > 1:
            out += [('D', 1, (a,), (b,)) for a, b in zip(up, lo)]
        else:
            out.append((n, r, up, lo))
    return out

def load(path):
    d = {}
    for sg, ts in parse(path):
        s2, k = canon(sg, unpack(ts)); d[k] = d.get(k, 0) + s2
    return {k: v for k, v in d.items() if v != 0}

FAC = {'A': (['p'], ['q']), 'B': (['r', 't'], ['s', 'u']),
       'C': (['v', 'x'], ['w', 'y'])}

def main():
    path, spec = sys.argv[1], sys.argv[2]
    ntrial = int(sys.argv[3]) if len(sys.argv) > 3 else 12
    M = int(sys.argv[4]) if len(sys.argv) > 4 else 8
    N = M // 2
    rng = np.random.default_rng(20260802)
    dim = 1 << M
    popc = np.array([bin(s).count("1") for s in range(dim)])
    psi = np.where(popc == N, rng.normal(size=dim), 0.0)
    F = Fock(M, psi)
    terms = load(path)
    print(f"{path}: {len(terms)} terms;  LHS = {' . '.join(spec)};  "
          f"M={M} spin orbitals, N={N} electrons")

    used = sorted({i for k in terms for t in k for i in t[2] + t[3]})
    worst = 0.0; scale = 0.0
    for trial in range(ntrial):
        sub = {a: int(rng.integers(M)) for a in
               ['p', 'q', 'r', 't', 's', 'u', 'v', 'x', 'w', 'y']}
        v = rng.normal(size=dim)
        # LHS: apply factors right to left
        lhs = v.copy()
        for ch in spec[::-1]:
            up, lo = FAC[ch]
            lhs = F.gno([sub[i] for i in up], [sub[i] for i in lo], lhs)
        # RHS: the derived term list
        rhs = np.zeros_like(v)
        for k, c in terms.items():
            coef = float(c); op = None
            for n, r, up, lo in k:
                idx_u = [sub[i] for i in up]; idx_l = [sub[i] for i in lo]
                if n == 'D':
                    coef *= 1.0 if idx_u[0] == idx_l[0] else 0.0
                elif n == 'G':
                    coef *= F.rdm(idx_u, idx_l)
                elif n == 'L':
                    coef *= F.lam(idx_u, idx_l)
                elif n == 'ta':
                    op = (idx_u, idx_l)
                if coef == 0.0: break
            if coef == 0.0: continue
            rhs = rhs + coef * (F.gno(op[0], op[1], v, tag=trial) if op else v)
        worst = max(worst, np.abs(lhs - rhs).max())
        scale = max(scale, np.abs(lhs).max())
    print(f"  {ntrial} random index tuples + random vectors:")
    print(f"  max |LHS - RHS| = {worst:.3e}   (typical |LHS| = {scale:.3f})")
    print("  VERDICT:", "EXACT OPERATOR IDENTITY" if worst < 1e-9 * max(scale, 1)
          else "MISMATCH")

if __name__ == "__main__":
    main()
