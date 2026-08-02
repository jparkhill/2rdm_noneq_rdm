#!/usr/bin/env python3
"""Stage 1 of the spectral attribution: classify every term of the full-cumulant
   1-2-2 commutator by (cumulant content, residual-operator rank).

   The residual normal-ordered operator {a_k} decides WHICH equation a term
   belongs to: taking the expectation value of the commutator against the
   reference kills every {a_k} with k>=1, so
       ta0  -> scalar   (energy-like)
       ta1  -> 1-body   channel
       ta2  -> 2-body   channel   <-- the 2-RDM equation
       ta3+ -> higher-body channel (closure territory)
   Cumulant content decides whether the term is a SOURCE (no lambda: drives the
   cumulant from the mean-field part) or a RESPONSE (carries lambda: the
   cumulant feeding back on its own rate).  The NoCumulants rewrite deleted the
   entire response set, so in the old generator d(lambda)/dt had no dependence
   on lambda at all.
"""
import sys, collections
sys.path.insert(0, 'tools')
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
        s2, k = canon(sg, unpack(ts))
        d[k] = d.get(k, 0) + s2
    return {k: v for k, v in d.items() if v != 0}


def profile(term):
    """(residual operator rank, sorted tuple of cumulant ranks, #gamma, #delta)"""
    ta = [t[1] for t in term if t[0] == 'ta']
    lam = sorted(t[1] for t in term if t[0] == 'L')
    ng = sum(1 for t in term if t[0] == 'G')
    nd = sum(1 for t in term if t[0] == 'D')
    assert len(ta) <= 1, f"more than one residual operator: {term}"
    return (ta[0] if ta else 0), tuple(lam), ng, nd


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else 'results/stage_joined_raw.m'
    base = load(path)
    print(f"{path}: {len(base)} surviving terms\n")

    tab = collections.Counter()
    for k in base:
        ta, lam, ng, nd = profile(k)
        tab[(ta, lam)] += 1

    ranks = sorted({t for t, _ in tab})
    lams = sorted({l for _, l in tab}, key=lambda x: (len(x), x))
    hdr = f"{'residual':>10} | " + " | ".join(f"{('none' if not l else '+'.join('L'+str(i) for i in l)):>8}" for l in lams) + " |    total"
    print(hdr); print("-" * len(hdr))
    for t in ranks:
        row = [tab.get((t, l), 0) for l in lams]
        print(f"{'ta'+str(t) if t else 'scalar':>10} | " +
              " | ".join(f"{v:>8}" for v in row) + f" | {sum(row):>8}")
    print("-" * len(hdr))
    tot = [sum(tab.get((t, l), 0) for t in ranks) for l in lams]
    print(f"{'total':>10} | " + " | ".join(f"{v:>8}" for v in tot) + f" | {sum(tot):>8}")
