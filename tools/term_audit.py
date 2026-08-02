#!/usr/bin/env python3
"""Audit a Wolfram term-list checkpoint (results/stage_*.m) for redundancy,
   WITHOUT occupying a Wolfram kernel.

   Answers three separate questions that are easy to conflate:

     1. like-term merging  -- are any two terms literally the same monomial?
                              (this is all AssociativeSimplify/Simplify can do)
     2. orbit factorization-- is the list a full SIGNED ORBIT under the index
                              antisymmetries the operators already have?  If so
                              only a fundamental domain is independent content.
     3. canonical packing  -- do equivalent terms end up in different
                              representations (delta packing) so they no longer
                              match?

   Conventions verified against MKNO.nb cell 145:
     {{d,2},{p,r},{q,s}}  ==  delta^p_q delta^r_s   -- a PLAIN PRODUCT, so it is
     canonicalized by sorting the (upper,lower) PAIRS, not by independent
     antisymmetry of the two index groups.  (Treating it as antisymmetric
     wrongly merges 159 distinct terms in stage_rv1v2.)

   Usage: python3 tools/term_audit.py results/stage_rv1v2.m [--swaps r:t,s:u,v:x,w:y]
"""
import re, sys, itertools
from collections import Counter

ANTISYM = {'G', 'L', 'ta', 'H'}      # gamma_k, lambda_k, {a_k}, eta_k


def parse(path):
    s = open(path).read()
    s = (s.replace('\\[Gamma]', 'G').replace('\\[Lambda]', 'L')
          .replace('\\[Eta]', 'H').replace('\\[Delta]', 'D'))
    s = re.sub(r'\s+', '', s)
    def split_top(body):
        out, depth, cur = [], 0, []
        for ch in body:
            if ch == '{': depth += 1
            if ch == '}': depth -= 1
            if ch == ',' and depth == 0:
                out.append(''.join(cur)); cur = []
            else:
                cur.append(ch)
        if cur: out.append(''.join(cur))
        return out
    terms = []
    for t in split_top(s[1:-1]):
        sign, tensors = 1, []
        for f in split_top(t[1:-1]):
            if f.startswith('{"Sign"'):
                sign = int(f.split(',')[1].rstrip('}'))
                continue
            m = re.match(r'^\{\{(\w+),(\d+)\},\{(.*?)\},\{(.*?)\}\}$', f)
            if not m:
                tensors.append(('??', 0, (f,), ())); continue
            tensors.append((m.group(1), int(m.group(2)),
                            tuple(x.strip('"') for x in m.group(3).split(',') if x),
                            tuple(x.strip('"') for x in m.group(4).split(',') if x)))
        terms.append((sign, tensors))
    return terms


def perm_sign(seq):
    seq = list(seq); sgn = 1
    for i in range(len(seq)):
        for j in range(len(seq) - 1 - i):
            if seq[j] > seq[j + 1]:
                seq[j], seq[j + 1] = seq[j + 1], seq[j]; sgn = -sgn
    return sgn, tuple(seq)


def canon(sign, tensors):
    s, cts = sign, []
    for t in tensors:
        n, r, up, lo = t
        if n == 'D' and r > 1:                       # packed product of deltas
            pr = sorted(zip(up, lo))
            up = tuple(a for a, _ in pr); lo = tuple(b for _, b in pr)
        elif n in ANTISYM and r > 1:
            su, up = perm_sign(up); sl, lo = perm_sign(lo); s *= su * sl
        cts.append((n, r, up, lo))
    cts.sort()
    return s, tuple(cts)


def relabel(ts, m):
    return [(n, r, tuple(m.get(i, i) for i in up), tuple(m.get(i, i) for i in lo))
            for n, r, up, lo in ts]


def audit(path, swaps):
    terms = parse(path)
    base = {}
    for sg, ts in terms:
        s2, k = canon(sg, ts); base[k] = base.get(k, 0) + s2
    base = {k: v for k, v in base.items() if v != 0}
    print(f"file              : {path}")
    print(f"terms as stored   : {len(terms)}")
    print(f"after like-merging: {len(base)}   "
          f"(gain {len(terms) - len(base)} terms = "
          f"{100*(1-len(base)/max(len(terms),1)):.1f}%)")
    print(f"terms with lambda : {sum(1 for k in base if any(t[0]=='L' for t in k))}")

    viol = {}
    for name, m in swaps.items():
        img = {}
        for k, v in base.items():
            s2, k2 = canon(v, relabel(list(k), m)); img[k2] = img.get(k2, 0) + s2
        bad = [k for k, v in base.items() if img.get(k, 0) != -v]
        viol[name] = bad
        print(f"  antisym {name:8s}: "
              f"{'OK' if not bad else str(len(bad)) + ' violations'}")

    group = []
    for bits in itertools.product([0, 1], repeat=len(swaps)):
        m = {}
        for b, kk in zip(bits, swaps):
            if b: m.update(swaps[kk])
        group.append(m)
    seen, orb = set(), 0
    reps = []
    for k in base:
        if k in seen: continue
        orb += 1; reps.append(k)
        for m in group:
            _, k2 = canon(1, relabel(list(k), m)); seen.add(k2)
    print(f"independent orbits: {orb}   "
          f"({len(base)/max(orb,1):.1f}x smaller than the stored list)")

    allbad = set().union(*viol.values()) if viol else set()
    if allbad:
        shapes = Counter(tuple(sorted((t[0], t[1]) for t in k)) for k in allbad)
        print("\nviolating-term shapes (all contain a packed delta -> "
              "non-canonical delta merge):")
        for s_, c in shapes.most_common():
            print("   ", c, " x ", " * ".join(f"{n}{r}" for n, r in s_))
    return base, reps


if __name__ == "__main__":
    path = sys.argv[1]
    sw = {'r<->t': {'r':'t','t':'r'}, 's<->u': {'s':'u','u':'s'},
          'v<->x': {'v':'x','x':'v'}, 'w<->y': {'w':'y','y':'w'}}
    audit(path, sw)
