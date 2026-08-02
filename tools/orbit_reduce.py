#!/usr/bin/env python3
"""Reduce a term list to one signed representative per symmetry orbit, and
   write it back out in Wolfram term-list format so the engine can multiply it.

   The reduction is EXACT: the full list is Sum over the group of sgn(g)*g(rep).
   Deltas are UNPACKED first -- packed deltas are not canonical (see docs/10).

   Usage:
     python3 tools/orbit_reduce.py IN.m OUT.m  --swaps r:t,s:u
"""
import sys, itertools
sys.path.insert(0, 'tools')
from term_audit import parse, canon, relabel

def unpack(ts):
    out = []
    for n, r, up, lo in ts:
        if n == 'D' and r > 1:
            out += [('D', 1, (a,), (b,)) for a, b in zip(up, lo)]
        else:
            out.append((n, r, up, lo))
    return out

def signed_group(swaps):
    g = []
    for bits in itertools.product([0, 1], repeat=len(swaps)):
        m, sgn = {}, 1
        for b, k in zip(bits, swaps):
            if b: m.update(swaps[k]); sgn = -sgn
        g.append((m, sgn))
    return g

WNAME = {'G': '\\[Gamma]', 'L': '\\[Lambda]', 'H': '\\[Eta]', 'D': '\\[Delta]'}
def to_wl(sign, tensors):
    facs = [f'{{"Sign", {sign}}}']
    for n, r, up, lo in tensors:
        nm = WNAME.get(n, n)
        u = ", ".join(f'"{i}"' for i in up); l = ", ".join(f'"{i}"' for i in lo)
        facs.append(f'{{{{{nm}, {r}}}, {{{u}}}, {{{l}}}}}')
    return "{" + ", ".join(facs) + "}"

def reduce_list(path, swaps):
    terms = parse(path)
    base = {}
    for sg, ts in terms:
        s2, k = canon(sg, unpack(ts)); base[k] = base.get(k, 0) + s2
    base = {k: v for k, v in base.items() if v != 0}
    group = signed_group(swaps)
    seen, reps = set(), []
    for k in base:
        if k in seen: continue
        reps.append(k)
        for m, _ in group:
            _, k2 = canon(1, relabel(list(k), m)); seen.add(k2)
    return base, reps, group

def expand(base, reps, group):
    out = {}
    for k in reps:
        img = {}
        for m, sgn in group:
            s2, k2 = canon(sgn, relabel(list(k), m)); img.setdefault(k2, s2)
        for k2, s2 in img.items():
            out[k2] = out.get(k2, 0) + s2 * base[k]
    return {k: v for k, v in out.items() if v != 0}

if __name__ == "__main__":
    inp, outp = sys.argv[1], sys.argv[2]
    spec = sys.argv[sys.argv.index('--swaps') + 1]
    swaps = {}
    for pair in spec.split(','):
        a, b = pair.split(':'); swaps[f'{a}<->{b}'] = {a: b, b: a}
    base, reps, group = reduce_list(inp, swaps)
    assert expand(base, reps, group) == base, "orbit expansion is not exact!"
    with open(outp, 'w') as f:
        f.write("{" + ",\n ".join(to_wl(base[k], k) for k in reps) + "}\n")
    print(f"{inp}: {len(base)} terms -> {len(reps)} representatives "
          f"({len(base)/len(reps):.1f}x)   round trip verified exact")
