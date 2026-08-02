"""Do the Wolfram combine and the pure-Python cancellation agree TERM BY TERM?"""
import sys
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
        s2, k = canon(sg, unpack(ts)); d[k] = d.get(k, 0) + s2
    return {k: v for k, v in d.items() if v != 0}

wl  = load('results/mkno_full_cumulant.m')          # Wolfram AssociativeSimplify
py  = load('results/stage_joined_raw.m')            # Python signed cancellation
print(f"Wolfram combine : {len(wl)} terms")
print(f"Python cancel   : {len(py)} terms")
missing = [k for k in py if k not in wl]
extra   = [k for k in wl if k not in py]
wrong   = [k for k in py if k in wl and wl[k] != py[k]]
print(f"IDENTICAL: {wl == py}   (missing {len(missing)}, extra {len(extra)}, wrong-coeff {len(wrong)})")

wln = load('results/mkno_full_nocum.m')
pyn = {k: v for k, v in py.items() if not any(t[0] == 'L' for t in k)}
print(f"nocum  Wolfram {len(wln)}  Python {len(pyn)}  identical: {wln == pyn}")
