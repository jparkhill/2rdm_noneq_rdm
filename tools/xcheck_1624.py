#!/usr/bin/env python3
"""Does the new full-cumulant derivation reproduce the notebook's stored
   96-term cell-1624 result?

   Cell 1624 was built with AllStepsOfProductRule1PV / OneParticleVaccuum1/2,
   which are undefined anywhere in the notebook or the refs (docs/09) -- so it
   was treated as unreproducible and the new derivation as not term-comparable.

   Claim tested here: cell 1624 is exactly the SCALAR, CUMULANT-FREE slice of
   the general product rule, i.e. what you get by keeping only terms with no
   residual normal-ordered operator and no lambda.  If so the missing functions
   were doing nothing the general rule cannot do, and the blocker dissolves.
"""
import sys, json
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

base = {}
for sg, ts in parse('results/stage_joined_raw.m'):
    s2, k = canon(sg, unpack(ts)); base[k] = base.get(k, 0) + s2
ours = {k: v for k, v in base.items()
        if v != 0 and not any(t[0] in ('ta', 'L') for t in k)}

NAME = {'gamma': 'G', 'delta': 'D', 'lambda': 'L', 'eta': 'H'}
theirs = {}
for t in json.load(open('docs/MKNO_terms.json'))['1624']:
    ts = [(NAME[x['name']], len(x['up']), tuple(x['up']), tuple(x['lo']))
          for x in t['tensors']]
    s2, k = canon(t['sign'], unpack(ts)); theirs[k] = theirs.get(k, 0) + s2
theirs = {k: v for k, v in theirs.items() if v != 0}

print(f"new derivation, scalar cumulant-free slice : {len(ours)} terms")
print(f"notebook cell 1624 as stored               : {len(theirs)} terms")
print(f"IDENTICAL (terms and signs)                : {ours == theirs}")
assert ours == theirs
