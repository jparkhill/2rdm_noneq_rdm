"""Particle-hole symmetry audit of an MK-derived term list.

The MKNO notebook reports terms after the substitution eta = delta - gamma,
so a term is a product of rank-1 factors drawn from {gamma, delta}.  Undoing
that substitution (delta -> eta + gamma) and expanding returns the list to the
gamma/eta representation in which particle-hole conjugation acts by the simple
relabelling gamma <-> eta.

For a degree-d term list, p-h symmetry implies the coefficient mass is
symmetric under (number of etas) k <-> d-k.  That is the weak test.  The strong
test is structural: the coefficient of every monomial must map to +/- the
coefficient of its conjugate.
"""
import json, collections, itertools

def load(path, cell):
    raw = json.load(open(path))[str(cell)]
    out = []
    for t in raw:
        fac = tuple((te['name'], te['up'][0], te['lo'][0]) for te in t['tensors'])
        out.append((t['sign'], fac))
    return out

def expand_delta(terms):
    """delta^a_b -> eta^a_b + gamma^a_b, expanded and collected."""
    coeff = collections.Counter()
    for sign, fac in terms:
        slots = [i for i, f in enumerate(fac) if f[0] == 'delta']
        for choice in itertools.product(('eta', 'gamma'), repeat=len(slots)):
            m = list(fac)
            for i, name in zip(slots, choice):
                m[i] = (name, m[i][1], m[i][2])
            coeff[tuple(sorted(m))] += sign
    return {k: v for k, v in coeff.items() if v != 0}

def conj(mono, transpose):
    swap = {'gamma': 'eta', 'eta': 'gamma'}
    out = [(swap[n], (b, a) if transpose else (a, b)) for n, a, b in mono]
    return tuple(sorted((n, i, j) for n, (i, j) in out))

def audit(coeff, label, transpose, phase):
    ok = bad = 0; worst = None
    for m, c in coeff.items():
        cm = conj(m, transpose)
        d = coeff.get(cm, 0)
        want = phase * c
        if d == want: ok += 1
        else:
            bad += 1
            if worst is None or abs(d - want) > abs(worst[1]): worst = (m, d - want)
    print(f"  {label:44s} match {ok:4d} / {ok+bad:4d}"
          + ("  <-- SYMMETRIC" if bad == 0 else f"   mismatch {bad}"))
    return bad == 0

if __name__ == '__main__':
    terms = load('docs/MKNO_terms.json', 1624)
    d = len(terms[0][1])
    coeff = expand_delta(terms)
    print(f"cell 1624: {len(terms)} terms, degree {d}"
          f" -> {len(coeff)} distinct gamma/eta monomials\n")

    mass = collections.Counter()
    for m, c in coeff.items():
        mass[sum(1 for n, _, _ in m if n == 'eta')] += c
    print("coefficient mass by eta-degree k (p-h symmetry needs k <-> %d-k):" % d)
    for k in range(d + 1):
        mirror = mass.get(d - k, 0)
        flag = '' if mass.get(k, 0) == mirror or k > d - k else \
               ('  == mirror' if mass.get(k,0) == mirror else '  != mirror %+d' % mirror)
        print(f"   k={k}: {mass.get(k,0):+6d}{flag}")
    print()
    print("structural test (coefficient of each monomial vs its conjugate):")
    for transpose in (False, True):
        for phase in (+1, -1):
            audit(coeff, f"gamma<->eta, transpose={transpose}, phase={phase:+d}",
                  transpose, phase)
