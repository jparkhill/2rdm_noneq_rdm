"""Particle-hole covariance audit of the SECULAR term list that became C++.

Source: the //NN{{Sign,...}} provenance comments in
refs/TightBindingTCL2/TCL_EE2.h (RhoDotSec, groups 61-64; RhoDotSecFull,
groups 71-74).  Every term has signature  Vs.Vt.gamma.gamma.eta.eta.

Particle-hole conjugation candidates tested:
  gamma <-> eta                                    (always)
  optional transpose of the gamma/eta index pair
  optional transpose of V upper/lower index pairs
  optional extra sign on Vs   -- physically motivated: Vs carries 1/(i*Delta)
                                 and Delta flips sign under p-h conjugation
  overall phase +-1

Terms are compared modulo relabelling of the summed (dummy) indices.
"""
import sys, itertools, collections
sys.path.insert(0, '/root/workspace/tools')
from cppterms import extract_file, fmt

FREE = ('p', 'q')


def norm_factor(n, u, l):
    """Quotient by the integral's own index symmetry.

    In the RI form used by the code, V^{pq}_{rs} = sum_R B(p,r) B(q,s) [/ i*Delta],
    and Delta = d2(p,r)+d2(q,s), so V is invariant under the simultaneous pair
    exchange (p,r) <-> (q,s), i.e. V^{pq}_{rs} == V^{qp}_{sr}.
    """
    if len(u) == 2:
        alt = ((u[1], u[0]), (l[1], l[0]))
        u, l = min(((u, l), alt))
    return (n, u, l)


def canon(facs, dummies):
    """Lexicographically smallest form over permutations of dummy labels."""
    best = None
    for perm in itertools.permutations(dummies):
        m = dict(zip(dummies, perm))
        f = tuple(sorted(
            norm_factor(n, tuple(m.get(i, i) for i in u),
                        tuple(m.get(i, i) for i in l))
            for n, u, l in facs))
        if best is None or f < best:
            best = f
    return best


def conjugate(facs, tr_g, tr_v, swap_pq=False):
    swap = {'gamma': 'eta', 'eta': 'gamma'}
    pq = {'p': 'q', 'q': 'p'} if swap_pq else {}
    rel = lambda t: tuple(pq.get(i, i) for i in t)
    out = []
    for n, u, l in facs:
        u, l = rel(u), rel(l)
        if n in swap:
            nu, nl = (l, u) if tr_g else (u, l)
            out.append((swap[n], nu, nl))
        else:
            nu, nl = (l, u) if tr_v else (u, l)
            out.append((n, nu, nl))
    return tuple(out)


def audit(terms, dummies, label):
    """terms: Counter of (canonical monomial) -> summed signed coefficient"""
    print(f"\n=== {label}: {sum(abs(v) for v in terms.values())} signed units, "
          f"{len(terms)} distinct monomials")
    rows = []
    for tr_g in (False, True):
      for tr_v in (False, True):
        for spq in (False, True):
            for vs_sign in (+1, -1):
                for phase in (+1, -1):
                    ok = 0
                    for m, c in terms.items():
                        cm = conjugate(m, tr_g, tr_v, spq)
                        nvs = sum(1 for n, _, _ in cm if n == 'Vs')
                        cm = canon(cm, dummies)
                        want = phase * (vs_sign ** nvs) * c
                        if terms.get(cm, 0) == want:
                            ok += 1
                    rows.append((ok, tr_g, tr_v, spq, vs_sign, phase))
    rows.sort(reverse=True)
    for ok, tg, tv, spq, vs, ph in rows[:5]:
        tag = '  <-- COVARIANT' if ok == len(terms) else ''
        print(f"   tr_gamma={int(tg)} tr_V={int(tv)} swap_pq={int(spq)} "
              f"sign(Vs)={vs:+d} phase={ph:+d}   match {ok:3d}/{len(terms)}{tag}")
    return rows[0][0] == len(terms)


def selftest(dummies):
    """A deliberately p-h covariant 4-term list must score 100%."""
    base = [(1, (('Vs', ('p', 'u'), ('p', 'r')), ('Vt', ('p', 's'), ('p', 't')),
                 ('gamma', ('p',), ('p',)), ('gamma', ('u',), ('t',)),
                 ('eta', ('p',), ('q',)), ('eta', ('s',), ('r',)))),
            (-1, (('Vs', ('q', 'u'), ('q', 'r')), ('Vt', ('q', 's'), ('q', 't')),
                  ('gamma', ('p',), ('q',)), ('gamma', ('u',), ('t',)),
                  ('eta', ('q',), ('q',)), ('eta', ('s',), ('r',))))]
    c = collections.Counter()
    for sg, f in base:
        c[canon(f, dummies)] += sg
        cm = conjugate(f, True, True, True)
        nvs = sum(1 for n, _, _ in cm if n == 'Vs')
        c[canon(cm, dummies)] += -sg * (1 ** nvs)
    c = collections.Counter({k: v for k, v in c.items() if v})
    ok = audit(c, dummies, "SELF-TEST (constructed p-h covariant list)")
    print("   self-test:", "PASS" if ok else "FAIL -- audit cannot detect covariance")
    return ok


if __name__ == '__main__':
    selftest(('r', 's', 't', 'u'))
    T = extract_file('/root/workspace/refs/TightBindingTCL2/TCL_EE2.h')
    dummies = ('r', 's', 't', 'u')
    groups = collections.defaultdict(collections.Counter)
    allc = collections.Counter()
    for ln, g, s, f in T:
        key = canon(f, dummies)
        groups[g][key] += s
        allc[key] += s
    allc = collections.Counter({k: v for k, v in allc.items() if v})
    audit(allc, dummies, "ALL 72 stamped terms, summed")
    sec = collections.Counter()
    for ln, g, sg, f in T:
        if g in (61, 62, 63, 64):
            sec[canon(f, dummies)] += sg
    sec = collections.Counter({k: v for k, v in sec.items() if v})
    audit(sec, dummies, "RhoDotSec only (groups 61-64)")
    full = collections.Counter()
    for ln, g, sg, f in T:
        if g in (71, 72, 73, 74):
            full[canon(f, dummies)] += sg
    full = collections.Counter({k: v for k, v in full.items() if v})
    audit(full, dummies, "RhoDotSecFull only (groups 71-74)")
    for g in sorted(groups):
        c = collections.Counter({k: v for k, v in groups[g].items() if v})
        if c:
            audit(c, dummies, f"group {g}")
