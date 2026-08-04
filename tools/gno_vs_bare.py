#!/usr/bin/env python3
"""Cross-validate the bare-Wick and GNO derivations of the 2-RDM EOM, then compare the
two ways of imposing lambda_3 = 0.

Part 1, equivalence with no truncation. The bare route (tools/eom2_coherent.py) anticommutes
to the vacuum and yields deltas times bare moments. The GNO route instead expands each bare
operator into generalized-normal-ordered operators,

    a^P_Q = sum_m sum_{|A|=|B|=m} eps(A) eps(B) gamma_m^{P_A}_{Q_B} {a^{P\\A}_{Q\\B}}

multiplies the GNO operators with the generalized Wick product rule (mkno.gno_product), and
takes the expectation, where every surviving term is one whose GNO factor has been fully
contracted away since <{a^..}> = 0 by construction. That path exercises the eps signs, the
recursive gno_to_bare, and the product rule -- code the bare route never touches -- so
numerical agreement is a genuine cross-validation rather than a tautology.

Part 2, the two truncations. Version A is what every closure tool in this repo does: take
the BARE EOM and substitute gamma_3 -> gamma_to_cumulants with the lambda_3 block deleted,
leaving gamma_2 alone. Version B is the MKNO-native one: take the GNO EOM, push it fully to
cumulant form with to_cumulant_form, and delete every term containing lambda_3. These are
the same function of (gamma_1, lambda_2) if both derivations are exact, so a disagreement
would mean one of the two truncations is not what it claims to be.
"""
import itertools
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from qdiss.symbolic import mkno
from tools.eom2_coherent import ONE as BARE_ONE, TWO as BARE_TWO

M, N = 6, 3
dim = 1 << M
D2 = ("p", "v"), ("q", "w")


def gno_expect_product(a, b):
    """<a^{a} a^{b}> with a, b bare operators given as (upper, lower)."""
    out = []
    for s1, f1 in mkno.bare_to_gno(*a):
        ga = tuple(f for f in f1 if f[0] != mkno.GNO)
        oa = [f for f in f1 if f[0] == mkno.GNO]
        for s2, f2 in mkno.bare_to_gno(*b):
            gb = tuple(f for f in f2 if f[0] != mkno.GNO)
            ob = [f for f in f2 if f[0] == mkno.GNO]
            if oa and ob:
                for s3, f3 in mkno.gno_product((oa[0][2], oa[0][3]), (ob[0][2], ob[0][3])):
                    if any(f[0] == mkno.GNO for f in f3):
                        continue
                    out.append((s1 * s2 * s3, ga + gb + tuple(f3)))
            elif oa or ob:
                continue
            else:
                out.append((s1 * s2, ga + gb))
    return mkno.collect(out)


def gno_commutator(x, y):
    fwd = gno_expect_product(x, y)
    rev = gno_expect_product(y, x)
    return mkno.collect(list(fwd) + [(-s, f) for s, f in rev])


GNO_ONE = gno_commutator((("r",), ("s",)), D2)
GNO_TWO = gno_commutator((("r", "t"), ("s", "u")), D2)
print("term counts")
print("   bare route : one-body %d, two-body %d" % (len(BARE_ONE), len(BARE_TWO)))
print("   GNO  route : one-body %d, two-body %d" % (len(GNO_ONE), len(GNO_TWO)))


def ranks(terms):
    d = {}
    for s, f in terms:
        for n, r, u, l in f:
            d[(n, r)] = d.get((n, r), 0) + 1
    return ", ".join("%s%d x%d" % (k[0], k[1], v) for k, v in sorted(d.items()))


print("   bare two-body factors : %s" % ranks(BARE_TWO))
print("   GNO  two-body factors : %s" % ranks(GNO_TWO))
print()


def ad(p, v):
    out = np.zeros_like(v); bit = 1 << p
    for s in np.nonzero(v)[0]:
        if s & bit: continue
        out[s | bit] += ((-1.0) ** bin(s & (bit - 1)).count("1")) * v[s]
    return out


def aa(p, v):
    out = np.zeros_like(v); bit = 1 << p
    for s in np.nonzero(v)[0]:
        if not (s & bit): continue
        out[s ^ bit] += ((-1.0) ** bin(s & (bit - 1)).count("1")) * v[s]
    return out


r = np.random.default_rng(11)
popc = np.array([bin(s).count("1") for s in range(dim)])
psi = np.where(popc == N, r.normal(size=dim) + 1j * r.normal(size=dim), 0.0)
psi = psi / np.linalg.norm(psi)
Q1 = np.array([aa(p, psi) for p in range(M)])
g1 = Q1.conj() @ Q1.T
Q2 = np.array([aa(v, aa(p, psi)) for p in range(M) for v in range(M)])
g2 = (Q2.conj() @ Q2.T).reshape(M, M, M, M)
Q3 = np.array([aa(x, aa(v, aa(p, psi)))
               for p in range(M) for v in range(M) for x in range(M)])
g3 = (Q3.conj() @ Q3.T).reshape((M,) * 6)

UP3, LO3 = ("p", "v", "x"), ("q", "w", "y")
FULL3 = [(s, ",".join("".join(u) + "".join(l) for n, rk, u, l in f) + "->pvxqwy",
          [(n, rk) for n, rk, u, l in f]) for s, f in mkno.gamma_to_cumulants(UP3, LO3)]
SEP3 = [t for t in FULL3 if ("L", 3) not in t[2]]
l2 = g2 - (np.einsum("pq,vw->pvqw", g1, g1) - np.einsum("pw,vq->pvqw", g1, g1))
tens = {("G", 1): g1, ("L", 2): l2}
g3sep = np.zeros((M,) * 6, dtype=complex)
for s, sub, keys in SEP3:
    g3sep = g3sep + s * np.einsum(sub, *[tens[k] for k in keys], optimize=True)
l3 = g3 - g3sep

A_ = r.normal(size=(M, M)) + 1j * r.normal(size=(M, M))
h = 0.25 * (A_ + A_.conj().T)
Bt = r.normal(size=(M,) * 4) + 1j * r.normal(size=(M,) * 4)
V = Bt - np.transpose(Bt, (1, 0, 2, 3))
V = V - np.transpose(V, (0, 1, 3, 2))
V = 0.25 * (V + np.transpose(V.conj(), (2, 3, 0, 1)))
EYE = np.eye(M, dtype=complex)


def evaluate(terms, Hop, hsub, pref, tabs):
    tot = np.zeros((M,) * 4, dtype=complex)
    for s, facs in terms:
        subs = [hsub]; ops = [Hop]
        ok = True
        for n, rk, u, l in facs:
            key = ("D", 1) if n == "D" else (n, rk)
            if key not in tabs:
                ok = False; break
            subs.append("".join(u) + "".join(l))
            ops.append(tabs[key])
        if not ok:
            continue
        tot = tot + s * pref * np.einsum(",".join(subs) + "->pvqw", *ops, optimize=True)
    return tot


TAB = {("D", 1): EYE, ("G", 1): g1, ("G", 2): g2, ("G", 3): g3,
       ("L", 2): l2, ("L", 3): l3}


def assemble(one, two, tabs):
    return 1j * (evaluate(one, h, "rs", 1.0, tabs) + evaluate(two, V, "rtsu", 0.25, tabs))


bare = assemble(BARE_ONE, BARE_TWO, TAB)
gno = assemble(GNO_ONE, GNO_TWO, TAB)
print("=== PART 1: untruncated equivalence ===")
print("   |bare EOM|                 = %.6f" % np.abs(bare).max())
print("   |GNO EOM|                  = %.6f" % np.abs(gno).max())
print("   |bare - GNO|               = %.3e" % np.abs(bare - gno).max())
print()

cum_one = mkno.to_cumulant_form(GNO_ONE)
cum_two = mkno.to_cumulant_form(GNO_TWO)
print("   GNO EOM pushed to cumulant form: one-body %d terms, two-body %d terms"
      % (len(cum_one), len(cum_two)))
print("   cumulant-form factors : %s" % ranks(cum_two))
cum = assemble(cum_one, cum_two, TAB)
print("   |cumulant-form - GNO|      = %.3e" % np.abs(cum - gno).max())
print()

print("=== PART 2: two ways of imposing lambda_3 = 0 ===")
TAB_A = dict(TAB); TAB_A[("G", 3)] = g3sep
verA = assemble(BARE_ONE, BARE_TWO, TAB_A)
noL3_one = [t for t in cum_one if not any(n == "L" and rk == 3 for n, rk, u, l in t[1])]
noL3_two = [t for t in cum_two if not any(n == "L" and rk == 3 for n, rk, u, l in t[1])]
print("   version B drops %d of %d one-body and %d of %d two-body cumulant terms"
      % (len(cum_one) - len(noL3_one), len(cum_one),
         len(cum_two) - len(noL3_two), len(cum_two)))
verB = assemble(noL3_one, noL3_two, TAB)
print("   |version A (bare + gamma_3 -> separable)| = %.6f" % np.abs(verA).max())
print("   |version B (GNO cumulant, drop lambda_3)| = %.6f" % np.abs(verB).max())
print("   |A - B|                                   = %.3e" % np.abs(verA - verB).max())
print("   |A - exact|                               = %.6f" % np.abs(verA - bare).max())
print("   |B - exact|                               = %.6f" % np.abs(verB - bare).max())
