#!/usr/bin/env python3
"""The coherent 1-body EOM by the Heisenberg route, and its TDHF limit.

Recipe (docs/CONTINUITY.md 3b): Heisenberg EOM for the operator whose
expectation IS the propagated quantity, then take the expectation.

    i d<a^p_q>/dt = <[H, a^p_q]>,   H = sum h^r_s a^r_s + (1/2) sum V a^{rs}_{tu}

The expectation of a normal-ordered BARE string is just gamma of that rank, so
no GNO machinery is needed here: anticommute, then read off gammas.

Verified numerically against explicit Fock-space matrices before any claim about
TDHF is made.
"""
import itertools
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from qdiss.symbolic import mkno
from gno import Fock


def bare_product(u1, l1, u2, l2):
    """a^u1_l1 . a^u2_l2 -> [(sign, deltas, up, lo)] with everything normal ordered."""
    word = mkno.ops_of(u1, l1) + mkno.ops_of(u2, l2)
    out = []
    for s, deltas, cre, ann in mkno.normal_order(word):
        if len(cre) != len(ann):
            continue
        out.append((s, deltas, tuple(cre), tuple(reversed(ann))))
    return out


def commutator_expectation(u1, l1, u2, l2):
    """<[a^u1_l1, a^u2_l2]> as a term list over gamma and delta factors."""
    terms = []
    for s, d, up, lo in bare_product(u1, l1, u2, l2):
        facs = d + ((("G", len(up), up, lo),) if up else ())
        terms.append((s, facs))
    for s, d, up, lo in bare_product(u2, l2, u1, l1):
        facs = d + ((("G", len(up), up, lo),) if up else ())
        terms.append((-s, facs))
    return mkno.collect(terms)


def show(label, terms):
    print("%s : %d terms" % (label, len(terms)))
    for s, f in terms:
        parts = []
        for n, r, u, l in f:
            sym = {"G": "g", "D": "d"}.get(n, n)
            parts.append("%s%d[%s;%s]" % (sym, r, "".join(u), "".join(l)))
        print("   %+d  %s" % (s, " ".join(parts)))


# --- one-body piece of H: gives the [h, gamma] commutator ---
one = commutator_expectation(("r",), ("s",), ("p",), ("q",))
show("<[a^r_s, a^p_q]>", one)
print()

# --- two-body piece of H: raises gamma_2 ---
two = commutator_expectation(("r", "t"), ("s", "u"), ("p",), ("q",))
show("<[a^rt_su, a^p_q]>", two)
print()

# --- numerical check of both, against explicit Fock-space matrices ---
M, N = 6, 3
rng = np.random.default_rng(11)
dim = 1 << M
popc = np.array([bin(s).count("1") for s in range(dim)])
psi = np.where(popc == N, rng.normal(size=dim), 0.0)
psi = psi / np.linalg.norm(psi)
F = Fock(M, psi)


def ev_bare(up, lo):
    return F.rdm(list(up), list(lo))


def evaluate(terms, sub):
    tot = 0.0
    for c, facs in terms:
        val = float(c)
        for n, r, u, l in facs:
            iu = [sub[i] for i in u]
            il = [sub[i] for i in l]
            if n == "D":
                val *= 1.0 if iu[0] == il[0] else 0.0
            elif n == "G":
                val *= ev_bare(iu, il)
            if val == 0.0:
                break
        tot += val
    return tot


def brute(u1, l1, u2, l2, sub):
    a = mkno.ops_of([sub[i] for i in u1], [sub[i] for i in l1])
    b = mkno.ops_of([sub[i] for i in u2], [sub[i] for i in l2])

    def apply(word, v):
        for kind, idx in reversed(word):
            v = F.ad(idx, v) if kind == "c" else F.aa(idx, v)
        return v

    ab = apply(a, apply(b, psi))
    ba = apply(b, apply(a, psi))
    return float(psi @ (ab - ba))


worst1 = worst2 = 0.0
for _ in range(300):
    sub = {k: int(rng.integers(M)) for k in "pqrstu"}
    worst1 = max(worst1, abs(evaluate(one, sub)
                            - brute(("r",), ("s",), ("p",), ("q",), sub)))
    worst2 = max(worst2, abs(evaluate(two, sub)
                            - brute(("r", "t"), ("s", "u"), ("p",), ("q",), sub)))
print("numerical check over 300 random index tuples, M=%d N=%d:" % (M, N))
print("  one-body commutator : max error %.3e" % worst1)
print("  two-body commutator : max error %.3e" % worst2)
print("  VERDICT:", "BOTH EXACT" if max(worst1, worst2) < 1e-12 else "MISMATCH")
