#!/usr/bin/env python3
"""Symbolic Mukherjee-Kutzelnigg generalized normal ordering.

Everything is derived from the one definition that tools/gno.py already
validates numerically, rather than transcribed from the GWT literature:

    a^P_Q = sum_{m=0..k} sum_{|A|=|B|=m} eps(A) eps(B) gamma_m^{P_A}_{Q_B}
                                         {a^{P\\A}_{Q\\B}}

Read forwards it expands a bare string into GNO operators; solved for the m=0
term it expands a GNO operator into bare strings. Composing the two with
ordinary anticommutation in between gives the product rule.

Term format matches tools/term_audit.canon so tools/verify_stage_operator.py
consumes our output unchanged:

    term   = (sign, (factor, ...))
    factor = (name, rank, upper_tuple, lower_tuple)
    name   = 'G' gamma | 'L' lambda | 'D' delta | 'A' the GNO operator
"""
import functools
import itertools
from collections import defaultdict

GNO = "A"


def eps(positions, k):
    """Sign of moving the given positions (increasing) to the front."""
    s = 1
    for rank, i in enumerate(positions):
        s *= (-1) ** (i - rank)
    return s


def perm_sign(seq):
    """Sign of the permutation that sorts seq, plus the sorted tuple."""
    xs = list(seq)
    s = 1
    for i in range(len(xs)):
        for j in range(len(xs) - 1 - i):
            if xs[j] > xs[j + 1]:
                xs[j], xs[j + 1] = xs[j + 1], xs[j]
                s = -s
    return s, tuple(xs)


def canon(sign, factors):
    """Canonical form: antisymmetric factors index-sorted, factors sorted."""
    s, out = sign, []
    for name, rank, up, lo in factors:
        if name == "D" and rank > 1:
            pr = sorted(zip(up, lo))
            up = tuple(a for a, _ in pr)
            lo = tuple(b for _, b in pr)
        elif rank > 1:
            su, up = perm_sign(up)
            sl, lo = perm_sign(lo)
            s *= su * sl
        out.append((name, rank, up, lo))
    out.sort()
    return s, tuple(out)


def collect(terms):
    """Canonicalize and sum like terms, dropping exact cancellations."""
    acc = defaultdict(int)
    for sign, factors in terms:
        s, key = canon(sign, factors)
        acc[key] += s
    return [(s, k) for k, s in sorted(acc.items()) if s != 0]


def gno_to_bare(up, lo):
    """{a^up_lo} as a flat sum of (sign, gamma_factors, bare_up, bare_lo)."""
    up, lo = tuple(up), tuple(lo)
    k = len(up)
    if k == 0:
        return [(1, (), (), ())]
    out = [(1, (), up, lo)]
    for m in range(1, k + 1):
        for A in itertools.combinations(range(k), m):
            for B in itertools.combinations(range(k), m):
                sgn = -eps(A, k) * eps(B, k)
                gfac = ("G", m, tuple(up[i] for i in A), tuple(lo[j] for j in B))
                ru = tuple(up[i] for i in range(k) if i not in A)
                rl = tuple(lo[j] for j in range(k) if j not in B)
                for s2, g2, bu, bl in gno_to_bare(ru, rl):
                    out.append((sgn * s2, (gfac,) + g2, bu, bl))
    return out


def bare_to_gno(up, lo):
    """a^up_lo as a sum of (sign, factors) over GNO operators times gammas."""
    up, lo = tuple(up), tuple(lo)
    k = len(up)
    out = []
    for m in range(0, k + 1):
        for A in itertools.combinations(range(k), m):
            for B in itertools.combinations(range(k), m):
                sgn = eps(A, k) * eps(B, k)
                facs = []
                if m:
                    facs.append(("G", m, tuple(up[i] for i in A),
                                 tuple(lo[j] for j in B)))
                ru = tuple(up[i] for i in range(k) if i not in A)
                rl = tuple(lo[j] for j in range(k) if j not in B)
                if ru:
                    facs.append((GNO, k - m, ru, rl))
                out.append((sgn, tuple(facs)))
    return out


def ops_of(up, lo):
    """a^up_lo as an operator word: creations then reversed annihilations."""
    return [("c", p) for p in up] + [("a", q) for q in reversed(lo)]


def normal_order(word):
    """Anticommute an operator word into (sign, deltas, creations, annihilations).

    Returns a list of such tuples. Uses a_q c_p = delta_pq - c_p a_q, applied to
    the leftmost creation that sits to the right of an annihilation.
    """
    out = []
    stack = [(1, (), list(word))]
    while stack:
        sign, deltas, w = stack.pop()
        swap = None
        for i in range(len(w) - 1):
            if w[i][0] == "a" and w[i + 1][0] == "c":
                swap = i
                break
        if swap is None:
            cre = tuple(p for kind, p in w if kind == "c")
            ann = tuple(q for kind, q in w if kind == "a")
            out.append((sign, deltas, cre, ann))
            continue
        q = w[swap][1]
        p = w[swap + 1][1]
        rest = w[:swap] + w[swap + 2:]
        stack.append((sign, deltas + (("D", 1, (p,), (q,)),), rest))
        swapped = w[:swap] + [w[swap + 1], w[swap]] + w[swap + 2:]
        stack.append((-sign, deltas, swapped))
    return out


def word_to_bare(cre, ann):
    """Creation/annihilation tuples back to a^up_lo, with the reordering sign."""
    up = tuple(cre)
    lo = tuple(reversed(ann))
    return 1, up, lo


def gno_product(a, b):
    """Product of two GNO operators, expanded back into GNO form.

    Each operand is (upper_tuple, lower_tuple). Returns a collected term list.
    """
    terms = []
    for s1, g1, u1, l1 in gno_to_bare(*a):
        for s2, g2, u2, l2 in gno_to_bare(*b):
            word = ops_of(u1, l1) + ops_of(u2, l2)
            for s3, deltas, cre, ann in normal_order(word):
                if len(cre) != len(ann):
                    continue
                s4, bu, bl = word_to_bare(cre, ann)
                for s5, facs in bare_to_gno(bu, bl):
                    terms.append((s1 * s2 * s3 * s4 * s5, g1 + g2 + deltas + facs))
    return collect(terms)


def to_wolfram(terms, path):
    """Write a term list in the .m format tools/term_audit.parse reads."""
    name_out = {"G": "\\[Gamma]", "L": "\\[Lambda]", "H": "\\[Eta]",
                "D": "\\[Delta]", GNO: "ta"}

    def fac(f):
        n, r, up, lo = f
        u = ",".join(str(x) for x in up)
        l = ",".join(str(x) for x in lo)
        return '{{%s,%d},{%s},{%s}}' % (name_out.get(n, n), r, u, l)

    body = ",".join(
        "{" + ",".join(['{"Sign",%d}' % s] + [fac(f) for f in facs]) + "}"
        for s, facs in terms
    )
    with open(path, "w") as fh:
        fh.write("{" + body + "}")
    return path

@functools.lru_cache(maxsize=None)
def set_partitions(k):
    """All set partitions of range(k), each a tuple of increasing-index blocks."""
    if k == 0:
        return [()]
    if k == 1:
        return [((0,),)]
    out = []
    for p in set_partitions(k - 1):
        for i in range(len(p)):
            out.append(p[:i] + (p[i] + (k - 1,),) + p[i + 1:])
        out.append(p + (((k - 1),),))
    return [tuple(sorted(pp, key=lambda b: b[0])) for pp in out]


def gamma_to_cumulants(up, lo):
    """gamma_k^up_lo as a list of (sign, factors) over gamma_1 and lambda_k."""
    up, lo = tuple(up), tuple(lo)
    k = len(up)
    if k == 1:
        return [(1, (("G", 1, up, lo),))]
    out = []
    for part in set_partitions(k):
        sizes = [len(b) for b in part]
        up_sgn = perm_sign(tuple(i for b in part for i in b))[0]
        for perm in itertools.permutations(range(k)):
            chunks, i = [], 0
            for s in sizes:
                chunks.append(perm[i:i + s])
                i += s
            if any(list(c) != sorted(c) for c in chunks):
                continue
            sgn = up_sgn * perm_sign(perm)[0]
            facs = []
            for blk, ch in zip(part, chunks):
                bu = tuple(up[i] for i in blk)
                bl = tuple(lo[j] for j in ch)
                facs.append(("G" if len(blk) == 1 else "L", len(blk), bu, bl))
            out.append((sgn, tuple(facs)))
    return out


def to_cumulant_form(terms):
    """Rewrite every gamma_k factor with k>1 in terms of gamma_1 and lambdas."""
    out = []
    for sign, facs in terms:
        pieces = [(sign, ())]
        for f in facs:
            name, rank, up, lo = f
            if name != "G" or rank == 1:
                pieces = [(s, fs + (f,)) for s, fs in pieces]
                continue
            expanded = []
            for s, fs in pieces:
                for s2, f2 in gamma_to_cumulants(up, lo):
                    expanded.append((s * s2, fs + f2))
            pieces = expanded
        out.extend(pieces)
    return collect(out)
