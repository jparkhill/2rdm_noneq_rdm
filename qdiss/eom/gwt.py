"""Derivation B: Mukherjee-Kutzelnigg generalized normal ordering and the generalized
Wick theorem.

Each bare operator is expanded over generalized-normal-ordered operators,

    a^P_Q = sum_m sum_{|A|=|B|=m} eps(A) eps(B) gamma_m^{P_A}_{Q_B} {a^{P\\A}_{Q\\B}}

the GNO operators are multiplied with the product rule in mkno.gno_product, and the
expectation keeps only terms whose GNO factor has been contracted away, since
<{a^..}> = 0 by construction. The eps factors are the graphical/combinatorial rule, so
this route can go wrong in ways the Wick route cannot, and vice versa.

Two representations are exposed. terms() gives the result in bare moments, directly
comparable to the Wick route. terms_cumulant() pushes it to gamma_1 and cumulants, which
is the MKNO-native form and the only one in which "lambda_3 = 0" is a statement about a
single factor.

Held against qdiss.eom.fullci, never against qdiss.eom.wick.
"""
from qdiss.symbolic import mkno
from .evaluate import assemble

METHOD = "MKNO generalized normal ordering / generalized Wick theorem"

OPERATOR = {1: (("p",), ("q",)), 2: (("p", "v"), ("q", "w")),
            3: (("p", "v", "x"), ("q", "w", "y"))}
OUT_SUBSCRIPT = {1: "pq", 2: "pvqw", 3: "pvxqwy"}


def _expect_product(a, b):
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


def _commutator_expectation(a, b):
    fwd = _expect_product(a, b)
    rev = _expect_product(b, a)
    return mkno.collect(list(fwd) + [(-s, f) for s, f in rev])


def terms(rank):
    op = OPERATOR[rank]
    return (_commutator_expectation((("r",), ("s",)), op),
            _commutator_expectation((("r", "t"), ("s", "u")), op))


def terms_cumulant(rank):
    one, two = terms(rank)
    return mkno.to_cumulant_form(one), mkno.to_cumulant_form(two)


def gamma_dot(rank, h, V, tensors, M, cumulant=False):
    one, two = terms_cumulant(rank) if cumulant else terms(rank)
    return assemble(one, two, h, V, tensors, OUT_SUBSCRIPT[rank], M)
