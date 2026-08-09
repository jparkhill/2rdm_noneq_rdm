"""Derivation A: ordinary Wick reordering with respect to the physical vacuum.

Uses only a_q c_p = delta_pq - c_p a_q. Every factor is produced by symbolic
anticommutation, so no combinatorial rule is applied by hand and no prefactor is
guessed: the terms carry the signs the algebra gives them.

Output is in BARE MOMENTS: deltas times gamma_k. No cumulant appears. The rank-k EOM
couples gamma_k to gamma_k and gamma_{k+1} and never to lower ranks, so for example
gamma_2-dot contains no gamma_1 at all.

Held against qdiss.eom.fullci, never against qdiss.eom.gwt.
"""
from qdiss.symbolic import mkno
from .evaluate import assemble

METHOD = "plain vacuum Wick, bare moments"

OPERATOR = {1: (("p",), ("q",)), 2: (("p", "v"), ("q", "w")),
            3: (("p", "v", "x"), ("q", "w", "y"))}
OUT_SUBSCRIPT = {1: "pq", 2: "pvqw", 3: "pvxqwy"}


def _bare_product(u1, l1, u2, l2):
    word = mkno.ops_of(u1, l1) + mkno.ops_of(u2, l2)
    out = []
    for s, deltas, cre, ann in mkno.normal_order(word):
        if len(cre) == len(ann):
            out.append((s, deltas, tuple(cre), tuple(reversed(ann))))
    return out


def _commutator_expectation(a, b):
    terms = []
    for s, d, up, lo in _bare_product(a[0], a[1], b[0], b[1]):
        terms.append((s, d + ((("G", len(up), up, lo),) if up else ())))
    for s, d, up, lo in _bare_product(b[0], b[1], a[0], a[1]):
        terms.append((-s, d + ((("G", len(up), up, lo),) if up else ())))
    return mkno.collect(terms)


def terms(rank):
    """(one_body_terms, two_body_terms) for the rank-k EOM, in bare moments."""
    op = OPERATOR[rank]
    return (_commutator_expectation((("r",), ("s",)), op),
            _commutator_expectation((("r", "t"), ("s", "u")), op))


def gamma_dot(rank, h, V, tensors, M):
    one, two = terms(rank)
    return assemble(one, two, h, V, tensors, OUT_SUBSCRIPT[rank], M)
