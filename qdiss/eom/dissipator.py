"""Rank-k contraction of a one-body GKSL (Davies) dissipator, by plain Wick reordering.

Generator on rho (docs/BATH-DERIVATION.md, the CP, quadratic form; qdiss.rdm1 is its rank-1, lambda_2 = 0 case):

    D[rho] = sum_j  ( L_j rho L_j^dag  -  1/2 { L_j^dag L_j , rho } ),      L_j = sum_pq l^j_pq a+_p a_q   (one-body)

For an observable O = a+_P a_Q (rank k),   d/dt <O>_D = tr(O D[rho]) = sum_j < L_j^dag O L_j - 1/2 (L_j^dag L_j O + O L_j^dag L_j) >.
Each of the three words has two one-body factors and one rank-k factor; anticommuting it to normal order gives
bare moments gamma_k, gamma_{k+1} (never higher, since two one-body operators add at most two indices and the
words reduce). The rates are absorbed by summing j with a weight: with L_j^dag L_j-type products we need the
two-index tensors

    R[p,q,r,s] = sum_j conj(l^j_pq) l^j_rs      so that  sum_j < (a+_q a_p)  O  (a+_r a_s) > -> R contracted

i.e. every term is  R[p,q,r,s] x (deltas) x gamma_m.  The term list is generated once per rank and evaluated with
einsum exactly as qdiss.eom.evaluate does for H. Validated against the Fock-space generator in tools/bath_rank2.py.
"""
import numpy as np
from qdiss.symbolic import mkno


def _words(rank):
    """Three operator words for tr(O D[rho]) with O = a+_P a_Q, in terms of generic one-body a+_a a_b (left factor,
    conjugated) and a+_c a_d (right factor). Indices: O uses (p,v,...)/(q,w,...); L^dag uses (b -> a): (a+_b a_a) is
    the Hermitian conjugate of l_ab a+_a a_b, so L^dag = sum conj(l_ab) a+_b a_a."""
    up = ("p", "v", "x")[:rank]
    lo = ("q", "w", "y")[:rank]
    O = mkno.ops_of(up, lo)
    Ld = mkno.ops_of(("b",), ("a",))          # conj(l_ab) a+_b a_a
    L = mkno.ops_of(("c",), ("d",))           # l_cd a+_c a_d
    LdL = Ld + L
    return [(1.0, Ld + O + L), (-0.5, LdL + O), (-0.5, O + LdL)], up, lo


def terms(rank):
    """List of (coef, deltas, gamma factor) with free indices a,b (conj l) and c,d (l) plus the observable's."""
    words, up, lo = _words(rank)
    acc = {}
    for coef, word in words:
        for s, deltas, cre, ann in mkno.normal_order(word):
            if len(cre) != len(ann):
                continue
            g = ("G", len(cre), tuple(cre), tuple(reversed(ann))) if cre else None
            sg, key = mkno.canon(1, deltas + ((g,) if g else ()))
            acc[key] = acc.get(key, 0.0) + coef * s * sg
    out = []
    for key, c in acc.items():
        if abs(c) < 1e-12:
            continue                      # the rank-(k+2) bare moments cancel identically among the three words
        deltas = tuple(f for f in key if f[0] == "D")
        gs = [f for f in key if f[0] == "G"]
        out.append((c, deltas, gs[0] if gs else None))
    return out, up, lo


def _sub(deltas):
    """Apply delta substitutions: return a dict index -> representative."""
    rep = {}
    for _, _, (i,), (j,) in deltas:
        ri, rj = rep.get(i, i), rep.get(j, j)
        if ri != rj:
            for k, v in list(rep.items()):
                if v == rj:
                    rep[k] = ri
            rep[rj] = ri
    return rep


def gamma_dot(rank, R, tensors, M):
    """d/dt gamma_rank under the dissipator. R[a,b,c,d] = sum_j conj(l^j_ab) l^j_cd. tensors: ("G",m) -> gamma_m.

    Each term is coef * (product of deltas) * R[abcd] * gamma_m[...]. Deltas are applied by index substitution; a
    delta that survives between two OUTPUT indices, or an output index absent from all operands, is supplied by an
    explicit identity matrix in the einsum."""
    tl, up, lo = terms(rank)
    out_sub = "".join(up) + "".join(lo)
    I = np.eye(M, dtype=complex)
    tot = np.zeros((M,) * (2 * rank), dtype=complex)
    for coef, deltas, g in tl:
        rep = _sub(deltas)
        r = lambda i: rep.get(i, i)
        subs, ops = ["".join(r(i) for i in "abcd")], [R]
        if g is not None:
            _, m, cre, ann = g
            subs.append("".join(r(i) for i in cre) + "".join(r(i) for i in ann))
            ops.append(tensors[("G", m)])
        present = set("".join(subs))
        # output indices: for each output letter o with representative r(o): if r(o) != o, add identity I[o, r(o)];
        # if r(o) == o but o is absent from operands (only possible if a delta among outputs was fully substituted),
        # it will have been mapped to its representative, so this case is covered by the first branch.
        for o in out_sub:
            if r(o) != o:
                subs.append(o + r(o)); ops.append(I)
            elif o not in present:
                # a free output index appearing in no operand and not substituted: impossible for a contraction of
                # this type unless a delta pairs it with another output index, handled above. Guard anyway.
                raise RuntimeError("free index %s unmatched in term %r" % (o, (coef, deltas, g)))
        tot = tot + coef * np.einsum(",".join(subs) + "->" + out_sub, *ops, optimize=True)
    return tot
