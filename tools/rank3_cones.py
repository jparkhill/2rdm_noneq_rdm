#!/usr/bin/env python3
r"""Rank 3: the five Gram cones as affine functions of gamma_1, gamma_2, gamma_3, verified against direct Grams.

gamma_k[P; Q] = <a+_{P1}..a+_{Pk} a_{Qk}..a_{Q1}> in this repo's convention (ket indices first, then bra).
Row/column labels are index triples. By anticommuting the strings to normal order (mkno.normal_order), each Gram is
a signed sum of deltas times gamma_m, m <= 3. This file derives the forms symbolically with the same Wick engine used
for the EOMs, evaluates them, and checks against the wavefunction Grams at 1e-15. Output: the term lists, for use in
the rank-3 barrier (cone-rate maps and adjoints follow by linearity).
"""
import sys, itertools
import numpy as np
sys.path.insert(0, "."); sys.path.insert(0, "tools")
import commutator_structure as cs
from qdiss.symbolic import mkno
from qdiss.eom import fullci

CONES = {   # name: (row string, column string) as ("c"/"a", index) lists; Gram[row, col] = <row^dag col>
    "D3  <a+a+a+ aaa>":   ([("a", "x"), ("a", "v"), ("a", "p")], [("a", "y"), ("a", "w"), ("a", "q")]),
    "Q3  <aaa a+a+a+>":   ([("c", "x"), ("c", "v"), ("c", "p")], [("c", "y"), ("c", "w"), ("c", "q")]),
    "T1  <(a+aa)^dag a+aa>":   ([("a", "x"), ("a", "v"), ("c", "p")], [("a", "y"), ("a", "w"), ("c", "q")]),
    "T2  <(a+a+a)^dag a+a+a>": ([("a", "x"), ("c", "v"), ("c", "p")], [("a", "y"), ("c", "w"), ("c", "q")]),
}


def dagger(word):
    return [("a" if k == "c" else "c", i) for k, i in reversed(word)]


def gram_terms(row, col):
    """<row^dag col> as (coef, deltas, gamma factor) with free indices p,v,x (row) and q,w,y (col)."""
    word = dagger(row) + col
    acc = {}
    for s, deltas, cre, ann in mkno.normal_order(word):
        if len(cre) != len(ann): continue
        g = ("G", len(cre), tuple(cre), tuple(reversed(ann))) if cre else None
        sg, key = mkno.canon(1, deltas + ((g,) if g else ()))
        acc[key] = acc.get(key, 0) + s * sg
    return [(c, tuple(f for f in k if f[0] == "D"), ([f for f in k if f[0] == "G"] or [None])[0]) for k, c in acc.items() if c != 0]


def evaluate(terms, tab, M):
    I = np.eye(M, dtype=complex); out = "pvxqwy"; tot = np.zeros((M,) * 6, dtype=complex)
    for c, deltas, g in terms:
        rep = {}
        for _, _, (i,), (j,) in deltas:
            ri, rj = rep.get(i, i), rep.get(j, j)
            if ri != rj:
                for k_, v_ in list(rep.items()):
                    if v_ == rj: rep[k_] = ri
                rep[rj] = ri
        r = lambda i: rep.get(i, i)
        subs, ops = [], []
        if g is not None:
            _, m, cre, ann = g; subs.append("".join(r(i) for i in cre) + "".join(r(i) for i in ann)); ops.append(tab[("G", m)])
        for o in out:
            if r(o) != o: subs.append(o + r(o)); ops.append(I)
        # deltas between two output indices where neither is substituted away cannot occur after rep; deltas among
        # summed-over indices do not exist here (no internal sums). Output indices absent everywhere: only if the
        # term is a pure product of deltas among outputs, covered above.
        present = set("".join(subs))
        for o in out:
            if o not in present: raise RuntimeError((c, deltas, g))
        tot += c * np.einsum(",".join(subs) + "->" + out, *ops, optimize=True)
    return tot


def main():
    atom, N, spin = cs.H3
    label, dm, dets, M, h, V = next(cs.states(atom, N, spin, "H3"))
    w, U = np.linalg.eigh(0.5 * (dm + dm.conj().T)); psi = np.zeros(1 << M, dtype=complex); psi[dets] = U[:, np.argmax(w)]
    tab = {("G", k): fullci.rdm(psi, M, k) for k in (1, 2, 3)}
    for name, (row, col) in CONES.items():
        terms = gram_terms(row, col)
        G = evaluate(terms, tab, M).reshape(M ** 3, M ** 3)
        ops = {"a": fullci.aa, "c": fullci.ad}
        def apply(word, idx):
            v = psi
            for (k, i) in word: v = ops[k](idx[i], v)
            return v
        # Gram[(p,v,x),(q,w,y)] = < row(p,v,x) psi | col(q,w,y) psi >, each string applied right-to-left as written
        # (the LAST element acts first), i.e. word [o1,o2,o3] means o1 o2 o3 |psi>.
        def apply_word(word, idx):
            v = psi
            for (k, i) in reversed(word): v = ops[k](idx[i], v)
            return v
        rows = np.array([apply_word(row, dict(zip("pvx", t))) for t in itertools.product(range(M), repeat=3)])
        cols = np.array([apply_word(col, dict(zip("qwy", t))) for t in itertools.product(range(M), repeat=3)])
        Gd = rows.conj() @ cols.T
        from collections import Counter
        print("%-26s terms: %s   |affine form - direct Gram| = %.1e   min eig %+.1e"
              % (name, dict(Counter((g[1] if g else 0) for _, _, g in terms)), np.abs(G - Gd).max(), np.linalg.eigvalsh(0.5 * (Gd + Gd.conj().T)).min()))
        for c, d, g in sorted(terms, key=lambda t: (-(t[2][1] if t[2] else 0), str(t))):
            print("      %+d  %s  %s" % (c, " ".join("d(%s%s)" % (a[0], b[0]) for _, _, a, b in d), ("gamma_%d[%s;%s]" % (g[1], "".join(g[2]), "".join(g[3]))) if g else "1"))


if __name__ == "__main__":
    main()
