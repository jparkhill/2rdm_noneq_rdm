#!/usr/bin/env python3
"""The MKNO flavor of lambda_3 = 0 versus ours: which variables are propagated.

The previous commit showed |A - B| = 2e-15 for two ways of imposing lambda_3 = 0. That was
weaker than advertised: once the GNO term list came out identical to the bare one, version B
was a rewriting of version A rather than an independent theory, so the agreement checked the
substitution bookkeeping and nothing more.

The flavors that actually differ are distinguished by what is propagated.

    ours  : propagate gamma_2 alone, with gamma_3 replaced by its separable part.
            gamma_1 is then whatever the contraction of gamma_2 says it is.
    MKNO  : propagate gamma_1 and lambda_2. gamma_1-dot needs NO closure, since the rank-1
            EOM depends only on h, V, gamma_1 and gamma_2. Only lambda_2-dot is truncated.

These coincide for all time only if the truncated gamma_2-dot contracts to the exact
gamma_1-dot. The contraction of the eight gamma_3 terms is not zero: it collapses through
sum_x gamma_3[pvx,qwx] = (N-2) gamma_2, and the separable reconstruction violates that sum
rule by 3.4e-2. So the two flavors should disagree, and the size of the disagreement is the
inconsistency between the propagated gamma_1 and the exact one.

An earlier tool, tools/pinning_gamma3.py, contracted gamma_2-dot and reported all closures
agreeing at 1e-16, concluding gamma_1-dot was gamma_3-blind. That conclusion is tested here
directly rather than assumed.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from qdiss.symbolic import mkno
from tools.eom1 import one as ONE1, two as TWO1
from tools.eom2_coherent import ONE as ONE2, TWO as TWO2

M, N = 6, 3
dim = 1 << M
UP3, LO3 = ("p", "v", "x"), ("q", "w", "y")
SEP3 = [(s, ",".join("".join(u) + "".join(l) for n, r, u, l in f) + "->pvxqwy",
         [(n, r) for n, r, u, l in f])
        for s, f in mkno.gamma_to_cumulants(UP3, LO3)
        if not any(n == "L" and r == 3 for n, r, u, l in f)]


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


EYE = np.eye(M, dtype=complex)


def evaluate(terms, Hop, hsub, pref, tabs, out_idx):
    tot = np.zeros((M,) * len(out_idx), dtype=complex)
    for s, facs in terms:
        subs = [hsub]; ops = [Hop]
        for n, rk, u, l in facs:
            subs.append("".join(u) + "".join(l))
            ops.append(EYE if n == "D" else tabs[rk])
        tot = tot + s * pref * np.einsum(",".join(subs) + "->" + out_idx, *ops, optimize=True)
    return tot


print("%4s | %11s | %11s %11s %11s"
      % ("H", "|g1dot|", "exact", "lambda_3=0", "gamma_3=0"))
acc = {"exact": [], "l3": [], "g3": []}
for tr in range(6):
    r = np.random.default_rng(400 + tr)
    A_ = r.normal(size=(M, M)) + 1j * r.normal(size=(M, M))
    h = 0.25 * (A_ + A_.conj().T)
    Bt = r.normal(size=(M,) * 4) + 1j * r.normal(size=(M,) * 4)
    V = Bt - np.transpose(Bt, (1, 0, 2, 3))
    V = V - np.transpose(V, (0, 1, 3, 2))
    V = 0.25 * (V + np.transpose(V.conj(), (2, 3, 0, 1)))
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
    l2 = g2 - (np.einsum("pq,vw->pvqw", g1, g1) - np.einsum("pw,vq->pvqw", g1, g1))
    t = {("G", 1): g1, ("L", 2): l2}
    g3s = np.zeros((M,) * 6, dtype=complex)
    for s, sub, keys in SEP3:
        g3s = g3s + s * np.einsum(sub, *[t[k] for k in keys], optimize=True)

    # MKNO side: gamma_1-dot is EXACT, needing only gamma_1 and gamma_2
    tb = {1: g1, 2: g2}
    d1_exact = 1j * (evaluate(ONE1, h, "rs", 1.0, tb, "pq")
                     + evaluate(TWO1, V, "rtsu", 0.25, tb, "pq"))
    row = []
    for key, g3u in (("exact", g3), ("l3", g3s), ("g3", np.zeros((M,) * 6, dtype=complex))):
        tb2 = {1: g1, 2: g2, 3: g3u}
        d2 = 1j * (evaluate(ONE2, h, "rs", 1.0, tb2, "pvqw")
                   + evaluate(TWO2, V, "rtsu", 0.25, tb2, "pvqw"))
        d1_from_d2 = np.einsum("pxqx->pq", d2) / (N - 1)
        err = np.abs(d1_from_d2 - d1_exact).max()
        acc[key].append(err); row.append(err)
    print("%4d | %11.4f | %11.3e %11.3e %11.3e"
          % (tr, np.abs(d1_exact).max(), row[0], row[1], row[2]))

print()
print("discrepancy between the contraction of gamma_2-dot and the exact gamma_1-dot:")
for k, nm in (("exact", "exact gamma_3"), ("l3", "lambda_3 = 0"), ("g3", "gamma_3 = 0")):
    a = np.array(acc[k])
    print("   %-14s worst %.4e   mean %.4e" % (nm, a.max(), a.mean()))
print()
print("A nonzero entry means propagating gamma_2 alone and propagating (gamma_1, lambda_2)")
print("are DIFFERENT dynamics under the same nominal lambda_3 = 0 truncation, because the")
print("truncated gamma_2-dot no longer contracts to the gamma_1-dot the theory requires.")
