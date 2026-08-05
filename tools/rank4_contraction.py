#!/usr/bin/env python3
"""The general-k cumulant contraction rule: derive rank 4 to see the pattern.

Two data points are not a pattern. Known:

    k=2   sum_x lambda_2[px,qx]     = -(gamma_1 - gamma_1^2) = -gamma eta      a PRODUCT
    k=3   sum_y lambda_3[pvy,qwy]   = E lambda_2 + lambda_2 E,  E = gamma - 1/2 = (gamma-eta)/2

both verified at 1e-16. The rank-2 form vanishes at both occupation boundaries, which is what
makes it a blocking factor; the rank-3 form is maximal there. So "blocking" is not
order-general and something else is. Rank 4 decides which.

Method. gamma_4 = sum over set partitions of products of gamma_1 and cumulants, plus lambda_4.
Contracting the fourth index pair,

    sum_z lambda_4[pvxz,qwyz] = (N-3) gamma_3 - sum_z (every orbit that is not lambda_4)

Each orbit is contracted directly by writing its einsum with the fourth upper and lower index
identified, so no rank-4 tensor is ever built. Candidate closed forms are then tested against
that, rather than guessed.

Ground truth: exact RDMs from qdiss.eom.fullci.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
from qdiss.eom import fullci
from qdiss.symbolic import mkno

UP2, LO2 = ("p", "v"), ("q", "w")
UP3, LO3 = ("p", "v", "x"), ("q", "w", "y")
UP4, LO4 = ("p", "v", "x", "z"), ("q", "w", "y", "u")


def wedge(a, b):
    return np.einsum("pq,vw->pvqw", a, b) - np.einsum("pw,vq->pvqw", a, b)


def orbit_terms(up, lo, drop_top):
    """(sign, subscript, factor keys) for the cumulant expansion, optionally dropping
    the top cumulant block."""
    k = len(up)
    out = []
    for s, facs in mkno.gamma_to_cumulants(up, lo):
        if drop_top and any(n == "L" and r == k for n, r, u, l in facs):
            continue
        out.append((s, [("".join(u), "".join(l)) for n, r, u, l in facs],
                    [(n, r) for n, r, u, l in facs]))
    return out


def build(terms, tens, out_sub, M, ident=None):
    """Evaluate a term list; ident=(a,b) identifies index a with b and sums it."""
    tot = np.zeros((M,) * len(out_sub), dtype=complex)
    for s, pairs, keys in terms:
        subs = []
        for u, l in pairs:
            su = u + l
            if ident:
                su = su.replace(ident[1], ident[0])
            subs.append(su)
        tot = tot + s * np.einsum(",".join(subs) + "->" + out_sub,
                                  *[tens[k] for k in keys], optimize=True)
    return tot


for M, N in ((6, 4), (6, 5), (8, 4)):
    rng = np.random.default_rng(9)
    psi = fullci.random_state(M, N, rng)
    g1 = fullci.rdm(psi, M, 1)
    g2 = fullci.rdm(psi, M, 2)
    g3 = fullci.rdm(psi, M, 3)
    l2 = g2 - wedge(g1, g1)
    T = {("G", 1): g1, ("L", 2): l2}
    l3 = g3 - build(orbit_terms(UP3, LO3, True), T, "pvxqwy", M)
    T[("L", 3)] = l3
    t4 = orbit_terms(UP4, LO4, True)
    print("=== M=%d N=%d ===  rank-4 expansion has %d non-lambda_4 terms" % (M, N, len(t4)))
    sub4 = build(t4, T, "pvxqwy", M, ident=("z", "u"))
    tr_l4 = (N - 3) * g3 - sub4

    I = np.eye(M, dtype=complex)
    E = g1 - 0.5 * I

    def anti(A, X, nup):
        """(A (+) A) acting on all slots of a rank-nup cumulant X."""
        if nup == 2:
            return (np.einsum("py,yvqw->pvqw", A, X) + np.einsum("vy,pyqw->pvqw", A, X)
                    + np.einsum("pvyw,yq->pvqw", X, A) + np.einsum("pvqy,yw->pvqw", X, A))
        return (np.einsum("pt,tvxqwy->pvxqwy", A, X)
                + np.einsum("vt,ptxqwy->pvxqwy", A, X)
                + np.einsum("xt,pvtqwy->pvxqwy", A, X)
                + np.einsum("pvxtwy,tq->pvxqwy", X, A)
                + np.einsum("pvxqty,tw->pvxqwy", X, A)
                + np.einsum("pvxqwt,ty->pvxqwy", X, A))

    cand = anti(E, l3, 3)
    resid = tr_l4 - cand
    print("   |sum_z lambda_4|        = %.5f" % np.abs(tr_l4).max())
    print("   |E-anticommutator|      = %.5f" % np.abs(cand).max())
    print("   |residual|              = %.5f" % np.abs(resid).max())
    best = None
    for c in np.linspace(0.0, 2.0, 21):
        r = np.abs(tr_l4 - c * cand).max()
        if best is None or r < best[0]:
            best = (r, c)
    print("   best scalar multiple    : c = %.2f leaves %.5f" % (best[1], best[0]))
    print("   rank-3 check: |sum_y lambda_3 - (E lam2 + lam2 E)| = %.2e"
          % np.abs(np.einsum("pvyqwy->pvqw", l3) - anti(E, l2, 2)).max())
    print()
