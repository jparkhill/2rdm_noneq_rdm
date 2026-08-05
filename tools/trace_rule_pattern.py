#!/usr/bin/env python3
"""The single-trace identity at ranks 2, 3, 4: find the rank-4 form and look for a pattern.

Known and verified:
    k=2   sum_x lambda_2[px,qx]   = -(g - g^2)                      = -g eta
    k=3   sum_y lambda_3[pvy,qwy] = E lambda_2 + lambda_2 E,  E = g - 1/2

Unknown: the rank-4 form. Rather than guess a single candidate, fit sum_z lambda_4 against the
structures that can produce a 6-index antisymmetric object from g, lambda_2 and lambda_3:

    lam3                     lambda_3 itself
    anti(g, lam3)            g contracted into each of the six slots of lambda_3
    chainA, chainB           two contraction topologies of lambda_2 with lambda_2,
                             fully antisymmetrised
    orbB(g), orbB(g^2)       the rank-3 gamma-wedge-lambda_2 orbit with g or g^2 in the
                             2-index slot
    orbA(g), orbA(g^2)       the rank-3 triple-wedge orbit

Coefficients are required to be COMMON across several (M,N), by stacking the equations before
solving. A fit that only works at one system size is an accident; one that holds across all of
them with the same numbers is the identity.

Ground truth: exact RDMs from qdiss.eom.fullci. Nothing here is imposed as a condition -- the
double-trace condition explored earlier is implied by the single trace and is redundant.
"""
import itertools
import sys
import numpy as np
sys.path.insert(0, ".")
from qdiss.eom import fullci
from qdiss.symbolic import mkno

UP3, LO3 = ("p", "v", "x"), ("q", "w", "y")
UP4, LO4 = ("p", "v", "x", "z"), ("q", "w", "y", "u")


def terms(up, lo, drop_top):
    k = len(up)
    out = []
    for s, facs in mkno.gamma_to_cumulants(up, lo):
        if drop_top and any(n == "L" and r == k for n, r, u, l in facs):
            continue
        out.append((s, [("".join(u), "".join(l)) for n, r, u, l in facs],
                    [(n, r) for n, r, u, l in facs]))
    return out


T3 = terms(UP3, LO3, True)
T4 = terms(UP4, LO4, True)


def build(tl, tens, out_sub, M, ident=None):
    tot = np.zeros((M,) * len(out_sub), dtype=complex)
    for s, pairs, keys in tl:
        subs = []
        for u, l in pairs:
            su = u + l
            if ident:
                su = su.replace(ident[1], ident[0])
            subs.append(su)
        tot = tot + s * np.einsum(",".join(subs) + "->" + out_sub,
                                  *[tens[k] for k in keys], optimize=True)
    return tot


def wedge(a, b):
    return np.einsum("pq,vw->pvqw", a, b) - np.einsum("pw,vq->pvqw", a, b)


def anti3(A, X):
    return (np.einsum("pt,tvxqwy->pvxqwy", A, X)
            + np.einsum("vt,ptxqwy->pvxqwy", A, X)
            + np.einsum("xt,pvtqwy->pvxqwy", A, X)
            + np.einsum("pvxtwy,tq->pvxqwy", X, A)
            + np.einsum("pvxqty,tw->pvxqwy", X, A)
            + np.einsum("pvxqwt,ty->pvxqwy", X, A))


def asym6(T):
    out = np.zeros_like(T)
    for pu in itertools.permutations(range(3)):
        su = np.linalg.det(np.eye(3)[list(pu)])
        for pl in itertools.permutations(range(3)):
            sl = np.linalg.det(np.eye(3)[list(pl)])
            out += su * sl * np.transpose(T, tuple(pu) + tuple(3 + i for i in pl))
    return out / 36.0


CASES = ((6, 3), (6, 4), (8, 3), (8, 4))
rows, rhs, labels = [], [], None
for M, N in CASES:
    rng = np.random.default_rng(9)
    psi = fullci.random_state(M, N, rng)
    g1 = fullci.rdm(psi, M, 1)
    g2 = fullci.rdm(psi, M, 2)
    g3 = fullci.rdm(psi, M, 3)
    l2 = g2 - wedge(g1, g1)
    tens = {("G", 1): g1, ("L", 2): l2}
    l3 = g3 - build(T3, tens, "pvxqwy", M)
    tens[("L", 3)] = l3
    target = (N - 3) * g3 - build(T4, tens, "pvxqwy", M, ident=("z", "u"))
    g1sq = g1 @ g1
    chainA = asym6(np.einsum("pvqt,txwy->pvxqwy", l2, l2))
    chainB = asym6(np.einsum("pvtw,xtqy->pvxqwy", l2, l2))
    cands = [("lam3", l3),
             ("anti(g,lam3)", anti3(g1, l3)),
             ("chainA", chainA),
             ("chainB", chainB),
             ("orbB(g)", build([t for t in T3 if ("L", 2) in t[2]],
                               {("G", 1): g1, ("L", 2): l2}, "pvxqwy", M)),
             ("orbB(g2)", build([t for t in T3 if ("L", 2) in t[2]],
                                {("G", 1): g1sq, ("L", 2): l2}, "pvxqwy", M)),
             ("orbA(g)", build([t for t in T3 if ("L", 2) not in t[2]],
                               {("G", 1): g1}, "pvxqwy", M)),
             ("orbA(g2)", build([t for t in T3 if ("L", 2) not in t[2]],
                                {("G", 1): g1sq}, "pvxqwy", M))]
    if labels is None:
        labels = [c[0] for c in cands]
    rows.append(np.stack([np.concatenate([c[1].real.ravel(), c[1].imag.ravel()])
                          for c in cands], axis=1))
    rhs.append(np.concatenate([target.real.ravel(), target.imag.ravel()]))
    print("(%d,%d)  |sum_z lambda_4| = %.5f" % (M, N, np.abs(target).max()))

print()
A = np.concatenate(rows, axis=0)
b = np.concatenate(rhs)
sol, res, rank, sv = np.linalg.lstsq(A, b, rcond=None)
pred = A @ sol
print("common-coefficient fit across %d systems, %d equations, basis rank %d"
      % (len(CASES), A.shape[0], rank))
for lb, c in zip(labels, sol):
    print("   %-14s %+12.6f" % (lb, c))
print()
print("   relative residual  %.4e   (|b| max %.5f)"
      % (np.abs(pred - b).max() / np.abs(b).max(), np.abs(b).max()))
print("   singular values    %s" % np.array2string(sv, precision=3))
print()
print("A residual at machine precision with these coefficients would BE the rank-4 identity.")
print("A large residual means the structure list is incomplete, not that no rule exists.")
