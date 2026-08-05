#!/usr/bin/env python3
"""The cumulant partial-trace identities, and which EOM terms violate them.

Part 1. Rank 2, restated and re-verified. From gamma_2 = gamma_1 wedge gamma_1 + lambda_2,
the contraction sum rule sum_x gamma_2[px,qx] = (N-1) gamma_1, and
sum_x (gamma_1 wedge gamma_1)[px,qx] = N gamma_1 - gamma_1^2,

    sum_x lambda_2[px,qx] = -(gamma_1 - gamma_1^2) = -gamma_1 eta_1

Part 2. Rank 3, by the same construction. gamma_3 = A + B + lambda_3 with A the gamma_1^3
orbit and B the gamma_1*lambda_2 orbit, and sum_y gamma_3[pvy,qwy] = (N-2) gamma_2, so

    sum_y lambda_3[pvy,qwy] = (N-2) gamma_2 - sum_y A - sum_y B

The orbit contractions are computed numerically from the orbit tensors, which is exact, and a
hand-derived closed form for sum_y A is then tested against that rather than assumed:

    sum_y A = N (g1 wedge g1) - (g1 wedge g1^2) - (g1^2 wedge g1)

Part 3. Which terms of the lambda_2 EOM break the rank-2 identity. Define the residual

    R = sum_x lambda_2[px,qx] + (gamma_1 - gamma_1^2)

which is zero on any physical state. Its rate is

    Rdot = sum_x lambda_2-dot[px,qx] + gamma_1-dot - (gamma_1-dot gamma_1 + gamma_1 gamma_1-dot)

with lambda_2-dot = gamma_2-dot - (gamma_1-dot wedge gamma_1 + gamma_1 wedge gamma_1-dot).
Exact dynamics keeps Rdot = 0. Splitting gamma_2-dot into its lambda_3-carrying terms and the
rest shows whether the violation comes exclusively from lambda_3, as expected.

Derivation: qdiss.eom.gwt cumulant form, validated against qdiss.eom.fullci at 1e-15.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
from qdiss.eom import fullci, gwt
from qdiss.eom.evaluate import assemble
from qdiss.eom.hamiltonian import random_V, random_h
from qdiss.symbolic import mkno

UP3, LO3 = ("p", "v", "x"), ("q", "w", "y")
ORB = {"A": [], "B": []}
for s, f in mkno.gamma_to_cumulants(UP3, LO3):
    if any(n == "L" and r == 3 for n, r, u, l in f):
        continue
    k = "B" if any(n == "L" and r == 2 for n, r, u, l in f) else "A"
    ORB[k].append((s, ",".join("".join(u) + "".join(l) for n, r, u, l in f) + "->pvxqwy",
                   [(n, r) for n, r, u, l in f]))
CUM1, CUM2 = gwt.terms_cumulant(1), gwt.terms_cumulant(2)
L3_TERMS = [t for t in CUM2[1] if any(n == "L" and r == 3 for n, r, u, l in t[1])]
NO_L3 = [t for t in CUM2[1] if not any(n == "L" and r == 3 for n, r, u, l in t[1])]


def wedge(a, b):
    return np.einsum("pq,vw->pvqw", a, b) - np.einsum("pw,vq->pvqw", a, b)


def orbit(key, g1, l2, M):
    t = {("G", 1): g1, ("L", 2): l2}
    out = np.zeros((M,) * 6, dtype=complex)
    for s, sub, keys in ORB[key]:
        out = out + s * np.einsum(sub, *[t[k] for k in keys], optimize=True)
    return out


print("=== Parts 1 and 2: the trace identities on exact states ===")
print("%8s | %14s | %14s %14s"
      % ("(M,N)", "rank-2 resid", "rank-3 resid", "sum_y A form"))
for M, N in ((4, 2), (6, 3), (8, 4)):
    rng = np.random.default_rng(5)
    psi = fullci.random_state(M, N, rng)
    g1 = fullci.rdm(psi, M, 1)
    g2 = fullci.rdm(psi, M, 2)
    g3 = fullci.rdm(psi, M, 3)
    l2 = g2 - wedge(g1, g1)
    A = orbit("A", g1, l2, M)
    B = orbit("B", g1, l2, M)
    l3 = g3 - A - B
    r2 = np.abs(np.einsum("pxqx->pq", l2) + (g1 - g1 @ g1)).max()
    sA, sB = np.einsum("pvyqwy->pvqw", A), np.einsum("pvyqwy->pvqw", B)
    r3 = np.abs(np.einsum("pvyqwy->pvqw", l3) - ((N - 2) * g2 - sA - sB)).max()
    g1sq = g1 @ g1
    predA = N * wedge(g1, g1) - wedge(g1, g1sq) - wedge(g1sq, g1)
    print("%8s | %14.3e | %14.3e %14.3e"
          % ("(%d,%d)" % (M, N), r2, r3, np.abs(sA - predA).max()))
print()
print("rank-3 resid checks sum_y lambda_3 = (N-2) gamma_2 - sum_y A - sum_y B, which is the")
print("definition rearranged. sum_y A form checks the hand-derived closed form.")
print()

M, N = 6, 3
EYE = np.eye(M, dtype=complex)
rng = np.random.default_rng(7)
h, V = random_h(M, rng), random_V(M, rng)
psi = fullci.random_state(M, N, rng)
g1 = fullci.rdm(psi, M, 1)
g2 = fullci.rdm(psi, M, 2)
g3 = fullci.rdm(psi, M, 3)
l2 = g2 - wedge(g1, g1)
l3_exact = g3 - orbit("A", g1, l2, M) - orbit("B", g1, l2, M)


def rdot(l3, which):
    t = {("D", 1): EYE, ("G", 1): g1, ("L", 2): l2, ("L", 3): l3}
    d1 = assemble(CUM1[0], CUM1[1], h, V, t, "pq", M)
    if which == "all":
        d2 = assemble(CUM2[0], CUM2[1], h, V, t, "pvqw", M)
    elif which == "no_l3":
        d2 = assemble(CUM2[0], NO_L3, h, V, t, "pvqw", M)
    else:
        d2 = assemble([], L3_TERMS, h, V, t, "pvqw", M)
    dl2 = d2 - (wedge(d1, g1) + wedge(g1, d1))
    return (np.einsum("pxqx->pq", dl2) + d1 - (d1 @ g1 + g1 @ d1)), d1


print("=== Part 3: which terms break the rank-2 identity? ===")
print("   %d of %d two-body terms in the lambda_2 EOM carry lambda_3"
      % (len(L3_TERMS), len(CUM2[1])))
print()
print("%-26s | %13s" % ("lambda_3 used", "|Rdot|"))
for nm, l3 in (("exact", l3_exact),
               ("zero (lambda_3 = 0)", np.zeros((M,) * 6, dtype=complex)),
               ("MK_gamma30, -(A+B)", -(orbit("A", g1, l2, M) + orbit("B", g1, l2, M)))):
    R, _ = rdot(l3, "all")
    print("%-26s | %13.3e" % (nm, np.abs(R).max()))
print()
print("%-26s | %13s %13s" % ("split, exact lambda_3", "no-l3 terms", "l3 terms only"))
Rn, _ = rdot(l3_exact, "no_l3")
Rl, _ = rdot(l3_exact, "l3")
print("%-26s | %13.3e %13.3e" % ("", np.abs(Rn).max(), np.abs(Rl).max()))
print()
print("If the no-l3 part of Rdot is nonzero on its own but the total vanishes with the exact")
print("lambda_3, the identity is maintained by a CANCELLATION between the lambda_3 terms and")
print("the rest, not by the lambda_3 terms being individually harmless.")
