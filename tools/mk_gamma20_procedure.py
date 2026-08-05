#!/usr/bin/env python3
"""MK_gamma20 by the stated gamma-zero procedure.

Procedure, as specified, and to be used for every gamma-zero closure from here on:

  1. Start from the MKNO 1-RDM EOM with NO approximation. In cumulant form
     (qdiss.eom.gwt.terms_cumulant(1)) it is a function of gamma_1 and lambda_2 only:
     2 + 12 terms with factors D1 x14, G1 x18, L2 x4.
  2. Solve the identity gamma_2 = 0 for the top cumulant. Since
     gamma_2 = gamma_1 wedge gamma_1 + lambda_2, this gives

         lambda_2 := -(gamma_1 wedge gamma_1)

     an explicit product of gamma_1's, not a tensor set to zero.
  3. Substitute that expression back into the unapproximated MKNO 1-RDM EOM and propagate
     gamma_1 alone.

This mirrors MK_gamma30 exactly: there the unapproximated MKNO 2-RDM EOM is a function of
gamma_1, lambda_2 and lambda_3, the identity gamma_3 = 0 is solved for lambda_3 = -(A + B),
and (gamma_1, lambda_2) are propagated.

The previous attempt reported "case (c) gamma_2 = 0" by passing a zero tensor to the
BARE-moment EOM. That is a different route to the same point in variable space, and this
file measures whether the two agree rather than asserting it. Also compared: TDHF, obtained
from the same unapproximated MKNO 1-RDM EOM with lambda_2 := 0, and the exact equation.

All term lists validated against qdiss.eom.fullci at 1e-15 by qdiss.eom.validate.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
from qdiss.eom import fullci, gwt, wick
from qdiss.eom.evaluate import assemble
from qdiss.eom.hamiltonian import random_V, random_h

M, N = 6, 3
EYE = np.eye(M, dtype=complex)
CUM1 = gwt.terms_cumulant(1)
inv = {}
for s, f in list(CUM1[0]) + list(CUM1[1]):
    for n, r, u, l in f:
        inv[(n, r)] = inv.get((n, r), 0) + 1
print("unapproximated MKNO 1-RDM EOM: %d + %d terms, factors %s"
      % (len(CUM1[0]), len(CUM1[1]),
         ", ".join("%s%d x%d" % (k[0], k[1], v) for k, v in sorted(inv.items()))))
print()


def wedge(a, b):
    return np.einsum("pq,vw->pvqw", a, b) - np.einsum("pw,vq->pvqw", a, b)


rng = np.random.default_rng(21)
h, V = random_h(M, rng), random_V(M, rng)
psi = fullci.random_state(M, N, rng)
g1 = fullci.rdm(psi, M, 1)
g2 = fullci.rdm(psi, M, 2)
l2_exact = g2 - wedge(g1, g1)


def mkno_d1(l2):
    t = {("D", 1): EYE, ("G", 1): g1, ("L", 2): l2}
    return assemble(CUM1[0], CUM1[1], h, V, t, "pq", M)


def bare_d1(g2u):
    return wick.gamma_dot(1, h, V, {("D", 1): EYE, ("G", 1): g1, ("G", 2): g2u}, M)


d_proc = mkno_d1(-wedge(g1, g1))
d_zero = bare_d1(np.zeros_like(g2))
d_free = 1j * (h.conj() @ g1 - g1 @ h.conj())
d_tdhf = mkno_d1(np.zeros_like(g2))
d_exact = mkno_d1(l2_exact)

print("=== the stated procedure, and what it equals ===")
print("   MK_gamma20 via MKNO EOM with lambda_2 := -(g1 wedge g1)   |d| = %.6f"
      % np.abs(d_proc).max())
print("      vs BARE EOM with gamma_2 := 0                          %.3e"
      % np.abs(d_proc - d_zero).max())
print("      vs i[h^T, gamma_1]  (free evolution)                   %.3e"
      % np.abs(d_proc - d_free).max())
print("      vs TDHF  (same EOM with lambda_2 := 0)                 %.3e"
      % np.abs(d_proc - d_tdhf).max())
print()
print("   TDHF via MKNO EOM with lambda_2 := 0                      |d| = %.6f"
      % np.abs(d_tdhf).max())
W = np.einsum("prqs,rs->pq", V, g1)
print("      vs i[(h+W)^T, gamma_1]                                 %.3e"
      % np.abs(d_tdhf - 1j * ((h + W).conj() @ g1 - g1 @ (h + W).conj())).max())
print("   exact MKNO EOM                                            |d| = %.6f"
      % np.abs(d_exact).max())
print("      vs FullCI                                              %.3e"
      % np.abs(d_exact - fullci.rdm_dot(psi, fullci.build_H(h, V, M), M, 1)).max())
print()

print("=== does the substitution cancel the interaction symbolically? ===")
n_l2 = sum(1 for s, f in CUM1[1] if any(n == "L" and r == 2 for n, r, u, l in f))
d_h_only = assemble(CUM1[0], [], h, V, {("D", 1): EYE, ("G", 1): g1,
                                        ("L", 2): -wedge(g1, g1)}, "pq", M)
d_v_only = assemble([], CUM1[1], h, V, {("D", 1): EYE, ("G", 1): g1,
                                        ("L", 2): -wedge(g1, g1)}, "pq", M)
print("   two-body terms carrying lambda_2 : %d of %d" % (n_l2, len(CUM1[1])))
print("   one-body piece under the substitution      |d| = %.6f" % np.abs(d_h_only).max())
print("   two-body piece under the substitution      |d| = %.3e" % np.abs(d_v_only).max())
print("   so the V contribution cancels to that residual once lambda_2 = -(g1 wedge g1).")
