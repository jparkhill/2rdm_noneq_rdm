#!/usr/bin/env python3
"""The connected-contraction closure for lambda_3, across systems and states, against baselines.

THE ANSATZ. lambda_3 is expanded in the CONNECTED span: antisymmetrised single contractions of
two 2-cumulants, optionally with gamma_1 or eta_1 = I - gamma_1 on the contraction line and with
E = gamma_1 - I/2 contracted into the free slots. Every member is a connected diagram, so for
non-interacting fragments block-diagonality of lambda_2 forces all six external indices into one
fragment and the cross-fragment leak is EXACTLY zero (verified for all 12 generators in
tools/connected_chains.py). Every member is a polynomial in gamma_1 and lambda_2 with no inverse
of any state-dependent operator, so screen 1 (Lipschitz) is passed structurally. Every member is
O(V^2), the order of lambda_3.

    lambda_2[p v, q w]   the 2-cumulant, p,v upper and q,w lower
    E = gamma_1 - I/2    particle-hole imbalance
    A3                   antisymmetriser over the three upper and three lower indices, 36 signed
                         permutations divided by 36
    Tr_3 T[pv,qw]        = sum_y T[p v y, q w y]
    anti(E, X)           E contracted into each of X's slots, summed

HOW THE COEFFICIENTS ARE FIXED, and the distinction that matters:
    ORACLE   least squares against lambda_3^exact. An UPPER BOUND on the span's reach. Not a
             closure: it has seen the answer.
    CLOSURE  least squares on the partial trace alone, Tr_3 lambda_3 = anti(E, lambda_2). This
             never looks at lambda_3, so it IS a closure. It is solved at each state from the
             current gamma_1 and lambda_2, which is what makes it state-adaptive rather than a
             fitted constant; the condition number of that small solve is reported, because an
             ill-conditioned solve is a Lipschitz hazard.

BASELINES, all from tools/ansatz_hierarchy.py on the same states: null (lambda_3 = 0), the
minimum-norm partial-trace inverse, and Mueller arXiv:2405.05093 Eq. 14 with the Lackner Eq.
43/44/45 correction. Of those only null is extensive.

Ground truth qdiss.eom.fullci in the N-electron sector.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from ansatz_hierarchy import fro, rank3_set, rdms
from molecular_2rdm import build_H_sector, molecule, sector
from qdiss import ansatz as az

H3_BENT = "H 0.00 0.00 0.0; H 1.00 0.00 0.0; H 0.35 0.90 0.0"
H4_TRAP = "H 0.00 0.00 0.0; H 0.80 0.00 0.0; H 1.05 0.72 0.0; H -0.15 0.72 0.0"
H4_CHAIN = "H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0"
LIH = "Li 0 0 0; H 0 0 1.6"
SYS = ((H3_BENT, "sto-3g", 3, 0, 1, "H3 bent"),
       (H4_TRAP, "sto-3g", 4, 0, 0, "H4 trapezoid"),
       (H4_CHAIN, "sto-3g", 4, 0, 0, "H4 linear"),
       (LIH, "sto-3g", 4, 0, 0, "LiH"))
KTS = (0.5, 0.1)

print("connected-span lambda_3 closure vs baselines.  err = |lambda_3^rec - lambda_3|_F / |lambda_3|_F")
print("the connected span is EXTENSIVE (leak exactly 0) and Lipschitz; null is too; the other")
print("two baselines are neither.")
print()
print("%-14s %-8s | %6s %6s %8s | %7s %7s %7s"
      % ("system", "state", "oracle", "closur", "cond", "null", "trace", "Mue+La"))
for atom, basis, N, chg, spin, label in SYS:
    mol, mf, M, h, V, enuc = molecule(atom, basis, chg, spin)
    dets = sector(M, N)
    es, vec = np.linalg.eigh(build_H_sector(h, V, M, dets))
    dm = np.zeros((len(dets),) * 2, dtype=complex)
    dm[0, 0] = 1.0
    states = [("pure", vec @ dm @ vec.conj().T)]
    for kt in KTS:
        w = np.exp(-(es - es[0]).real / kt)
        w /= w.sum()
        states.append(("kT=%.1f" % kt, vec @ np.diag(w.astype(complex)) @ vec.conj().T))
    for sname, dmat in states:
        g1, g2, g3 = rdms(dmat, M, N, dets)
        l2 = az.cumulant(2, g2, {("G", 1): g1})
        l3 = az.cumulant(3, g3, {("G", 1): g1, ("L", 2): l2})
        base = az.connected_basis(g1, l2, M)
        B = np.array([T.reshape(-1) for _, T in base])
        S = az.trace_target(3, g1, l2)
        co, *_ = np.linalg.lstsq(B.T, l3.reshape(-1), rcond=None)
        Tr = np.array([np.einsum("pvyqwy->pvqw", b.reshape(l3.shape)).reshape(-1) for b in B])
        cc, *_ = np.linalg.lstsq(Tr.T, S.reshape(-1), rcond=None)
        sv = np.linalg.svd(Tr.T, compute_uv=False)
        sv = sv[sv > 1e-12 * sv[0]]
        rk3 = rank3_set(g1, l2, N, M)
        print("%-14s %-8s | %6.4f %6.4f %8.1e | %7.4f %7.4f %7.4f"
              % (label, sname,
                 fro((co @ B).reshape(l3.shape) - l3) / fro(l3),
                 fro((cc @ B).reshape(l3.shape) - l3) / fro(l3),
                 sv[0] / sv[-1],
                 fro(rk3["null"] - l3) / fro(l3),
                 fro(rk3["trace"] - l3) / fro(l3),
                 fro(rk3["scaled+trace"] - l3) / fro(l3)))
