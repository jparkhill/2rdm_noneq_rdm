#!/usr/bin/env python3
"""One connected structure for lambda_3: the E-dressed pair chain. Is one coefficient universal?

MOTIVATION, from tools/connected_chains.py. Fitting lambda_3^exact inside the 4-dimensional
connected span returned coefficients whose ratio c(eta)/c(bare) was -2.000 at H3 and -1.998 at H4,
with the two anti(E)-dressed coefficients near zero. Since the two chains carry I and
eta_1 = I - gamma_1 on the contraction line, that combination is I - 2(I - gamma_1) = 2 gamma_1 - I
= 2E. So the oracle collapses onto a single structure. This file tests that directly.

THE STRUCTURE. With lambda_2[p v, q w] the 2-cumulant (p,v upper; q,w lower), E = gamma_1 - I/2
the particle-hole imbalance, and A3 the full antisymmetriser over the three upper and three lower
indices (36 signed permutations, divided by 36):

    CHAIN_E[p v x, q w y] = A3[ sum_{r,s} lambda_2[p v, q r] E[r,s] lambda_2[s x, w y] ]
                            then Hermitised in the pair sense, X -> (X + X^dag)/2

It is CONNECTED: the two lambda_2 factors share the summed index through E, so for
non-interacting fragments block-diagonality of lambda_2 forces all six external indices into one
fragment. Hence leak is exactly zero, verified in tools/connected_chains.py for every member of
this family. It is a polynomial in gamma_1 and lambda_2 with no inverse of anything, so it is
Lipschitz with a fixed constant. It is O(V^2), which is the order of lambda_3.

THE TWO WAYS TO FIX THE SINGLE COEFFICIENT c, which must not be confused:
    ORACLE   c from least squares against lambda_3^exact. An upper bound on the structure's
             reach; NOT a closure, since it has seen the answer.
    CLOSURE  c from least squares on the partial trace alone,
                 Tr_3 (c CHAIN_E) = anti(E, lambda_2),   Tr_3 T[pv,qw] = sum_y T[p v y, q w y]
             This never looks at lambda_3, so it IS a closure with no free parameters. With
             C(M,2)^2 equations for one unknown it is massively overdetermined.

STATES. Pure FCI ground states and thermal mixtures at kT = 0.5 hartree, since CLAUDE.md requires
a mixed row and the target regime is non-equilibrium and mixed. Ground truth qdiss.eom.fullci in
the N-electron sector.
"""
import itertools
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
       (H4_CHAIN, "sto-3g", 4, 0, 0, "H4 linear chain"),
       (LIH, "sto-3g", 4, 0, 0, "LiH"))
KT = 0.5

_P3 = [(p, (-1) ** sum(1 for i in range(3) for j in range(i + 1, 3) if p[i] > p[j]))
       for p in itertools.permutations(range(3))]


def A3(T):
    out = np.zeros_like(T)
    for pu, su in _P3:
        for pl, sl in _P3:
            out += (su * sl) * T.transpose(tuple(pu) + tuple(3 + i for i in pl))
    return out / 36.0


def chain_E(g1, l2, M):
    E = g1 - 0.5 * np.eye(M, dtype=complex)
    X = A3(np.einsum("pvqr,rs,sxwy->pvxqwy", l2, E, l2, optimize=True))
    return 0.5 * (X + np.conj(X.transpose(3, 4, 5, 0, 1, 2)))


print("CHAIN_E = A3[ lambda_2 . E . lambda_2 ], Hermitised.  Connected, so leak is exactly zero.")
print("oracle  c fitted to lambda_3^exact (upper bound, not a closure)")
print("closure c fitted to Tr_3 = anti(E, lambda_2) only (no free parameters)")
print()
print("%-18s %-9s %5s | %8s %8s | %8s %8s | %9s %9s"
      % ("system", "state", "M", "c orac", "err", "c clos", "err", "trace res", "|lam3|"))
for atom, basis, N, chg, spin, label in SYS:
    mol, mf, M, h, V, enuc = molecule(atom, basis, chg, spin)
    dets = sector(M, N)
    es, vec = np.linalg.eigh(build_H_sector(h, V, M, dets))
    dm = np.zeros((len(dets),) * 2, dtype=complex)
    dm[0, 0] = 1.0
    w = np.exp(-(es - es[0]).real / KT)
    w /= w.sum()
    for sname, dmat in (("pure", vec @ dm @ vec.conj().T),
                        ("kT=%.1f" % KT, vec @ np.diag(w.astype(complex)) @ vec.conj().T)):
        g1, g2, g3 = rdms(dmat, M, N, dets)
        l2 = az.cumulant(2, g2, {("G", 1): g1})
        l3 = az.cumulant(3, g3, {("G", 1): g1, ("L", 2): l2})
        C = chain_E(g1, l2, M)
        S = az.trace_target(3, g1, l2)
        TrC = np.einsum("pvyqwy->pvqw", C)
        co = float(np.real(np.vdot(C, l3) / np.vdot(C, C)))
        cc = float(np.real(np.vdot(TrC, S) / np.vdot(TrC, TrC)))
        print("%-18s %-9s %5d | %8.3f %8.4f | %8.3f %8.4f | %9.2e %9.5f"
              % (label, sname, M, co, fro(co * C - l3) / fro(l3),
                 cc, fro(cc * C - l3) / fro(l3),
                 az.trace_residual(cc * C, S, 3), fro(l3)))
print()
print("baselines on the same states (tools/ansatz_hierarchy.py): null 1.0000,")
print("min-norm trace 0.26-0.82, Mueller+Lackner 1.58-2.73.  All of those are non-extensive")
print("except null; CHAIN_E is extensive by construction.")
