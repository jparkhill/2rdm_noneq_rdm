#!/usr/bin/env python3
"""Extensivity (CLAUDE.md screen 2) of EVERY ansatz at ranks 2 and 3, with a real second fragment.

WHY THIS SUPERSEDES THE ZERO-PADDING TEST in tools/ansatz_screens.py. Padding with empty orbitals
leaves N unchanged, so it cannot detect a reconstruction whose coefficients depend on the particle
number -- which is exactly how Mueller et al. arXiv:2405.05093 Eq. 14 is built, its prefactors
being (N-2)/N and 2(N-1)(N-2)/N^2. Here the spectator is a real fragment carrying electrons, so
N genuinely changes and the N-dependence is exposed.

THE CONSTRUCTION, exact and free of any basis-matching question.
    Fragment A: H3 bent / sto-3g, M_A = 6 spin orbitals, N_A = 3, at its FCI ground state.
    Fragment B: H2 / sto-3g at 0.9 Angstrom, M_B = 4, N_B = 2, at its FCI ground state.
    The composite is the PRODUCT state |Psi_A> tensor |Psi_B> in the M = 10 Fock space, built by
    multiplying determinant amplitudes; since every A index is below every B index the fermionic
    reordering sign is +1, and this is checked rather than assumed by comparing the composite's
    exact cumulants against the block direct sums below.

WHAT EXTENSIVITY DEMANDS, for a composite of two non-interacting fragments:
    gamma_1^AB = gamma_1^A (+) gamma_1^B                 block direct sum, verified here
    lambda_k^AB = lambda_k^A (+) lambda_k^B              BLOCK ADDITIVE with EXACTLY zero
                                                          cross-fragment elements, verified here
Note this is a statement about CUMULANTS, not moments: gamma_2^AB does have cross-fragment
elements (one particle in A, one in B) and must, since gamma_2 = gamma_1 wedge gamma_1 + lambda_2.

COLUMNS, for a reconstruction fed the composite gamma_1 and lambda_2 with the composite N = 5:
    leak      max |lambda_k^rec| over index tuples touching BOTH fragments, divided by
              max |lambda_k^rec| over tuples wholly inside A. Must be exactly 0.
    drift     max change of the A block against the same ansatz run on fragment A alone with
              N = 3, divided by the magnitude of that A-block answer. Must be exactly 0.
    exact     the same two quantities for the exact lambda_k, which are the control and must
              come out at round-off.

Symbols as in qdiss.ansatz. Ground truth qdiss.eom.fullci.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from ansatz_hierarchy import fro, rank2_set, rank3_set
from molecular_2rdm import build_H_sector, molecule, sector
from qdiss import ansatz as az
from qdiss.eom import fullci

H3_BENT = "H 0.00 0.00 0.0; H 1.00 0.00 0.0; H 0.35 0.90 0.0"
H2 = "H 0 0 0; H 0 0 0.9"
MA, NA = 6, 3
MB, NB = 4, 2
M, N = MA + MB, NA + NB


def ground(atom, basis, chg, spin, Mx, Nx):
    mol, mf, Mchk, h, V, enuc = molecule(atom, basis, chg, spin)
    assert Mchk == Mx, "expected M = %d, got %d" % (Mx, Mchk)
    dets = sector(Mx, Nx)
    es, vec = np.linalg.eigh(build_H_sector(h, V, Mx, dets))
    psi = np.zeros(1 << Mx, dtype=complex)
    psi[dets] = vec[:, 0]
    return psi, es[0].real + enuc


psiA, eA = ground(H3_BENT, "sto-3g", 0, 1, MA, NA)
psiB, eB = ground(H2, "sto-3g", 0, 0, MB, NB)
# product state: all A labels are below all B labels, so |d_A>|d_B> = |d_A | d_B<<MA> with sign +1
psi = np.zeros(1 << M, dtype=complex)
for dA in np.nonzero(psiA)[0]:
    for dB in np.nonzero(psiB)[0]:
        psi[dA | (dB << MA)] = psiA[dA] * psiB[dB]
print("composite = H3 bent (M=6,N=3) x H2/0.9A (M=4,N=2)  ->  M=%d N=%d" % (M, N))
print("  E_A %.9f   E_B %.9f   norm of product state %.12f"
      % (eA, eB, float(np.linalg.norm(psi))))


def cums(p, Mx):
    g1 = fullci.rdm(p, Mx, 1)
    g2 = fullci.rdm(p, Mx, 2)
    g3 = fullci.rdm(p, Mx, 3)
    l2 = az.cumulant(2, g2, {("G", 1): g1})
    l3 = az.cumulant(3, g3, {("G", 1): g1, ("L", 2): l2})
    return g1, l2, l3


g1A, l2A, l3A = cums(psiA, MA)
g1, l2, l3 = cums(psi, M)


def blocks(T, k):
    """(max inside A, max touching both fragments) for a 2k-index tensor on M orbitals."""
    A = np.abs(T)
    inA = float(A[tuple(slice(0, MA) for _ in range(2 * k))].max())
    mask = np.zeros((M,) * (2 * k), dtype=bool)
    for ax in range(2 * k):
        sl = [slice(None)] * (2 * k)
        sl[ax] = slice(MA, M)
        mask[tuple(sl)] = True
    inside = np.zeros_like(mask)
    inside[tuple(slice(0, MA) for _ in range(2 * k))] = True
    allB = np.ones_like(mask)
    allB[tuple(slice(MA, M) for _ in range(2 * k))] = False
    cross = mask & allB & ~inside
    return inA, float(A[cross].max())


print()
print("CONTROL: the exact composite cumulants must be block additive")
for k, lam, lamA in ((2, l2, l2A), (3, l3, l3A)):
    inA, cross = blocks(lam, k)
    dA = float(np.abs(lam[tuple(slice(0, MA) for _ in range(2 * k))] - lamA).max())
    print("  rank %d  |A block| %.6f   cross-fragment %.2e   |A block - fragment alone| %.2e"
          % (k, inA, cross, dA))
g1blk = float(np.abs(g1[:MA, :MA] - g1A).max())
print("  gamma_1 block additive: |A block - fragment alone| %.2e   cross %.2e"
      % (g1blk, float(np.abs(g1[:MA, MA:]).max())))

print()
print("SCREEN 2 for every ansatz.  leak and drift must both be exactly 0.")
print("  %-4s %-20s %12s %12s %12s" % ("rank", "ansatz", "|A block|", "leak", "drift"))
for k, setter_big, setter_small, lam, lamA in (
        (2, rank2_set(g1, N, M), rank2_set(g1A, NA, MA), l2, l2A),
        (3, rank3_set(g1, l2, N, M), rank3_set(g1A, l2A, NA, MA), l3, l3A)):
    for nm in setter_big:
        big, small = setter_big[nm], setter_small[nm]
        inA, cross = blocks(big, k)
        blk = big[tuple(slice(0, MA) for _ in range(2 * k))]
        ref = max(float(np.abs(small).max()), 1e-30)
        print("  %-4d %-20s %12.6f %12.3e %12.3e"
              % (k, nm, inA, cross / max(inA, 1e-30),
                 float(np.abs(blk - small).max()) / ref))
    inA, cross = blocks(lam, k)
    print("  %-4d %-20s %12.6f %12.3e %12.3e"
          % (k, "exact (control)", inA, cross / inA,
             float(np.abs(lam[tuple(slice(0, MA) for _ in range(2 * k))] - lamA).max())
             / float(np.abs(lamA).max())))
    print()
