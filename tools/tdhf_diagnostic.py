#!/usr/bin/env python3
"""Diagnose the pole-extraction path by routing TDHF through it.

    TDHF(pyscf)   pyscf.tdscf.TDHF on the UHF reference. Reference implementation.
    tdhf(claude)  the MKNO 1-RDM EOM, qdiss.eom.gwt.terms_cumulant(1), a function of gamma_1 and
                  lambda_2 only, with the two-particle cumulant TRUNCATED, lambda_2 := 0.
                  Linearised at its own HF fixed point by the same finite-difference Jacobian
                  and the same imaginary-eigenvalue readout used for the gamma_2 closures.

Note these are different SPECIES of approximation. tdhf(claude) truncates a cumulant to zero;
lam3_0 and muller RECONSTRUCT gamma_3 from lower RDMs. The point of the row is not to compete
but to test the extraction: if tdhf(claude) reproduces TDHF(pyscf), the Jacobian-plus-imaginary-
eigenvalue readout is sound and any trouble lies in the closures; if it does not, the readout is
the problem.

The gamma_2-flow rows are carried over for context and are extracted at the FCI reference over
the pair basis, which is a larger space with many more modes -- a caveat, not a result.
"""
import sys
import numpy as np
sys.path.insert(0, "."); sys.path.insert(0, "tools")
from pyscf import gto, tdscf
from molecular_2rdm import build_H_sector, sector
from qdiss.eom import fullci, gwt
from qdiss.eom.evaluate import assemble
from unrestricted_spectra import spin_orbital_uhf, uhf_stable

CUM1 = gwt.terms_cumulant(1)
inv = {}
for s, f in list(CUM1[0]) + list(CUM1[1]):
    for n, r, u, l in f:
        inv[(n, r)] = inv.get((n, r), 0) + 1
print("MKNO 1-RDM EOM: %d + %d terms, factors %s"
      % (len(CUM1[0]), len(CUM1[1]),
         ", ".join("%s%d x%d" % (k[0], k[1], v) for k, v in sorted(inv.items()))))

H4_TRAP = "H 0.00 0.00 0.0; H 0.80 0.00 0.0; H 1.05 0.72 0.0; H -0.15 0.72 0.0"
mol = gto.M(atom=H4_TRAP, basis="sto-3g", verbose=0)
N = mol.nelectron
mfu = uhf_stable(mol)
M, h, V = spin_orbital_uhf(mol, mfu)[0:3]
es, vec = np.linalg.eigh(build_H_sector(h, V, M, sector(M, N)))
exc = np.sort((es - es[0]).real); exc = exc[exc > 1e-8][:10]

td = tdscf.TDHF(mfu); td.nstates = 14; td.verbose = 0; td.kernel()
tdp = np.sort(np.asarray(td.e).real)[:10]

EYE = np.eye(M, dtype=complex)
ZERO4 = np.zeros((M,) * 4, dtype=complex)


def d1_mkno_lam2zero(g1):
    """MKNO 1-RDM EOM with the two-particle cumulant truncated to zero."""
    return assemble(CUM1[0], CUM1[1], h, V,
                    {("D", 1): EYE, ("G", 1): g1, ("L", 2): ZERO4}, "pq", M)


def fock(g):
    return h + np.einsum("prqs,rs->pq", V, g)


g = np.zeros((M, M), dtype=complex)
for _ in range(800):
    _, vv = np.linalg.eigh(0.5 * (fock(g).conj() + fock(g).conj().conj().T))
    new = vv[:, :N] @ vv[:, :N].conj().T
    if np.abs(new - g).max() < 1e-13:
        g = new; break
    g = 0.6 * new + 0.4 * g
f0 = d1_mkno_lam2zero(g)
print("tdhf(claude) fixed point: |gamma_1-dot| = %.2e   |g^2-g| = %.2e"
      % (np.abs(f0).max(), np.abs(g @ g - g).max()))

B = []
for i in range(M):
    E = np.zeros((M, M), dtype=complex); E[i, i] = 1.0; B.append(E)
for i in range(M):
    for j in range(i + 1, M):
        E = np.zeros((M, M), dtype=complex); E[i, j] = 1.0; E[j, i] = 1.0; B.append(E)
        E = np.zeros((M, M), dtype=complex); E[i, j] = 1j; E[j, i] = -1j; B.append(E)
n = len(B); J = np.zeros((n, n)); eps = 1e-6
for k in range(n):
    dp = (d1_mkno_lam2zero(g + eps * B[k]) - f0) / eps
    for m in range(n):
        J[m, k] = float(np.real(np.trace(B[m].conj().T @ dp))) / float(
            np.real(np.trace(B[m].conj().T @ B[m])))
w = np.linalg.eigvals(J)
tdc = np.sort(np.abs(w.imag[w.imag > 1e-6]))[:10]

print()
print("Excitation energies, Hartree, sorted ascending. UHF reference, H4 asymmetric trapezoid.")
print("%-14s | %s" % ("method", "  ".join("%7d" % (i + 1) for i in range(10))))
for k, a in (("FCI", exc), ("TDHF(pyscf)", tdp), ("tdhf(claude)", tdc)):
    print("%-14s | %s" % (k, "  ".join("%7.4f" % x for x in a)))
print()
nmin = min(len(tdp), len(tdc))
print("tdhf(claude) vs TDHF(pyscf), elementwise over the first %d: max |diff| = %.3e"
      % (nmin, np.abs(tdc[:nmin] - tdp[:nmin]).max()))
