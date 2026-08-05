#!/usr/bin/env python3
"""Physical spin-orbital Hamiltonians from PySCF, validated against FCI.

Everything downstream (poles, functionals) rests on these integrals, so the only thing this
file does is build them in the convention of qdiss.eom.hamiltonian,

    H = sum_pq h[p,q] a+_p a_q  +  (1/4) sum V[r,t,s,u] a+_r a+_t a_u a_s ,  V = <rt||su>

and check that the Fock-space matrix built by qdiss.eom.fullci.build_H reproduces the PySCF
FCI total energy. If that fails, nothing after it means anything.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
from pyscf import ao2mo, fci, gto, scf
from qdiss.eom import fullci


def spin_orbital_H(atom, basis, charge=0, verbose=False):
    mol = gto.M(atom=atom, basis=basis, charge=charge, verbose=0)
    mf = scf.RHF(mol).run()
    C = mf.mo_coeff
    nmo = C.shape[1]
    hcore = C.T @ mf.get_hcore() @ C
    eri = ao2mo.restore(1, ao2mo.kernel(mol, C), nmo)
    M = 2 * nmo
    h = np.zeros((M, M), dtype=complex)
    for p in range(nmo):
        for q in range(nmo):
            h[2 * p, 2 * q] = hcore[p, q]
            h[2 * p + 1, 2 * q + 1] = hcore[p, q]
    phys = np.zeros((M,) * 4)
    for p in range(nmo):
        for q in range(nmo):
            for r in range(nmo):
                for s in range(nmo):
                    val = eri[p, r, q, s]
                    for sp in (0, 1):
                        for sq in (0, 1):
                            phys[2 * p + sp, 2 * q + sq, 2 * r + sp, 2 * s + sq] += val
    V = (phys - np.transpose(phys, (0, 1, 3, 2))).astype(complex)
    return mol, mf, M, h, V, float(mol.energy_nuc())


CASES = (("H 0 0 0; H 0 0 0.74", "6-31g", 2, 0),
         ("H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0", "sto-3g", 4, 0),
         ("He 0 0 0; H 0 0 0.9", "sto-3g", 2, 1))

print("%-34s %6s %5s %5s | %14s %14s %11s"
      % ("system", "basis", "M", "N", "our ground E", "PySCF FCI", "diff"))
for atom, basis, N, chg in CASES:
    mol, mf, M, h, V, enuc = spin_orbital_H(atom, basis, chg)
    if M > 12:
        print("%-34s %6s %5d %5d | skipped, Fock space too large" % (atom, basis, M, N))
        continue
    Hm = fullci.build_H(h, V, M)
    dim = 1 << M
    popc = np.array([bin(s).count("1") for s in range(dim)])
    sub = np.where(popc == N)[0]
    e0 = np.linalg.eigvalsh(Hm[np.ix_(sub, sub)])[0].real + enuc
    cis = fci.FCI(mf)
    efci = cis.kernel()[0]
    print("%-34s %6s %5d %5d | %14.9f %14.9f %11.2e"
          % (atom.replace(";", ","), basis, M, N, e0, efci, abs(e0 - efci)))
