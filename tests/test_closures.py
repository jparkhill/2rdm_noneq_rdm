#!/usr/bin/env python3
"""Validation of qdiss.closures against the properties each reconstruction is meant to have.

Checks, on real PySCF Hamiltonians at the exact FCI ground-state gamma_1:
  1. antisymmetry of every reconstruction in both index pairs
  2. the partial-trace identity, before and after trace_correct
  3. that hf reproduces TDHF exactly, i[(h+W)^T, gamma_1]
  4. that occupation-diagonal reconstructions move no occupation, and mp2 does
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from physical_hamiltonian import spin_orbital_H
from qdiss import closures
from qdiss.eom import fullci, gwt
from qdiss.eom.evaluate import assemble

CUM1 = gwt.terms_cumulant(1)

for atom, basis, N, chg, label in (
        ("H 0 0 0; H 0 0 0.74", "6-31g", 2, 0, "H2 / 6-31g"),
        ("H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0", "sto-3g", 4, 0, "H4 / sto-3g")):
    mol, mf, M, h, V, enuc = spin_orbital_H(atom, basis, chg)
    Hm = fullci.build_H(h, V, M)
    popc = np.array([bin(s).count("1") for s in range(1 << M)])
    sub = np.where(popc == N)[0]
    es, vec = np.linalg.eigh(Hm[np.ix_(sub, sub)])
    psi0 = np.zeros(1 << M, dtype=complex)
    psi0[sub] = vec[:, 0]
    g1 = fullci.rdm(psi0, M, 1)
    rr = np.random.default_rng(5)
    kick = np.zeros(1 << M, dtype=complex)
    kick[sub] = rr.normal(size=len(sub)) + 1j * rr.normal(size=len(sub))
    psik = psi0 + 0.3 * kick / np.linalg.norm(kick)
    psik /= np.linalg.norm(psik)
    g1k = fullci.rdm(psik, M, 1)
    natk, Unok = np.linalg.eigh(0.5 * (g1k + g1k.conj().T))
    nat, Uno = np.linalg.eigh(0.5 * (g1 + g1.conj().T))
    W = np.einsum("prqs,rs->pq", V, g1)
    print("=== %s   M=%d N=%d ===" % (label, M, N))
    print("%-9s | %9s %9s | %11s %11s | %11s"
          % ("recon", "antisym u", "antisym l", "trace raw", "trace tc", "max|ndot| kicked"))
    for name in closures.NAMES:
        l2 = closures.lambda2(g1, name, h, V)
        au = np.abs(l2 + np.transpose(l2, (1, 0, 2, 3))).max()
        al = np.abs(l2 + np.transpose(l2, (0, 1, 3, 2))).max()
        traw = closures.trace_residual(l2, g1)
        ttc = closures.trace_residual(closures.trace_correct(l2, g1), g1)
        l2k = closures.lambda2(g1k, name, h, V)
        t = {("D", 1): np.eye(M, dtype=complex), ("G", 1): g1k, ("L", 2): l2k}
        d = assemble(CUM1[0], CUM1[1], h, V, t, "pq", M)
        nd = np.abs(np.real(np.diag(Unok.conj().T @ d @ Unok))).max()
        print("%-9s | %9.1e %9.1e | %11.3e %11.3e | %11.3e"
              % (name, au, al, traw, ttc, nd))
    t = {("D", 1): np.eye(M, dtype=complex), ("G", 1): g1,
         ("L", 2): closures.lambda2(g1, "hf")}
    d = assemble(CUM1[0], CUM1[1], h, V, t, "pq", M)
    f = (h + W).conj()
    print("   hf equals TDHF i[(h+W)^T,g1] : %.3e"
          % np.abs(d - 1j * (f @ g1 - g1 @ f)).max())
    print()
