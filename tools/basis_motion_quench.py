#!/usr/bin/env python3
r"""Does the auxiliary basis need memory? Same test, in the regimes where one-shot closures FAIL.

WHY THIS RERUN. tools/basis_motion.py ran the test on a dipole-kicked ground state with eps = 0.05.
The coefficients drifted 2 percent over t <= 1: a near-adiabatic trajectory, precisely where one-shot
closures already work (0.03-0.12 relative error at pure states). It therefore said nothing about
whether the ancilla is needed. The memory expansion is meant to matter where the one-shot closures
FAILED, which the earlier measurements localise exactly:

  - low rank(rho): at Neff ~ 2-3 the connected span's ORACLE error was 0.51-0.61, i.e. the span itself
    failed, not just its coefficients;
  - violent dynamics: the regime Mueller et al. probe with a quench, where their closure degrades.

THREE REGIMES, all propagated EXACTLY so lambda_3(t) is known:

  A  adiabatic   ground state of H, dipole-kicked, propagated under H          (control, as before)
  B  quench      ground state of H(compressed geometry), propagated under H(stretched geometry).
                 The state is far from any eigenstate of the propagating Hamiltonian, so the dynamics
                 is violent density sloshing rather than a small oscillation.
  C  low rank    equal mixture of the two lowest eigenstates of H(compressed), propagated under
                 H(stretched). rank(rho) = 2, so Neff = 2 -- the regime where the connected span's
                 own reach collapsed.

REPORTED, at sample times: how far the state has moved (|gamma_2(t) - gamma_2(0)|, so "violent" is
quantified rather than asserted); the span error of the moving connected basis; the basis-motion
fraction |sum_a c_a S_a-dot| / |lambda_3-dot|; and the coefficient drift |c(t) - c(0)| / |c(0)|, which
is the quantity that decides whether an ancilla is needed at all.

  c drift small  -> lambda_3 is a fixed function of the instantaneous (gamma_1, lambda_2); the one-shot
                    connected closure is the right answer and no ancilla is needed.
  c drift O(1)   -> the coefficients carry genuine memory and the auxiliary EOM is justified.

Ground truth qdiss.eom.fullci. Symbols as in qdiss.ansatz.
"""
import sys
import numpy as np
import scipy.linalg as sla
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from ansatz_hierarchy import rdms
from molecular_2rdm import build_H_sector, molecule, sector
from qdiss import ansatz as az
from qdiss.eom import fullci

COMPRESSED = "H 0.00 0.00 0.0; H 0.60 0.00 0.0; H 0.79 0.54 0.0; H -0.11 0.54 0.0"
STRETCHED = "H 0.00 0.00 0.0; H 1.20 0.00 0.0; H 1.58 1.08 0.0; H -0.22 1.08 0.0"
N, EPS, DT = 4, 0.05, 1e-5
TIMES = (0.0, 0.25, 0.5, 1.0, 2.0)


def setup(atom):
    mol, mf, M, h, V, enuc = molecule(atom, "sto-3g", 0, 0)
    dets = sector(M, N)
    es, vec = np.linalg.eigh(build_H_sector(h, V, M, dets))
    return mol, mf, M, h, V, dets, es, vec


molA, mfA, M, hA, VA, dets, esA, vecA = setup(COMPRESSED)
molB, mfB, M2, hB, VB, dets2, esB, vecB = setup(STRETCHED)
assert M == M2 and dets == dets2

dxm = mfA.mo_coeff.T @ molA.intor("int1e_r")[0] @ mfA.mo_coeff
assert np.abs(dxm).max() > 1e-8, "null dipole"
d = np.zeros((M, M), dtype=complex)
for p in range(dxm.shape[0]):
    for q in range(dxm.shape[0]):
        d[2 * p, 2 * q] = d[2 * p + 1, 2 * q + 1] = dxm[p, q]
Mu = np.zeros((len(dets),) * 2, dtype=complex)
idx = {dd: i for i, dd in enumerate(dets)}
for j, dd in enumerate(dets):
    b = np.zeros(1 << M, dtype=complex)
    b[dd] = 1.0
    col = np.zeros(1 << M, dtype=complex)
    for p in range(M):
        for q in range(M):
            if abs(d[p, q]) > 1e-14:
                col += d[p, q] * fullci.ad(p, fullci.aa(q, b))
    for d2, i in idx.items():
        if col[d2] != 0:
            Mu[i, j] = col[d2]
Mu = 0.5 * (Mu + Mu.conj().T)

K = sla.expm(-1j * EPS * Mu)
REGIMES = [("A adiabatic", esA, vecA),      # kicked ground state of H_A, propagated under H_A
           ("B quench", esB, vecB),          # ground state of H_A, propagated under H_B
           ("C quench rank2", esB, vecB)]    # rank-2 mixture of H_A states, propagated under H_B


def dm_of(tag):
    if tag.startswith("A"):
        v = K @ vecA[:, 0]
        v /= np.linalg.norm(v)
        return np.outer(v, v.conj())
    if tag.startswith("B"):
        v = vecA[:, 0]
        return np.outer(v, v.conj())
    v0, v1 = vecA[:, 0], vecA[:, 1]
    return 0.5 * (np.outer(v0, v0.conj()) + np.outer(v1, v1.conj()))


def cums_of(dm, es, vec, t, M):
    ph = np.exp(-1j * es.real * t)
    dmt = vec @ (ph[:, None] * (vec.conj().T @ dm @ vec) * ph.conj()[None, :]) @ vec.conj().T
    g = rdms(dmt, M, N, dets)
    l2 = az.cumulant(2, g[1], {("G", 1): g[0]})
    l3 = az.cumulant(3, g[2], {("G", 1): g[0], ("L", 2): l2})
    return g[0], g[1], l2, l3


print("Does the auxiliary basis need memory? Regimes where one-shot closures fail.")
print("compressed -> stretched H4 quench; rank-2 mixture for regime C.")
print()
for tag, es, vec in REGIMES:
    dm = dm_of(tag)
    print("=== %s ===" % tag)
    print("  %6s | %11s | %9s %9s %9s"
          % ("t", "|dg2 from 0|", "span err", "motion", "c drift"))
    c_ref, g2_0 = None, None
    for t in TIMES:
        g1, g2, l2, l3 = cums_of(dm, es, vec, t, M)
        g1b, g2b, l2b, l3b = cums_of(dm, es, vec, t + DT, M)
        if g2_0 is None:
            g2_0 = g2
        MOV = np.array([X.reshape(-1) for _, X in az.connected_basis(g1, l2, M)])
        MOVb = np.array([X.reshape(-1) for _, X in az.connected_basis(g1b, l2b, M)])
        c, *_ = np.linalg.lstsq(MOV.T, l3.reshape(-1), rcond=None)
        se = np.linalg.norm(MOV.T @ c - l3.reshape(-1)) / np.linalg.norm(l3)
        if c_ref is None:
            c_ref = c
        dl3 = (l3b - l3) / DT
        motion = np.linalg.norm(c @ ((MOVb - MOV) / DT)) / max(np.linalg.norm(dl3.reshape(-1)), 1e-30)
        print("  %6.2f | %11.5f | %9.4f %9.4f %9.4f"
              % (t, np.abs(g2 - g2_0).max(), se, motion,
                 np.linalg.norm(c - c_ref) / max(np.linalg.norm(c_ref), 1e-30)))
    print()
