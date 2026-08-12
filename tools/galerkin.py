#!/usr/bin/env python3
r"""The Galerkin machinery for the coefficient flow, validated, with its truncation error budget.

THE SCHEME. lambda_3 = sum_a c_a S_a with S_a the connected structures. Writing
lambda_3 = sum_a c_a S_a + r with r orthogonal to the span,

    lambda_3-dot = sum_a c_a-dot S_a + sum_a c_a S_a-dot + r-dot

and projecting onto S_b gives the GALERKIN coefficient flow

    Gram c-dot = <S, lambda_3-dot>  -  <S, sum_a c_a S_a-dot>          Gram_ab = <S_a, S_b>

The exact coefficient rate, meanwhile, is c-dot_exact = d/dt of the fitted coefficients. Subtracting,

    c-dot_Galerkin - c-dot_exact = Gram^{-1} <S, r-dot>

so the discrepancy is exactly the rate at which the OUT-OF-SPAN residual feeds back into the span.
That is the truncation error of the scheme, and it is the number that decides whether a Galerkin
projection on this basis is worth integrating. It is measured here BEFORE any integrator is written.

VALIDATIONS FIRST, as the project rules require:
  V1  the MKNO rank-3 term list reproduces fullci.rdm_dot for gamma_3
  V2  the analytic cumulant rate, lambda_3-dot = gamma_3-dot - separable_dot(3, ...), matches a
      finite difference of lambda_3 along the exact trajectory
  V3  the analytic structure rates from connected_basis_dot match finite differences of the
      structures themselves
only then
  E1  |c-dot_Galerkin - c-dot_exact| / |c-dot_exact|, the leakage, in three dynamical regimes.

Symbols as in qdiss.ansatz. The basis is orthonormalised and frame-aligned (see
tools/basis_motion_quench.py) so the coefficients are unique; the raw generators span a
4-to-6-dimensional space depending on threshold, so a raw lstsq coefficient is not well defined.

Ground truth qdiss.eom.fullci.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from ansatz_hierarchy import rdms
from molecular_2rdm import build_H_sector, molecule, sector
from qdiss import ansatz as az
from qdiss.eom import fullci

COMPRESSED = "H 0.00 0.00 0.0; H 0.60 0.00 0.0; H 0.79 0.54 0.0; H -0.11 0.54 0.0"
STRETCHED = "H 0.00 0.00 0.0; H 1.20 0.00 0.0; H 1.58 1.08 0.0; H -0.22 1.08 0.0"
N, DT = 4, 1e-5


def setup(atom):
    mol, mf, M, h, V, enuc = molecule(atom, "sto-3g", 0, 0)
    dets = sector(M, N)
    es, vec = np.linalg.eigh(build_H_sector(h, V, M, dets))
    return mol, mf, M, h, V, dets, es, vec


molA, mfA, M, hA, VA, dets, esA, vecA = setup(COMPRESSED)
molB, mfB, _, hB, VB, _, esB, vecB = setup(STRETCHED)
HB = fullci.build_H(hB, VB, M)
EYE = np.eye(M, dtype=complex)


def cums(dm, es, vec, t):
    ph = np.exp(-1j * es.real * t)
    dmt = vec @ (ph[:, None] * (vec.conj().T @ dm @ vec) * ph.conj()[None, :]) @ vec.conj().T
    g = rdms(dmt, M, N, dets)
    l2 = az.cumulant(2, g[1], {("G", 1): g[0]})
    l3 = az.cumulant(3, g[2], {("G", 1): g[0], ("L", 2): l2})
    return g[0], g[1], g[2], l2, l3


def frame(g1, l2, ref=None, kfix=None):
    B = np.array([X.reshape(-1) for _, X in az.connected_basis(g1, l2, M)])
    _, sv, Vt = np.linalg.svd(B, full_matrices=False)
    k = int((sv > 1e-9 * sv[0]).sum()) if kfix is None else kfix
    F = Vt[:k]
    if ref is not None and ref.shape[0] == k:
        W, _, Zt = np.linalg.svd(F @ ref.conj().T)
        F = (W @ Zt).conj().T @ F
    return F, B


# ---------------------------------------------------------------- V1, V2, V3 on a pure-state trajectory
psi = np.zeros(1 << M, dtype=complex)
psi[dets] = vecA[:, 0]
dm0 = np.outer(vecA[:, 0], vecA[:, 0].conj())
g1, g2, g3, l2, l3 = cums(dm0, esB, vecB, 0.0)
g1b, g2b, g3b, l2b, l3b = cums(dm0, esB, vecB, DT)

d1 = fullci.rdm_dot(psi, HB, M, 1)
d2 = fullci.rdm_dot(psi, HB, M, 2)
d3 = fullci.rdm_dot(psi, HB, M, 3)
dl2 = d2 - (az.wedge(d1, g1) + az.wedge(g1, d1))
print("VALIDATIONS (compressed ground state propagated under the stretched Hamiltonian)")
print("  V1 |gamma_3-dot fullci vs finite difference|     = %.3e   (|d3| = %.4f)"
      % (np.abs(d3 - (g3b - g3) / DT).max(), np.abs(d3).max()))
dl3 = d3 - az.separable_dot(3, {("G", 1): g1, ("L", 2): l2}, {("G", 1): d1, ("L", 2): dl2})
print("  V2 |lambda_3-dot analytic vs finite difference|  = %.3e   (|dl3| = %.4f)"
      % (np.abs(dl3 - (l3b - l3) / DT).max(), np.abs(dl3).max()))
Sd = az.connected_basis_dot(g1, l2, d1, dl2, M)
S0 = [X for _, X in az.connected_basis(g1, l2, M)]
S1 = [X for _, X in az.connected_basis(g1b, l2b, M)]
worst = max(np.abs(Sd[i] - (S1[i] - S0[i]) / DT).max() / max(np.abs(Sd[i]).max(), 1e-30)
            for i in range(len(Sd)))
print("  V3 worst relative error of the %d structure rates = %.3e" % (len(Sd), worst))
print()

# ---------------------------------------------------------------- E1 the leakage, three regimes
REG = [("A adiabatic", dm0, esA, vecA),
       ("B quench", dm0, esB, vecB),
       ("C quench rank2",
        0.5 * (np.outer(vecA[:, 0], vecA[:, 0].conj()) + np.outer(vecA[:, 1], vecA[:, 1].conj())),
        esB, vecB)]
print("E1  leakage = |c-dot_Galerkin - c-dot_exact| / |c-dot_exact| = |Gram^-1 <S, r-dot>| / |c-dot|")
print("    %-16s %6s | %10s %10s %10s %10s"
      % ("regime", "t", "span err", "|cdot ex|", "|cdot Gal|", "leakage"))
for tag, dm, es, vec in REG:
    pf, pk = None, None
    for t in (0.0, 0.25, 0.5, 1.0):
        g1, g2, g3, l2, l3 = cums(dm, es, vec, t)
        g1b, g2b, g3b, l2b, l3b = cums(dm, es, vec, t + DT)
        F, _ = frame(g1, l2, pf, pk)
        pf, pk = F, F.shape[0]
        Fb, _ = frame(g1b, l2b, F, F.shape[0])
        c = F.conj() @ l3.reshape(-1)
        cb = Fb.conj() @ l3b.reshape(-1)
        se = np.linalg.norm(F.T @ c - l3.reshape(-1)) / np.linalg.norm(l3)
        cdot_ex = (cb - c) / DT
        # analytic rates at this state
        d1 = np.einsum("pxqx->pq", (g2b - g2) / DT) / (N - 1)   # gamma_1-dot by contraction
        d1 = (g1b - g1) / DT
        dl2 = (l2b - l2) / DT
        dl3 = (l3b - l3) / DT
        Fdot = (Fb - F) / DT
        # Galerkin: Gram is identity in the orthonormal frame, so c-dot = <F, dl3> - <F, Fdot c>
        cdot_gal = F.conj() @ dl3.reshape(-1) - (F.conj() @ (Fdot.T @ c))
        leak = np.linalg.norm(cdot_gal - cdot_ex) / max(np.linalg.norm(cdot_ex), 1e-30)
        print("    %-16s %6.2f | %10.4f %10.4f %10.4f %10.4f"
              % (tag, t, se, np.linalg.norm(cdot_ex), np.linalg.norm(cdot_gal), leak))
    print()
