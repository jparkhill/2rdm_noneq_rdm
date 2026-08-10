#!/usr/bin/env python3
"""Self-validation of qdiss.ansatz before any of it is used.

Four checks, at ranks 2 AND 3 from identical code, on exact FCI states so every quantity has a
known value:

  1. round trip  coord_to_tensor(tensor_to_coord(L)) == L for the exact cumulants, confirming the
     sorted-tuple representation loses nothing.
  2. trace target  Tr_2 lambda_2 == -(gamma_1 - gamma_1^2) and Tr_3 lambda_3 == anti(E, lambda_2),
     the order-general single-trace identity, recomputed here through the new anti().
  3. min_norm_trace  reproduces the target trace, is Hermitian and antisymmetric, and has norm
     no larger than the exact cumulant with the same trace (a necessary condition for being the
     minimum-norm solution) -- plus the sharp test, that the residual lambda_exact - lambda_minnorm
     is orthogonal to the min-norm solution, which holds iff min_norm sits in the row space.
  4. agreement with the existing closed forms already in the codebase: closures.trace_correct's
     Corr map at k = 2 and qdiss.gpu.trace_inverse at k = 3. Those are particular solutions, not
     necessarily minimum-norm; the check is that they have the SAME trace, and the norm gap is
     reported.

Symbols as in qdiss.ansatz. Ground truth is qdiss.eom.fullci in the N-electron sector.
"""
import sys
import numpy as np
import torch
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from molecular_2rdm import ORB_B, build_H_sector, molecule, sector
from qdiss import ansatz as az, closures, gpu
from qdiss.eom import fullci

H4_TRAP = "H 0.00 0.00 0.0; H 0.80 0.00 0.0; H 1.05 0.72 0.0; H -0.15 0.72 0.0"
H3_BENT = "H 0.00 0.00 0.0; H 1.00 0.00 0.0; H 0.35 0.90 0.0"
SYS = ((H3_BENT, "sto-3g", 3, 0, 1, "H3 bent"), (H4_TRAP, "sto-3g", 4, 0, 0, "H4 trapezoid"))
dev = gpu.device()


def fro(x):
    return float(np.sqrt(np.abs(np.vdot(x, x)).real))


for atom, basis, N, chg, spin, label in SYS:
    mol, mf, M, h, V, enuc = molecule(atom, basis, chg, spin)
    dets = sector(M, N)
    es, vec = np.linalg.eigh(build_H_sector(h, V, M, dets))
    psi = np.zeros(1 << M, dtype=complex)
    psi[dets] = vec[:, 0]
    g1 = fullci.rdm(psi, M, 1)
    g2 = fullci.rdm(psi, M, 2)
    g3 = fullci.rdm(psi, M, 3)
    l2 = az.cumulant(2, g2, {("G", 1): g1})
    l3 = az.cumulant(3, g3, {("G", 1): g1, ("L", 2): l2})
    print("=== %s   M=%d N=%d   E %.9f ===" % (label, M, N, es[0].real + enuc))
    print("    %-4s %-14s %11s %11s %11s %11s %11s"
          % ("k", "object", "roundtrip", "trace res", "|lam_k|", "|minnorm|", "orthog"))
    for k, lam, low in ((2, l2, g1), (3, l3, l2)):
        rt = float(np.abs(az.coord_to_tensor(az.tensor_to_coord(lam, M, k), M, k) - lam).max())
        S = az.trace_target(k, g1, low)
        tr = az.trace_residual(lam, S, k)
        mn = az.min_norm_trace(S, M, k)
        trm = az.trace_residual(mn, S, k)
        # orthogonality of the residual to the min-norm solution
        orth = abs(complex(np.vdot(mn, lam - mn))) / max(fro(mn) * fro(lam - mn), 1e-30)
        print("    %-4d %-14s %11.2e %11.2e %11.5f %11.5f %11.2e"
              % (k, "lambda_%d exact" % k, rt, tr, fro(lam), fro(mn), orth))
        print("    %-4s %-14s %11s %11.2e" % ("", "min_norm_trace", "", trm))

        # structural checks on the min-norm solution
        ax_u = np.abs(mn + np.swapaxes(mn, 0, 1)).max()
        herm = np.abs(mn - np.conj(np.moveaxis(mn, range(2 * k), list(range(k, 2 * k)) + list(range(k))))).max()
        print("    %-4s %-14s antisym %8.2e   herm %8.2e" % ("", "", ax_u, herm))

        # comparison with the existing closed forms
        if k == 2:
            I = np.eye(M, dtype=complex)
            Sc = np.einsum("pxqx->pq", np.zeros((M,) * 4, dtype=complex))  # placeholder
            other = (1.0 / (M - 2)) * closures.sym_wedge(S, I) \
                - (np.trace(S) / ((M - 2) * (M - 1))) * closures.wedge(I, I)
            nm = "closures Corr(S)"
        else:
            other = gpu.down(gpu.trace_inverse(gpu.up(S, dev), M, ORB_B))
            nm = "gpu.trace_inverse"
        print("    %-4s %-14s %11s %11.2e %11s %11.5f  gap %+.3f%%"
              % ("", nm, "", az.trace_residual(other, S, k), "", fro(other),
                 100.0 * (fro(other) - fro(mn)) / fro(mn)))
    print()
