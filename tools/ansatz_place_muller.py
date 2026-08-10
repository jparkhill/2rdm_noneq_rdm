#!/usr/bin/env python3
"""Where Mueller-2405 sits in the cumulant hierarchy: is it lambda_3 = 0 with rescaled
separable parts, or does it carry a genuine 3-cumulant?

APPROXIMATION UNDER EXAMINATION. None is propagated here; this is a static decomposition at the
exact FCI ground state. The object decomposed is the reconstructed 3-RDM
    gamma_3[p v x, q w y] = <a+_p a+_v a+_x a_y a_w a_q>,   Tr gamma_3 = N(N-1)(N-2),
produced by qdiss.gpu.muller2405_gamma3, which is Eq. 14 of Mueller, Luoma and Schaefer
arXiv:2405.05093 in the closed-system limit, optionally repaired by the Lackner Eq. 43/44/45
contraction correction. NOT the RDMFT Mueller exchange functional (that is muller_rdmft).

THE DECOMPOSITION. The MKNO generalized normal ordering (qdiss.symbolic.mkno) writes gamma_3 as a
sum over set partitions of the six indices into cumulant blocks. Grouping by which cumulants
appear, and writing "wedge" for the antisymmetrised product,
    gamma_3  =  S3[gamma_1]            the 3! = 6 terms built from gamma_1 alone,
              + S12[gamma_1, lambda_2] the terms with exactly one lambda_2 and one gamma_1,
              + lambda_3               the single fully-connected term.
Symbols: gamma_1[p,q] = <a+_p a_q>, Tr gamma_1 = N; lambda_2 = gamma_2 - gamma_1 wedge gamma_1 is
the 2-cumulant; lambda_3 is the 3-cumulant. "Tr_3" means contraction of the third index pair,
sum_x T[p v x, q w x].

THE TEST. Least squares for scalars (a, b, c) in
    gamma_3[recon]  ~=  a S3 + b S12 + c lambda_3[exact]
over all M^6 components. If Mueller is "lambda_3 = 0 with N-dependent rescaling" then c comes out
0 and (a, b) differ from (1, 1). Residual is reported relative to |gamma_3[exact]|, so a small
residual means the three structures span the reconstruction.

Ground truth: FCI in the N-electron sector, gamma_3 by qdiss.eom.fullci.rdm, whose rank-3 EOM
was validated against the same module at 8.7e-16.
"""
import sys
import numpy as np
import torch
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from molecular_2rdm import ORB_B, SEP3, build_H_sector, molecule, sector
from qdiss import gpu
from qdiss.eom import fullci
from qdiss.symbolic import mkno

UP3, LO3 = ("p", "v", "x"), ("q", "w", "y")
# The three MKNO groups, by which cumulants each set-partition term contains.
_ALL3 = [(s, ",".join("".join(u) + "".join(l) for n, r, u, l in f) + "->pvxqwy",
          [(n, r) for n, r, u, l in f])
         for s, f in mkno.gamma_to_cumulants(UP3, LO3)]
GRP_S3 = [t for t in _ALL3 if all(k == ("G", 1) for k in t[2])]
GRP_S12 = [t for t in _ALL3 if ("L", 2) in t[2] and ("L", 3) not in t[2]]
GRP_L3 = [t for t in _ALL3 if ("L", 3) in t[2]]

H4_TRAP = "H 0.00 0.00 0.0; H 0.80 0.00 0.0; H 1.05 0.72 0.0; H -0.15 0.72 0.0"
H3_BENT = "H 0.00 0.00 0.0; H 1.00 0.00 0.0; H 0.35 0.90 0.0"
SYS = ((H3_BENT, "sto-3g", 3, 0, 1, "H3 bent / sto-3g"),
       (H4_TRAP, "sto-3g", 4, 0, 0, "H4 asymmetric trapezoid / sto-3g"))

dev = gpu.device()


def build(group, tab):
    out = None
    for s, sub, keys in group:
        piece = s * torch.einsum(sub, *[tab[k] for k in keys])
        out = piece if out is None else out + piece
    return out


for atom, basis, N, chg, spin, label in SYS:
    mol, mf, M, h, V, enuc = molecule(atom, basis, chg, spin)
    dets = sector(M, N)
    Hs = build_H_sector(h, V, M, dets)
    es, vec = np.linalg.eigh(Hs)
    psi = np.zeros(1 << M, dtype=complex)
    psi[dets] = vec[:, 0]
    g1 = fullci.rdm(psi, M, 1)
    g2 = fullci.rdm(psi, M, 2)
    g3 = fullci.rdm(psi, M, 3)

    g1t, g2t, g3t = (gpu.up(x, dev) for x in (g1, g2, g3))
    l2t = g2t - gpu.wedge(g1t, g1t)
    tab = {("G", 1): g1t, ("L", 2): l2t}
    S3 = build(GRP_S3, tab)
    S12 = build(GRP_S12, tab)
    # lambda_3 exact, by subtracting the separable groups from the exact gamma_3.
    l3t = g3t - S3 - S12
    # Self-check: the MKNO L3 group must be the identity on lambda_3.
    chk = build(GRP_L3, {("L", 3): l3t}) - l3t

    print("=== %s   M=%d N=%d   E_FCI %.9f ===" % (label, M, N, es[0].real + enuc))
    print("    MKNO group sizes: S3 %d terms, S12 %d terms, L3 %d terms"
          % (len(GRP_S3), len(GRP_S12), len(GRP_L3)))
    print("    identity check |L3group(lam3) - lam3| = %.2e" % float(chk.abs().max()))
    nrm = float(g3t.abs().max())
    print("    |gamma_3|_max %.4f   |S3|_max %.4f   |S12|_max %.4f   |lambda_3|_max %.4f"
          % (nrm, float(S3.abs().max()), float(S12.abs().max()), float(l3t.abs().max())))

    basis_t = [S3, S12, l3t]
    A = np.array([[float(torch.vdot(x.reshape(-1), y.reshape(-1)).real) for y in basis_t]
                  for x in basis_t])
    print()
    print("    least squares  gamma_3[recon] ~ a*S3 + b*S12 + c*lambda_3[exact]")
    print("    %-16s %9s %9s %9s | %10s %10s"
          % ("reconstruction", "a", "b", "c", "resid/|g3|", "|g3rec-g3|"))
    for nm, T in (("exact", g3t),
                  ("muller2405", gpu.muller2405_gamma3(g1t, g2t, N, M, True, ORB_B)),
                  ("muller2405_raw", gpu.muller2405_gamma3(g1t, g2t, N, M, False, ORB_B)),
                  ("lam3_0", S3 + S12)):
        b = np.array([float(torch.vdot(x.reshape(-1), T.reshape(-1)).real) for x in basis_t])
        cf = np.linalg.solve(A, b)
        fit = cf[0] * S3 + cf[1] * S12 + cf[2] * l3t
        print("    %-16s %9.5f %9.5f %9.5f | %10.2e %10.2e"
              % (nm, cf[0], cf[1], cf[2],
                 float((T - fit).abs().max()) / nrm, float((T - g3t).abs().max())))
    print()
