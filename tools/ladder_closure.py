#!/usr/bin/env python3
r"""Level-1 ladder closure WITH the contraction identity imposed as a linear condition on G.

WHY THIS SUPERSEDES tools/ladder_optimise.py. That run imposed only the number condition
sum_x frak_a_x^dag frak_a_x = (N-1) I, which commit 34183c7 showed does NOT imply the contraction
identity for a MODELLED ladder: random Stiefel ladders gave tr Mat_2 ~ 16 against N(N-1) = 12, and
the variational objective ran away to E = -4.87 against E_FCI = -1.943. The identity must be imposed.

PARAMETERISATION. frak_A is a free complex chi' x (M r_1) matrix whose column blocks are the ladders
frak_a_x of size chi' x r_1. Then

    G_xy = frak_a_x^dag frak_a_y,    assembled G = frak_A^dag frak_A  is (M r_1) x (M r_1)

so G >= 0 and rank G <= chi' are STRUCTURAL -- hence gamma_2 >= 0 is structural, and at N = k the
number condition would force every ladder to vanish. The remaining conditions are LINEAR in G and
are imposed by a ramped quadratic penalty, with residuals REPORTED so the reader can judge:

    (C1)  sum_x G_xx = (N-1) I_{r_1}                     r_1^2 conditions
    (C2)  Tr_2 gamma_2 = (N-1) gamma_1                   M^2 conditions

RECONSTRUCTION, verified to 1.4e-15 with the exact ladder in tools/ladder_level1.py:
    Mat_2[(p,v),(q,w)] = 2 Rc(:,p)^dag G_vw Rc(:,q)      over sorted pairs
    gamma_2[p,v,q,w]   = sign * Mat_2[sort(p,v), sort(q,w)] / 2,  zero when p = v or q = w
    E = sum_pq h gamma_1 + sum_{a,b sorted} V[a,b] Mat_2[a,b] / 2 + E_nuc

TWO OBJECTIVES, distinct questions, not to be confused:
  ORACLE       min |Mat_2 - Mat_2^exact|_F at fixed exact gamma_1. The EXPRESSIVE POWER of the family
               at a given chi'. Zero means the family contains the truth.
  VARIATIONAL  min E at fixed exact gamma_1. The family lies inside the D-positive,
               contraction-consistent set, which is LARGER than the N-representable set, so this is a
               LOWER bound on E_FCI. The gap measures how far the imposed conditions fall short of
               pinning the energy -- the same quantity a DQG semidefinite calculation reports.

Ground truth qdiss.eom.fullci and the PySCF FCI energy.
"""
import sys
import numpy as np
import torch
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from molecular_2rdm import build_H_sector, molecule, sector
from pyscf import fci
from qdiss import ansatz as az
from qdiss.eom import fullci

torch.set_default_dtype(torch.float64)
CASES = (("H 0.00 0.00 0.0; H 0.80 0.00 0.0; H 1.05 0.72 0.0; H -0.15 0.72 0.0",
          "sto-3g", 4, 0, "H4 trapezoid"),
         ("Li 0 0 0; H 0 0 1.6", "sto-3g", 4, 0, "LiH"))
STEPS, LR = 3000, 0.03


def run(atom, basis, N, spin, label):
    mol, mf, M, h, V, enuc = molecule(atom, basis, 0, spin)
    dets = sector(M, N)
    es, vec = np.linalg.eigh(build_H_sector(h, V, M, dets))
    psi = np.zeros(1 << M, dtype=complex)
    psi[dets] = vec[:, 0]
    g1, g2 = fullci.rdm(psi, M, 1), fullci.rdm(psi, M, 2)
    e_fci, e_hf = float(fci.FCI(mf).kernel()[0]), float(mf.e_tot)
    PAIR = az.tuples(M, 2)
    pos = {t: i for i, t in enumerate(PAIR)}
    Vp = np.array([[V[p, v, q, w] for (q, w) in PAIR] for (p, v) in PAIR])
    Mat2ex = az.tensor_to_coord(g2, M, 2) * 2.0
    e1 = float(np.real(np.einsum("pq,pq->", h, g1))) + enuc

    AI = np.zeros((M, M), dtype=int)
    SA = np.zeros((M, M))
    for p in range(M):
        for v in range(M):
            if p == v:
                continue
            AI[p, v] = pos[(min(p, v), max(p, v))]
            SA[p, v] = 1.0 if p < v else -1.0

    U1 = np.array([fullci.aa(p, psi)[sector(M, N - 1)] for p in range(M)]).T
    Q, sv, Vh = np.linalg.svd(U1, full_matrices=False)
    kp = sv > 1e-10 * sv[0]
    r1 = int(kp.sum())
    Rc = np.diag(sv[kp]) @ Vh[kp, :]

    Rt, Vt, Mex = torch.tensor(Rc), torch.tensor(Vp), torch.tensor(Mat2ex)
    g1t = torch.tensor(g1)
    pidx = torch.tensor([p for (p, v) in PAIR])
    vidx = torch.tensor([v for (p, v) in PAIR])
    AIt, SAt = torch.tensor(AI), torch.tensor(SA)
    Ir = torch.eye(r1, dtype=torch.complex128)

    print("=== %s  M=%d N=%d  r_1=%d | E_HF %.6f  E_FCI %.6f  E_corr %.6f ==="
          % (label, M, N, r1, e_hf, e_fci, e_fci - e_hf))

    def mat2_of(A, chi):
        a = A.reshape(chi, M, r1)
        G = torch.einsum("cxr,cys->xyrs", a.conj(), a)
        left = Rt[:, pidx]
        Gsel = G[vidx, :][:, vidx]
        return 2.0 * torch.einsum("ra,abrs,sb->ab", left.conj(), Gsel, left), G

    def resid(Mat2, G):
        c1 = torch.einsum("xxrs->rs", G) - (N - 1.0) * Ir
        g2t = 0.5 * SAt.unsqueeze(-1).unsqueeze(-1) * SAt * Mat2[AIt][:, :, AIt]
        return c1, torch.einsum("pxqx->pq", g2t) - (N - 1.0) * g1t

    print("    %5s | %11s %8s | %13s %11s | %9s %9s"
          % ("chi'", "oracle", "rel", "var E", "E-E_FCI", "|C1|", "|C2|"))
    for chi in (2, 4, 6, 8, 12, 16, 28):
        row = {}
        for mode in ("oracle", "var"):
            g = torch.Generator().manual_seed(11)
            A = ((torch.randn(chi, M * r1, generator=g)
                  + 1j * torch.randn(chi, M * r1, generator=g)) * 0.5).requires_grad_(True)
            opt = torch.optim.Adam([A], lr=LR)
            for it in range(STEPS):
                opt.zero_grad()
                Mat2, G = mat2_of(A, chi)
                c1, c2 = resid(Mat2, G)
                pen = 10.0 ** (1.0 + 4.0 * it / STEPS)
                base = (torch.linalg.matrix_norm(Mat2 - Mex) ** 2 if mode == "oracle"
                        else torch.real(torch.einsum("ab,ab->", Vt, Mat2)) / 2.0)
                (base + pen * (torch.linalg.matrix_norm(c1) ** 2
                               + torch.linalg.matrix_norm(c2) ** 2)).backward()
                opt.step()
            with torch.no_grad():
                Mat2, G = mat2_of(A, chi)
                c1, c2 = resid(Mat2, G)
                row[mode] = (float(torch.linalg.matrix_norm(Mat2 - Mex)),
                             float(torch.real(torch.einsum("ab,ab->", Vt, Mat2)) / 2.0),
                             float(torch.linalg.matrix_norm(c1)),
                             float(torch.linalg.matrix_norm(c2)))
        ev = e1 + row["var"][1]
        print("    %5d | %11.3e %8.4f | %13.6f %11.6f | %9.1e %9.1e"
              % (chi, row["oracle"][0], row["oracle"][0] / float(np.linalg.norm(Mat2ex)),
                 ev, ev - e_fci, row["var"][2], row["var"][3]))
    print()


for a, b, n, s, l in (CASES if __name__ == "__main__" else ()):
    run(a, b, n, s, l)
