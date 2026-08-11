#!/usr/bin/env python3
r"""Optimise the level-1 ladder on its Stiefel manifold: expressive power and variational bound.

WHY NOT AN SVD TRUNCATION. tools/ladder_level1.py bounded the achievable error by truncating the
exact Mat_2 to rank chi'. That was worthless: truncation drops positive eigenvalues, so it violates
tr gamma_2 = N(N-1), and at chi' = 1 it returned E = -4.21, which is the one-body energy with
gamma_2 annihilated. The ladder family satisfies the contraction identity EXACTLY, so the truncation
is not in the family and bounds nothing. The family must be optimised over directly.

THE MANIFOLD. Stacking the M ladders vertically into frak_A of size (M chi') x r_1, the constraint
sum_x frak_a_x^dag frak_a_x = (N-1) I_{r_1} reads frak_A^dag frak_A = (N-1) I, i.e.

    frak_A = sqrt(N-1) * Theta ,      Theta^dag Theta = I ,      Theta in St(M chi', r_1)

So the closure is an unconstrained point on a Stiefel manifold. Parameterise by a free complex
matrix and retract with QR; every iterate satisfies the constraint exactly, hence every iterate
gives a gamma_2 that is D-positive, contraction-consistent with gamma_1, and obeys the CAR hole
condition. No penalty, no projection, no repair.

    Mat_2[(p,v),(q,w)] = 2 Rc(:,p)^dag G_vw Rc(:,q) ,    G_vw = frak_a_v^dag frak_a_w
    gamma_2 on sorted pairs = Mat_2 / 2
    E = sum_pq h gamma_1 + sum_{a,b over sorted pairs} V[a,b] Mat_2[a,b] / 2 + E_nuc

(the two-body prefactor 0.25 of qdiss.eom.hamiltonian becomes 1 on sorted pairs because the four
sign-equivalent index orders each contribute equally; this is checked against the exact energy.)

TWO OBJECTIVES, which answer different questions and must not be confused:

  ORACLE       minimise |Mat_2(Theta) - Mat_2^exact|_F with gamma_1 fixed exact. Measures the
               EXPRESSIVE POWER of the family at a given chi'. Zero means the family contains the
               truth.
  VARIATIONAL  minimise E with gamma_1 fixed exact. The family is contained in the D-positive,
               contraction-consistent, hole-consistent set, which is LARGER than the
               N-representable set, so this gives a LOWER bound on E_FCI. How far below measures
               how much the imposed constraints fail to pin the energy.

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
STEPS, LR = 4000, 0.02


def run(atom, basis, N, spin, label):
    mol, mf, M, h, V, enuc = molecule(atom, basis, 0, spin)
    dets = sector(M, N)
    es, vec = np.linalg.eigh(build_H_sector(h, V, M, dets))
    psi = np.zeros(1 << M, dtype=complex)
    psi[dets] = vec[:, 0]
    g1, g2 = fullci.rdm(psi, M, 1), fullci.rdm(psi, M, 2)
    e_fci, e_hf = float(fci.FCI(mf).kernel()[0]), float(mf.e_tot)
    PAIR = az.tuples(M, 2)
    NP = len(PAIR)
    Vp = np.array([[V[p, v, q, w] for (q, w) in PAIR] for (p, v) in PAIR])
    Mat2ex = az.tensor_to_coord(g2, M, 2) * 2.0
    e1 = float(np.real(np.einsum("pq,pq->", h, g1))) + enuc
    # convention check: the sorted-pair contraction must reproduce the exact energy
    e_chk = e1 + float(np.real(np.einsum("ab,ab->", Vp, Mat2ex)) / 2.0)
    U1 = np.array([np.array([fullci.aa(p, psi)[d] for d in sector(M, N - 1)])
                   for p in range(M)]).T
    Q, sv, Vh = np.linalg.svd(U1, full_matrices=False)
    kp = sv > 1e-10 * sv[0]
    r1 = int(kp.sum())
    Rc = np.diag(sv[kp]) @ Vh[kp, :]
    print("=== %s   M=%d N=%d   r_1=%d   E_HF %.6f  E_FCI %.6f  E_corr %.6f ==="
          % (label, M, N, r1, e_hf, e_fci, e_fci - e_hf))
    print("    convention check: |E(sorted-pair contraction) - E_FCI| = %.2e"
          % abs(e_chk - e_fci))

    Rt = torch.tensor(Rc)
    Vt = torch.tensor(Vp)
    Mex = torch.tensor(Mat2ex)
    idx = torch.tensor([[p, v] for (p, v) in PAIR])

    def build(Theta, chi):
        A = np.sqrt(N - 1.0) * Theta.reshape(M, chi, r1)
        G = torch.einsum("xar,yas->xyrs", A.conj(), A)
        # Mat_2[(p,v),(q,w)] = 2 Rc(:,p)^dag G_vw Rc(:,q)
        left = Rt[:, idx[:, 0]]                       # r1 x NP, column a -> p of pair a
        Gsel = G[idx[:, 1], :][:, idx[:, 1]]          # NP x NP x r1 x r1 by v,w
        return 2.0 * torch.einsum("ra,abrs,sb->ab", left.conj(), Gsel, left)

    print("    %6s | %13s %13s | %13s %13s %9s"
          % ("chi'", "oracle |dM|", "rel", "var E", "E-E_FCI", "%Ecorr"))
    for chi in (2, 3, 4, 6, 8, 12, 16, 28):
        if M * chi < r1:
            continue
        out = {}
        for mode in ("oracle", "var"):
            g = torch.Generator().manual_seed(11)
            Z = (torch.randn(M * chi, r1, generator=g)
                 + 1j * torch.randn(M * chi, r1, generator=g)).requires_grad_(True)
            opt = torch.optim.Adam([Z], lr=LR)
            best = None
            for it in range(STEPS):
                opt.zero_grad()
                Theta, _ = torch.linalg.qr(Z)
                Mat2 = build(Theta, chi)
                if mode == "oracle":
                    loss = torch.linalg.matrix_norm(Mat2 - Mex) ** 2
                else:
                    loss = torch.real(torch.einsum("ab,ab->", Vt, Mat2)) / 2.0
                loss.backward()
                opt.step()
                v = float(loss.detach())
                best = v if best is None else min(best, v)
            out[mode] = best
        eo = np.sqrt(out["oracle"])
        ev = e1 + out["var"]
        print("    %6d | %13.3e %13.4f | %13.6f %13.6f %8.1f%%"
              % (chi, eo, eo / float(np.linalg.norm(Mat2ex)), ev, ev - e_fci,
                 100.0 * (ev - e_hf) / (e_fci - e_hf)))
    print()


for a, b, n, s, l in (CASES if __name__ == "__main__" else ()):
    run(a, b, n, s, l)
