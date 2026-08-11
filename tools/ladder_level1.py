#!/usr/bin/env python3
"""Level-1 ladder closure: validate the construction, then bound what it can achieve vs chi'.

THE SCHEME AT LEVEL k = 1. The propagated object is Mat_1 = gamma_1 = R_c^dag R_c with R_c of size
r_1 x M, r_1 = rank gamma_1. The closure is one ladder frak_a_x : C^{r_1} -> C^{chi'}, x = 1..M,
subject to sum_x frak_a_x^dag frak_a_x = (N-1) I_{r_1}, and it supplies gamma_2 through

    Mat_2[(p,v),(q,w)] = 2 R_c(:,p)^dag G_vw R_c(:,q),     G_vw = frak_a_v^dag frak_a_w
    gamma_2 = coord_to_tensor(Mat_2 / 2)

Then E = sum_pq h[p,q] gamma_1[p,q] + (1/4) sum V[p,v,q,w] gamma_2[p,v,q,w] + E_nuc, with the
Hamiltonian convention of qdiss.eom.hamiltonian (TWO_BODY_PREFACTOR = 0.25).

THE CONSTRAINT SET IS A STIEFEL MANIFOLD. Stacking the ladders vertically into a single
(M chi') x r_1 matrix frak_A, the condition reads frak_A^dag frak_A = (N-1) I, i.e.
frak_A = sqrt(N-1) * Theta with Theta^dag Theta = I. So the closure is an unconstrained point on
St(M chi', r_1) and can be optimised by ordinary retraction methods -- no penalty, no projection.

WHAT IS COMPUTED HERE, in the order the project rules demand (validate conventions first):

  V1. energy from the EXACT gamma_1, gamma_2 must reproduce the FCI energy. Checks the
      contraction convention and the prefactor before anything is built on them.
  V2. the EXACT level-1 ladder, built from the FCI state, must give back the exact gamma_2 and
      hence E_FCI. Checks the G-transfer formula and the sorted-tuple bookkeeping.
  V3. the MEAN-FIELD ladder. For a determinant |O>, V^{(k)} = span{|O\I> : I subset O}, the ladder
      is the free-fermion annihilator on hole configurations, and sum_x frak_a_x^dag frak_a_x
      counts the remaining occupied orbitals = (N-k) exactly. Fed the HF gamma_1 it must return the
      HF energy, i.e. gamma_2 = gamma_1 wedge gamma_1. This is the anchor the closure dresses.

  B1. A BOUND on what any chi' can achieve. Since Mat_2 = R_c^dag G R_c with rank G <= chi', the
      family's Mat_2 has rank <= chi'. The best rank-chi' approximation to the exact Mat_2 is its
      SVD truncation, so the energy error of that truncation is a LOWER BOUND on the error of any
      ladder at that chi'. Reported against E_FCI and against the HF and MP2-scale correlation
      energy for context. Caveat stated rather than hidden: the truncation is not itself of
      G-transfer form and need not be contraction-consistent, so it bounds but does not attain.

Ground truth qdiss.eom.fullci and the PySCF FCI energy.
"""
import sys
from math import factorial
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from molecular_2rdm import build_H_sector, molecule, sector
from pyscf import fci
from qdiss import ansatz as az
from qdiss.eom import fullci

CASES = (("H 0.00 0.00 0.0; H 0.80 0.00 0.0; H 1.05 0.72 0.0; H -0.15 0.72 0.0",
          "sto-3g", 4, 0, "H4 trapezoid"),
         ("Li 0 0 0; H 0 0 1.6", "sto-3g", 4, 0, "LiH"))


def annihilate(p, tup):
    out = p
    for t in tup:
        out = fullci.aa(t, out)
    return out


def energy(g1, g2, h, V, enuc):
    return float(np.real(np.einsum("pq,pq->", h, g1)
                         + 0.25 * np.einsum("pvqw,pvqw->", V, g2)) + enuc)


def gamma2_from_ladder(Rc, G, M):
    """Mat_2[(p,v),(q,w)] = 2 Rc(:,p)^dag G_vw Rc(:,q) over sorted pairs, then antisymmetrise."""
    PAIR = az.tuples(M, 2)
    Y = np.zeros((len(PAIR),) * 2, dtype=complex)
    for a, (p, v) in enumerate(PAIR):
        for b, (q, w) in enumerate(PAIR):
            Y[a, b] = 2.0 * (Rc[:, p].conj() @ (G[v][w] @ Rc[:, q]))
    return az.coord_to_tensor(Y / 2.0, M, 2), Y


for atom, basis, N, spin, label in CASES:
    mol, mf, M, h, V, enuc = molecule(atom, basis, 0, spin)
    dets = sector(M, N)
    es, vec = np.linalg.eigh(build_H_sector(h, V, M, dets))
    psi = np.zeros(1 << M, dtype=complex)
    psi[dets] = vec[:, 0]
    g1, g2 = fullci.rdm(psi, M, 1), fullci.rdm(psi, M, 2)
    e_fci_ref = fci.FCI(mf).kernel()[0]
    e_hf = float(mf.e_tot)
    print("=== %s / %s   M=%d N=%d ===" % (label, basis, M, N))
    print("    E_HF %.9f   E_FCI(pyscf) %.9f   E_corr %.6f"
          % (e_hf, e_fci_ref, e_fci_ref - e_hf))

    # V1
    e_exact = energy(g1, g2, h, V, enuc)
    print("    V1  E[gamma_1, gamma_2 exact]      = %.9f   dev from FCI %.2e"
          % (e_exact, abs(e_exact - e_fci_ref)))

    # the exact level-1 tower
    S1, S2 = sector(M, N - 1), sector(M, N - 2)
    U1 = np.array([annihilate(psi, (p,))[S1] for p in range(M)]).T
    Q, sv, Vh = np.linalg.svd(U1, full_matrices=False)
    kp = sv > 1e-10 * sv[0]
    r1 = int(kp.sum())
    P1 = Q[:, kp]
    Rc = np.diag(sv[kp]) @ Vh[kp, :]                # gamma_1 = Rc^dag Rc  (Mat_1 = 1! gamma_1)
    U2 = np.array([annihilate(psi, I)[S2] for I in az.tuples(M, 2)]).T
    Q2, sv2, _ = np.linalg.svd(U2, full_matrices=False)
    kp2 = sv2 > 1e-10 * sv2[0]
    P2 = Q2[:, kp2]
    chi_ex = int(kp2.sum())
    lad = np.zeros((M, chi_ex, r1), dtype=complex)
    for x in range(M):
        for al in range(r1):
            b = np.zeros(1 << M, dtype=complex)
            b[S1] = P1[:, al]
            lad[x, :, al] = P2.conj().T @ fullci.aa(x, b)[S2]
    Gex = [[lad[x].conj().T @ lad[y] for y in range(M)] for x in range(M)]
    print("    r_1 = rank gamma_1 = %d of %d,  chi'(exact) = rank Mat_2 = %d"
          % (r1, M, chi_ex))
    print("        |gamma_1 - Rc^dag Rc| = %.2e   |sum_x a^dag a - (N-1)I| = %.2e"
          % (np.abs(g1 - Rc.conj().T @ Rc).max(),
             np.abs(sum(Gex[x][x] for x in range(M)) - (N - 1) * np.eye(r1)).max()))

    # V2
    g2L, Mat2L = gamma2_from_ladder(Rc, Gex, M)
    print("    V2  exact ladder: |gamma_2 - exact| = %.2e   E = %.9f  dev %.2e"
          % (np.abs(g2L - g2).max(), energy(g1, g2L, h, V, enuc),
             abs(energy(g1, g2L, h, V, enuc) - e_fci_ref)))

    # V3  mean-field ladder from the HF determinant
    occ = [i for i in range(M) if i < N]           # HF occupies the N lowest spin orbitals
    det = 0
    for i in occ:
        det |= (1 << i)
    phi = np.zeros(1 << M, dtype=complex)
    phi[det] = 1.0
    g1H, g2H = fullci.rdm(phi, M, 1), fullci.rdm(phi, M, 2)
    U1H = np.array([annihilate(phi, (p,))[S1] for p in range(M)]).T
    Qh, svh, Vhh = np.linalg.svd(U1H, full_matrices=False)
    kph = svh > 1e-10 * svh[0]
    r1H = int(kph.sum())
    P1H = Qh[:, kph]
    RcH = np.diag(svh[kph]) @ Vhh[kph, :]
    U2H = np.array([annihilate(phi, I)[S2] for I in az.tuples(M, 2)]).T
    Q2h, sv2h, _ = np.linalg.svd(U2H, full_matrices=False)
    kp2h = sv2h > 1e-10 * sv2h[0]
    P2H = Q2h[:, kp2h]
    ladH = np.zeros((M, int(kp2h.sum()), r1H), dtype=complex)
    for x in range(M):
        for al in range(r1H):
            b = np.zeros(1 << M, dtype=complex)
            b[S1] = P1H[:, al]
            ladH[x, :, al] = P2H.conj().T @ fullci.aa(x, b)[S2]
    GH = [[ladH[x].conj().T @ ladH[y] for y in range(M)] for x in range(M)]
    g2HL, _ = gamma2_from_ladder(RcH, GH, M)
    e_mf = energy(g1H, g2HL, h, V, enuc)
    print("    V3  mean-field ladder: r_1 = %d (= N), chi' = %d (= C(N,2) = %d)"
          % (r1H, int(kp2h.sum()), N * (N - 1) // 2))
        # gamma_2 of a determinant is exactly gamma_1 wedge gamma_1
    print("        |gamma_2^MF - gamma_1 wedge gamma_1| = %.2e   E = %.9f  dev from E_HF %.2e"
          % (np.abs(g2HL - az.wedge(g1H, g1H)).max(), e_mf, abs(e_mf - e_hf)))

    # B1  the rank bound on achievable accuracy
    print()
    print("    B1  LOWER BOUND on the energy error of any ladder at a given chi':")
    print("        best rank-chi' approximation to the exact Mat_2, by SVD truncation")
    Mat2 = az.tensor_to_coord(g2, M, 2) * 2.0
    w, Uw = np.linalg.eigh(0.5 * (Mat2 + Mat2.conj().T))
    order = np.argsort(-np.abs(w))
    print("        %6s %14s %14s %10s" % ("chi'", "E(trunc)", "E-E_FCI", "frac of Ecorr"))
    for chi in sorted(set([1, 2, 3, 4, 6, 8, 10, 14, 20, chi_ex, len(w)])):
        if chi > len(w):
            continue
        idx = order[:chi]
        Mt = (Uw[:, idx] * w[idx]) @ Uw[:, idx].conj().T
        g2t = az.coord_to_tensor(Mt / 2.0, M, 2)
        et = energy(g1, g2t, h, V, enuc)
        print("        %6d %14.9f %14.6f %9.1f%%"
              % (chi, et, et - e_fci_ref,
                 100.0 * (et - e_hf) / (e_fci_ref - e_hf) if abs(e_fci_ref - e_hf) > 1e-12 else 0.0))
    print()
