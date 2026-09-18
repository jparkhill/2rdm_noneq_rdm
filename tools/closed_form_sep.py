#!/usr/bin/env python3
r"""Route F, episode 5c: the whole lambda_3 = 0 rank-2 flow in closed algebraic form, verified directly.

Symbols (docs/COMMUTATOR-STRUCTURE.md). Pair space: X = Mat_2, N = g1 (+) g1, K_h, K_V, K_WT = W^T (+) W^T,
L = Mat(lambda_2), Xd = X - L. One-body C[p,q] = sum_rtu V[r,t,p,u] l2[r,t,q,u] (rank-1 collision matrix);
K_C = C (+) C, K_gC = (g1 C) (+) (g1 C), K_Cg = (C g1) (+) (C g1).
ph space (M^2 x M^2, pairs ordered (ket, bra)): ph(T)[(p,q),(v,w)] = T[p,v,q,w]; Vd = ph(V), Ld = ph(lambda_2);
Sw swaps the two entries of a pair index; gI = g1 (x) I, gTI = g1^T (x) I, Ig = I (x) g1, IgT = I (x) g1^T.
back(A)[p,v,q,w] = A[(p,q),(v,w)], antisymmetrised over (p,v) and (q,w), then to Mat_2 (x2), /2.

Identified in tools/phi_identify.py (Phi_ggg) and tools/phi_gl_sparse.py (Phi_gl, held-out 8e-15):

    Phi_ggg = i(S^dag Xd - Xd S),                                S = K_WT - K_V N
    Phi_gl  = i(S^dag L  - L  S) + i(Q - Q^dag) + Phi_ph,          Q = (1/2){N, K_C} - K_gC - K_Cg   (wait: see code; the
              identified combination is (1/2)(N K_C - K_gC) + (1/2)(K_C N - K_Cg) in the i(P - P^dag) form)
    Phi_ph  = 2 i [ (P1 + P2 - P3 - P4) - h.c. ],  P1 = Sw gTI Vd Ld,  P4 = Sw Ig Vd Ld,
                                                   P2 = Ld Vd (Sw gI)^dag,  P3 = Ld Vd (Sw IgT)^dag
so that, with the drift D = i[K_h + K_V, X],

    D + Phi_sep  =  i[K_h + K_V, X] + i(S^dag X - X S) + i(Q - Q^dag) + Phi_ph .

This file evaluates the right-hand side from (h, V, gamma_1, gamma_2) alone and compares it with the wick
term list evaluated at gamma_3 = separable_3(gamma_1, lambda_2), at eight exact states. No fitting.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
import commutator_structure as cs
from commutator_k3 import pair_k2
from phi_identify_ph import ph_d, ph_d_inv, antisym4
from qdiss import ansatz as az
from qdiss.eom import fullci, wick


def closed_form(h, V, g1, g2, M):
    l2 = az.cumulant(2, g2, {("G", 1): g1})
    X = cs.to_mat2(g2, M)
    L = cs.to_mat2(l2, M)
    Np = pair_k2(g1, 0 * V, M)
    Kh, KV = pair_k2(h, 0 * V, M), pair_k2(0 * h, V, M)
    W = np.einsum("prqs,rs->pq", V, g1)
    KWT = pair_k2(W.T, 0 * V, M)
    C = np.einsum("rtpu,rtqu->pq", V, l2)
    KC, KgC, KCg = pair_k2(C, 0 * V, M), pair_k2(g1 @ C, 0 * V, M), pair_k2(C @ g1, 0 * V, M)
    S = KWT - KV @ Np
    Q = 0.5 * (Np @ KC - KgC) + 0.5 * (KC @ Np - KCg)
    # ph channel
    I2 = np.eye(M, dtype=complex)
    Vd, Ld = ph_d(V, M), ph_d(l2, M)
    Sw = np.eye(M * M, dtype=complex).reshape(M, M, M, M).transpose(0, 1, 3, 2).reshape(M * M, M * M)
    gI, gTI, Ig, IgT = np.kron(g1, I2), np.kron(g1.T, I2), np.kron(I2, g1), np.kron(I2, g1.T)
    P1 = Sw @ gTI @ Vd @ Ld
    P4 = Sw @ Ig @ Vd @ Ld
    P2 = Ld @ Vd @ (Sw @ gI).conj().T
    P3 = Ld @ Vd @ (Sw @ IgT).conj().T
    back = lambda A: cs.to_mat2(antisym4(ph_d_inv(A, M)), M) / 2.0
    Pph = back(P1 + P2 - P3 - P4)
    Phi_ph = 2j * (Pph - Pph.conj().T)
    drift = 1j * ((Kh + KV) @ X - X @ (Kh + KV))
    sandwich = 1j * (S.conj().T @ X - X @ S)
    qterm = 1j * (Q - Q.conj().T)
    return drift, sandwich, qterm, Phi_ph, dict(S=S, Q=Q, Gamma=1j * (S - S.conj().T))


def main():
    for atom, N, spin in (cs.H3, cs.H4):
        tag = "H3" if N == 3 else "H4"
        for label, dm, dets, M, h, V in cs.states(atom, N, spin, tag, quench_times=(1.0, 2.0)):
            if "t=0" in label:
                continue
            H = fullci.build_H(h, V, M)
            (g1, g2, g3, d1, d2), rk = cs.ensemble_rdms(dm, dets, M, H)
            l2 = az.cumulant(2, g2, {("G", 1): g1})
            EYE = np.eye(M, dtype=complex)
            tab = {("D", 1): EYE, ("G", 1): g1, ("G", 2): g2, ("G", 3): az.separable(3, {("G", 1): g1, ("L", 2): l2})}
            ref = cs.to_mat2(wick.gamma_dot(2, h, V, tab, M), M)
            drift, sandwich, qterm, Phi_ph, parts = closed_form(h, V, g1, g2, M)
            tot = drift + sandwich + qterm + Phi_ph
            Xd_exact = cs.to_mat2(d2, M)
            print("%-42s |closed form - wick(lambda_3=0)| = %.1e   |wick|=%.3f   pieces: drift %.3f  sandwich %.3f  Q %.3f  ph %.3f   |Gamma|_2=%.3f   rel err vs exact %.4f"
                  % (label, np.abs(tot - ref).max(), np.linalg.norm(ref), np.linalg.norm(drift), np.linalg.norm(sandwich),
                     np.linalg.norm(qterm), np.linalg.norm(Phi_ph), np.linalg.norm(parts["Gamma"], 2),
                     np.linalg.norm(tot - Xd_exact) / np.linalg.norm(Xd_exact)))


if __name__ == "__main__":
    main()
