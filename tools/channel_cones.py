#!/usr/bin/env python3
r"""Route F, episode 6a: which piece of the rank-2 flow threatens which cone (D, Q, G)?

Symbols (docs/COMMUTATOR-STRUCTURE.md). gamma_2[p,v,q,w] = <a+_p a+_v a_w a_q>, gamma_1[p,q] = <a+_p a_q>.
The three Gram matrices (all PSD for a representable state):
    D[(pv),(qw)] = <a+_p a+_v a_w a_q>                     = gamma_2[p,v,q,w]          (pp; X = Mat_2 = 2D on sorted pairs)
    Q[(pv),(qw)] = <a_p a_v a+_w a+_q>                     = d_vw d_pq - d_vq d_pw - d_vw g1[q,p] + d_vq g1[w,p]
                                                              + d_pw g1[q,v] - d_pq g1[w,v] + gamma_2[w,q,v,p]   (hh)
    G[(pq),(rs)] = <a+_p a_q a+_s a_r>                     = d_qs g1[p,r] - gamma_2[p,s,r,q]                     (ph, M^2 x M^2)
Each is verified here against the Gram matrix built directly from the (ensemble of) wavefunction(s).
Because Q and G are affine in (gamma_1, gamma_2), any rank-2 flow gamma_2-dot (with gamma_1-dot by contraction)
induces Qdot and Gdot linearly.

The lambda_3 = 0 flow in closed form (episode 5c) is
    drift  i[K_h + K_V, X]  +  pp sandwich  i(S^dag X - X S)  +  i(Q_C - Q_C^dag)  +  Phi_ph      (Q_C = C(x)g1 + g1(x)C)
and the exact flow adds the lambda_3 remainder R3 = Phi[lambda_3]. For each piece and each cone with a
kernel we report the LEAK  max_{v in ker M} |v^dag Mdot v| / |Mdot|_F : zero means the piece is tangent to
that cone at this state (the exact total must be tangent to all three, XR theorem). The kernel is taken at
eigenvalue < 1e-8. For thermal states there is no kernel and only |Mdot| is reported.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
import commutator_structure as cs
from commutator_k3 import from_mat2
from closed_form_sep import closed_form
from qdiss import ansatz as az
from qdiss.eom import fullci, wick

KER = 1e-8


def q_matrix(g1, g2, M):
    I = np.eye(M)
    T = (np.einsum("vw,pq->pvqw", I, I) - np.einsum("vq,pw->pvqw", I, I)
         - np.einsum("vw,qp->pvqw", I, g1) + np.einsum("vq,wp->pvqw", I, g1)
         + np.einsum("pw,qv->pvqw", I, g1) - np.einsum("pq,wv->pvqw", I, g1)
         + g2.transpose(3, 2, 1, 0))                 # gamma_2[w,q,v,p]
    return cs.to_mat2(T, M) / 2.0                     # sorted pairs, no factor 2


def g_matrix(g1, g2, M):
    I = np.eye(M)
    T = np.einsum("qs,pr->pqrs", I, g1) - g2.transpose(0, 3, 2, 1)   # gamma_2[p,s,r,q] -> index order p,q,r,s
    return T.reshape(M * M, M * M)


def gram_direct(dm, dets, M, kind):
    """Ensemble-averaged Gram matrices from the wavefunctions: kind in {'D','Q','G'}."""
    w, U = np.linalg.eigh(0.5 * (dm + dm.conj().T))
    P = cs.pairs(M)
    out = 0
    for i in np.where(w > 1e-12)[0]:
        psi = np.zeros(1 << M, dtype=complex)
        psi[dets] = U[:, i]
        if kind == "D":
            vecs = [fullci.aa(v, fullci.aa(p, psi)) for (p, v) in P]           # a_v a_p |psi>
        elif kind == "Q":
            vecs = [fullci.ad(v, fullci.ad(p, psi)) for (p, v) in P]           # a+_v a+_p |psi>
        else:
            vecs = [fullci.ad(q, fullci.aa(p, psi)) for p in range(M) for q in range(M)]   # a+_q a_p |psi>
        Vm = np.array(vecs)
        out = out + w[i] * (Vm.conj() @ Vm.T)
    return out


def leak(Mmat, Mdot):
    x, U = np.linalg.eigh(0.5 * (Mmat + Mmat.conj().T))
    ker = U[:, x < KER]
    if ker.shape[1] == 0:
        return float("nan"), 0
    B = ker.conj().T @ (0.5 * (Mdot + Mdot.conj().T)) @ ker
    return np.abs(B).max() / max(np.linalg.norm(Mdot), 1e-30), ker.shape[1]


def main():
    for atom, N, spin in (cs.H3, cs.H4):
        tag = "H3" if N == 3 else "H4"
        for label, dm, dets, M, h, V in cs.states(atom, N, spin, tag, quench_times=(1.0,)):
            if "t=0" in label:
                continue
            H = fullci.build_H(h, V, M)
            (g1, g2, g3, d1, d2), rk = cs.ensemble_rdms(dm, dets, M, H)
            Dm_, Qm, Gm = cs.to_mat2(g2, M) / 2.0, q_matrix(g1, g2, M), g_matrix(g1, g2, M)
            chk = [np.abs(A - gram_direct(dm, dets, M, k)).max() for A, k in ((Dm_, "D"), (Qm, "Q"), (Gm, "G"))]
            mins = [np.linalg.eigvalsh(0.5 * (A + A.conj().T)).min() for A in (Dm_, Qm, Gm)]
            print("=== %s  rank(rho)~%d   D,Q,G vs direct Gram: %.1e %.1e %.1e   min eig: %.1e %.1e %.1e"
                  % (label, rk, *chk, *mins))
            drift, sandwich, qterm, Phi_ph, parts = closed_form(h, V, g1, g2, M)
            Xdot = cs.to_mat2(d2, M)
            R3 = Xdot - (drift + sandwich + qterm + Phi_ph)
            pieces = {"drift i[K_2,X]": drift, "pp sandwich i(S^dag X - X S)": sandwich, "i(Q_C - Q_C^dag)": qterm,
                      "Phi_ph (ring)": Phi_ph, "lambda_3 remainder": R3, "lambda_3 = 0 total": drift + sandwich + qterm + Phi_ph,
                      "exact total": Xdot}
            print("    %-30s  %-26s %-26s %-26s" % ("piece", "D: leak (n_ker) |Ddot|", "Q: leak (n_ker) |Qdot|", "G: leak (n_ker) |Gdot|"))
            for nm, F in pieces.items():
                g2dot = from_mat2(F, M)
                g1dot = np.einsum("pvqv->pq", g2dot) / (N - 1)
                Ddot = F / 2.0
                Qdot = q_matrix(g1dot, g2dot, M) - q_matrix(0 * g1, 0 * g2, M)     # affine part removed
                Gdot = g_matrix(g1dot, g2dot, M)
                row = []
                for Mm, Md in ((Dm_, Ddot), (Qm, Qdot), (Gm, Gdot)):
                    lk, nk = leak(Mm, Md)
                    row.append("%.1e (%2d) %.3f" % (lk, nk, np.linalg.norm(Md)) if nk else "  ---  (0)  %.3f" % np.linalg.norm(Md))
                print("    %-30s  %-26s %-26s %-26s" % (nm, *row))
            print()


if __name__ == "__main__":
    main()
