#!/usr/bin/env python3
r"""Route F, episode 4d: the geminal mean-field flow, derived in closed form, tested as a closure.

Symbols (docs/COMMUTATOR-STRUCTURE.md). Pair-space matrices over sorted pairs:
    X = Mat_2 = 2 gamma_2,  Xd = Mat(g1 ^ g1),  L = X - Xd = Mat(lambda_2),  N = g1 (+) g1,
    K_h = h (+) h,  K_V[(pv),(qw)] = V[pvqw],  K_WT = W^T (+) W^T with W[p,q] = sum_rs V[p,r,q,s] g1[r,s].
Identity (tools/phi_identify.py, residual 1e-14):

    Phi[A_3(g1 g1 g1)]  =  i ( S^dag Xd - Xd S ),     S = K_WT - K_V N .

Hence with the drift D = i[K_h + K_V, X] the lambda_2 = lambda_3 = 0 flow is
    D + Phi_ggg  =  i[K_h + K_V, X]  +  i(S^dag Xd - Xd S).
Since the exact-flow object is X, two closures suggest themselves:
    GMF-d : Xdot = i[K_h + K_V, X] + i(S^dag Xd - Xd S)        (as derived; Xd = X - L recomputed each step)
    GMF-x : Xdot = i(S_tot^dag X - X S_tot),  S_tot = K_h + K_WT + K_V (I - N)   (sandwich in X: XR-structural)
GMF-x differs from GMF-d by -i(S^dag L - L S). GMF-x preserves X >= 0 by construction (XR, Rdot = -i R S_tot);
its K/Gamma split is  K = K_h + K_WT + K_V - (K_V N + N K_V)/2,  Gamma = i[N, K_V].

Reported per state: the identity check; rank-2 error of GMF-d, GMF-x, D, lambda_3 = 0; contraction error;
dE/dt; spectral fraction; and tr(Gamma X) (trace conservation requires it to vanish).
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
import commutator_structure as cs
from commutator_k3 import contract, pair_k2
from hamiltonian_pair_flow import energy
from phi_anatomy import spectral_split
from qdiss import ansatz as az
from qdiss.eom import fullci, wick


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
            X = cs.to_mat2(g2, M)
            Xd = cs.to_mat2(az.wedge(g1, g1), M)
            L = X - Xd
            Xdot = cs.to_mat2(d2, M)
            tab = {("D", 1): EYE, ("G", 1): g1, ("G", 2): g2, ("G", 3): np.zeros((M,) * 6, dtype=complex)}
            Dm = cs.to_mat2(wick.gamma_dot(2, h, V, tab, M), M)
            tab[("G", 3)] = az.separable(3, {("G", 1): g1, ("L", 2): 0 * l2})
            Phi_ggg = cs.to_mat2(wick.gamma_dot(2, h, V, tab, M), M) - Dm
            tab[("G", 3)] = az.separable(3, {("G", 1): g1, ("L", 2): l2})
            Xl = cs.to_mat2(wick.gamma_dot(2, h, V, tab, M), M)
            Np = pair_k2(g1, 0 * V, M)
            KV = pair_k2(0 * h, V, M)
            Kh = pair_k2(h, 0 * V, M)
            W = np.einsum("prqs,rs->pq", V, g1)
            KWT = pair_k2(W.T, 0 * V, M)
            S = KWT - KV @ Np
            ident = np.abs(1j * (S.conj().T @ Xd - Xd @ S) - Phi_ggg).max()
            Stot = Kh + KWT + KV @ (np.eye(len(Np)) - Np)
            gmf_d = Dm + 1j * (S.conj().T @ Xd - Xd @ S)
            gmf_x = 1j * (Stot.conj().T @ X - X @ Stot)
            Gam = 1j * (Np @ KV - KV @ Np)
            print("=== %s  rank(rho)~%d  |Xdot|=%.4f   identity |Phi_ggg - i(S^dag Xd - Xd S)| = %.1e   |Gamma|_2=%.3f  tr(Gamma X)=%+.1e"
                  % (label, rk, np.linalg.norm(Xdot), ident, np.linalg.norm(Gam, 2), np.trace(Gam @ X).real))
            for nm, F in (("D alone", Dm), ("GMF-d  D + i(S^dag Xd - Xd S)", gmf_d),
                          ("GMF-x  i(Stot^dag X - X Stot)", gmf_x), ("lambda_3 = 0", Xl), ("exact", Xdot)):
                sp = spectral_split(X, F)[0]
                print("    %-34s rank-2 err %.4f   contraction err %.4f   dE/dt %+.1e   spectral %.4f   tr Xdot %+.1e"
                      % (nm, np.linalg.norm(F - Xdot) / np.linalg.norm(Xdot),
                         np.linalg.norm(contract(F, M, N) - d1) / max(np.linalg.norm(d1), 1e-30),
                         energy(F, h, V, M, N), np.linalg.norm(sp) / np.linalg.norm(F), np.trace(F).real))
            print()


if __name__ == "__main__":
    main()
