#!/usr/bin/env python3
r"""Route F, episode 7, step 1-2: validation of the F0 + delta closure before any trajectory.

Checks (symbols in tools/f0delta.py):
  V1  exact flow is a fixed point of the correction: |delta(Xdot_exact)| / |Xdot_exact| = O(sigma).
  V2  V = 0: F0 is exact and delta = 0.
  V3  at the eight standard states: kernel leaks of F0 + delta on D, Q, G (must be ~0 where kernels exist),
      contraction error vs the exact rank-1 EOM (must be ~1e-14), dE/dt (must be ~1e-14), rank-2 error vs exact,
      alongside F0, muller2405 + Lackner, F0 + delta_D (D only) and muller + delta.
  V4  sigma / mu sensitivity on one state.
  V5  extensivity: H3 x H2 product state (block-diagonal h, V); cross-fragment norm of F0 + delta vs exact (0).
  Timing per flow evaluation.
"""
import sys
import time
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
import commutator_structure as cs
from commutator_k3 import from_mat2
from channel_cones import q_matrix, g_matrix, KER
from hamiltonian_pair_flow import energy
from f0delta import F0Delta
from molecular_2rdm import build_H_sector, molecule, sector
from qdiss.eom import fullci


def leaks(g1, g2, F, M, N):
    g2F = from_mat2(F, M)
    g1F = np.einsum("pvqv->pq", g2F) / (N - 1)
    out = []
    for Pm, Pd in ((cs.to_mat2(g2, M) / 2, F / 2),
                   (q_matrix(g1, g2, M), q_matrix(g1F, g2F, M) - q_matrix(0 * g1, 0 * g2, M)),
                   (g_matrix(g1, g2, M), g_matrix(g1F, g2F, M))):
        x, U = np.linalg.eigh(0.5 * (Pm + Pm.conj().T))
        K = U[:, x < KER]
        out.append(np.abs(K.conj().T @ Pd @ K).max() / np.linalg.norm(Pd) if K.shape[1] else float("nan"))
    return out


def main():
    for atom, N, spin in (cs.H3, cs.H4):
        tag = "H3" if N == 3 else "H4"
        first = True
        for label, dm, dets, M, h, V in cs.states(atom, N, spin, tag, quench_times=(1.0, 2.0)):
            if "t=0" in label:
                continue
            H = fullci.build_H(h, V, M)
            (g1, g2, g3, d1, d2), rk = cs.ensemble_rdms(dm, dets, M, H)
            cl = F0Delta(h, V, M, N)
            Xd = cs.to_mat2(d2, M)
            if first:
                print("=== %s  M=%d N=%d  n2=%d  rank(C)=%d  null dim=%d" % (tag, M, N, cl.n2, cl.rankC, cl.Z.shape[1]))
                # V1 fixed point, V2 V=0
                dx = cl.delta(g1, g2, Xd)
                print("  V1 |delta(exact)|/|Xdot| = %.2e  (sigma=%.0e)" % (np.linalg.norm(dx) / np.linalg.norm(Xd), cl.sigma))
                cl0 = F0Delta(h, 0 * V, M, N)
                H0 = fullci.build_H(h, 0 * V, M)
                (_, _, _, _, d20), _ = cs.ensemble_rdms(dm, dets, M, H0)
                F00 = cl0.F0(g1, g2)
                print("  V2 V=0: |F0 - exact| = %.1e   |delta| = %.1e" % (np.abs(F00 - cs.to_mat2(d20, M)).max(), np.linalg.norm(cl0.delta(g1, g2, F00))))
                # V4 sensitivity
                F = cl.F0(g1, g2)
                for sg, mu in ((1e-2, 1e2), (1e-3, 1e2), (1e-2, 1e4), (1e-1, 1e2)):
                    clx = F0Delta(h, V, M, N, sigma=sg, mu=mu)
                    Fx = F + clx.delta(g1, g2, F)
                    lk = leaks(g1, g2, Fx, M, N)
                    print("  V4 sigma=%.0e mu=%.0e : rank-2 err %.4f  leaks D %.1e Q %.1e G %.1e  |delta|=%.4f  |delta(exact)|/|Xdot|=%.1e"
                          % (sg, mu, np.linalg.norm(Fx - Xd) / np.linalg.norm(Xd), *lk, np.linalg.norm(Fx - F),
                             np.linalg.norm(clx.delta(g1, g2, Xd)) / np.linalg.norm(Xd)))
                first = False
            print("--- %s  rank(rho)~%d  |Xdot|=%.4f" % (label, rk, np.linalg.norm(Xd)))
            d1x = cl.exact_g1dot(g1, g2)
            print("    [check] exact rank-1 EOM from gamma_2 vs fullci: %.1e" % np.abs(d1x - d1).max())
            for mode in ("f0", "muller", "f0delta_D", "f0delta", "muller_delta"):
                t0 = time.time()
                F = cl.flow(g2, mode)
                dt = time.time() - t0
                g2F = from_mat2(F, M)
                cerr = np.linalg.norm(np.einsum("pvqv->pq", g2F) / (N - 1) - d1x) / np.linalg.norm(d1x)
                lk = leaks(g1, g2, F, M, N)
                print("    %-13s rank-2 err %.4f   contraction err %.1e   dE/dt %+.1e   leaks D %.1e Q %.1e G %.1e   herm %.0e   %.2fs"
                      % (mode, np.linalg.norm(F - Xd) / np.linalg.norm(Xd), cerr, energy(F, h, V, M, N), *lk,
                         np.abs(F - F.conj().T).max(), dt))
        print()

    # V5 extensivity: H3 (M=6,N=3) x H2 (M=4,N=2) with block-diagonal integrals
    molA, mfA, MA, hA, VA, _ = molecule(cs.H3[0], "sto-3g", 0, 1)
    molB, mfB, MB, hB, VB, _ = molecule("H 0 0 0; H 0.74 0 0", "sto-3g", 0, 0)
    M, N = MA + MB, 3 + 2
    h = np.zeros((M, M), dtype=complex); h[:MA, :MA] = hA; h[MA:, MA:] = hB
    V = np.zeros((M,) * 4, dtype=complex); V[:MA, :MA, :MA, :MA] = VA; V[MA:, MA:, MA:, MA:] = VB
    dets = sector(M, N)
    # product state: kicked ground state of each fragment
    def frag_state(hf, Vf, Mf, Nf, mol, mf):
        d = sector(Mf, Nf)
        es, vec = np.linalg.eigh(build_H_sector(hf, Vf, Mf, d))
        import scipy.linalg as sla
        Uk = sla.expm(-1j * 0.3 * cs.sector_operator(cs.dipole(mol, mf, Mf), Mf, d))
        psi = np.zeros(1 << Mf, dtype=complex); psi[d] = Uk @ vec[:, 0]
        return psi
    pA = frag_state(hA, VA, MA, 3, molA, mfA)
    pB = frag_state(hB, VB, MB, 2, molB, mfB)
    # tensor product in occupation-number encoding: fragment A occupies bits 0..MA-1, B bits MA..M-1
    psi = np.zeros(1 << M, dtype=complex)
    for a in np.nonzero(pA)[0]:
        for b in np.nonzero(pB)[0]:
            # fermionic sign for interleaving: B operators pass over A's electrons; ordering convention of
            # fullci is by bit position, so a plain product is consistent when all B bits are above all A bits.
            psi[a | (b << MA)] += pA[a] * pB[b]
    H = fullci.build_H(h, V, M)
    g1 = fullci.rdm(psi, M, 1); g2 = fullci.rdm(psi, M, 2); d2 = fullci.rdm_dot(psi, H, M, 2)
    cl = F0Delta(h, V, M, N)
    Xd = cs.to_mat2(d2, M)
    P = cs.pairs(M)
    cross = np.array([[(p < MA) != (v < MA) or (q < MA) != (w < MA) or ((p < MA) != (q < MA))
                       for (q, w) in P] for (p, v) in P])
    # extensivity is a statement about the CUMULANT rate: dlambda_2/dt = dgamma_2/dt - (dg1 ^ g1 + g1 ^ dg1) must have
    # zero cross-fragment elements (the gamma_2 rate itself has gamma_1^A ^ gamma_1^B cross terms that move).
    from qdiss import ansatz as az
    def dl2(F):
        g2F = from_mat2(F, M); g1F = np.einsum("pvqv->pq", g2F) / (N - 1)
        return cs.to_mat2(g2F - az.wedge(g1F, g1) - az.wedge(g1, g1F), M)
    print("=== V5 extensivity H3 x H2, M=%d N=%d   |Xdot_exact|=%.4f  exact cross-fragment norm of dlambda_2/dt %.1e"
          % (M, N, np.linalg.norm(Xd), np.linalg.norm(dl2(Xd)[cross])))
    for mode in ("f0", "f0delta"):
        t0 = time.time(); F = cl.flow(g2, mode); dt = time.time() - t0
        print("    %-9s cross-fragment norm of dlambda_2/dt %.2e   (|dlambda_2/dt| = %.3f)   rank-2 err %.4f   %.1fs"
              % (mode, np.linalg.norm(dl2(F)[cross]), np.linalg.norm(dl2(F)), np.linalg.norm(F - Xd) / np.linalg.norm(Xd), dt))


if __name__ == "__main__":
    main()
