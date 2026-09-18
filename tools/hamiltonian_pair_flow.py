#!/usr/bin/env python3
r"""Route F, episode 3: the parameter-free Hamiltonian pair flow  Xdot = i[grad_X E, X].

Symbols as in docs/COMMUTATOR-STRUCTURE.md. X = Mat_2 over sorted pairs, X[(p,v),(q,w)] =
2 gamma_2[pvqw]. Physical energy

    E[X] = sum_pq h[p,q] gamma_1[p,q] + (1/4) sum_pvqw V[p,v,q,w] gamma_2[p,v,q,w],
    gamma_1 = Tr_2 gamma_2 / (N-1),

a LINEAR functional of X, so its gradient G = grad_X E (the Hermitian pair-space matrix with
dE = Re tr(G^dag dX)) is state-independent. It is K_2 with different weights on the h- and V-parts;
both weights are determined here numerically by matching dE against tr(G^dag dX) on random dX, and
then checked in closed form.

With the u(C(M,2)) Lie-Poisson bracket {F, G}(X) = i tr(X [grad F, grad G]) the Hamiltonian flow is

    Xdot_H = i [G, X]

which conserves the spectrum of X (positivity structural), E (by antisymmetry of the bracket), and
has real linearised frequencies at any stable point. No closure parameter enters. The question is
only its accuracy. Reported per state: |Xdot_H - Xdot_exact| / |Xdot_exact| against the same for the
drift D = i[K_2, X] alone and for D + Phi[gamma_3^recon] with lambda_3 = 0 and Mueller-2405 + Lackner;
dE/dt of each; and the induced 1-RDM flow error |c(Xdot) - gamma_1-dot| / |gamma_1-dot|,
c(T) = Tr_2 T / (N-1).

EPISODE-3 IDENTITY (measured 1e-16): c(i[K_h, X]) = i[h, gamma_1] and c(Phi) = (N-2) c(i[K_V, X]), so
gamma_1-dot = i[h, gamma_1] + (N-1) c(i[K_V, X]) exactly and c(i[grad E, X]) = gamma_1-dot / (2(N-1)).
Hence K* = K_h + (N-1) K_V = 2(N-1) grad E gives a commutator flow with EXACT induced rank-1 EOM and
exact energy conservation; its rank-2 accuracy is reported alongside.

Also reported, because it is the first thing anyone will ask: does the Hamiltonian pair flow
reproduce TDHF at a determinant? At the HF determinant (lambda_2 = 0) the exact Xdot is the exact
t = 0 rate, D + Phi[gamma_3 = g1^g1^g1]; we print |Xdot_H - Xdot_exact| there too.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
import commutator_structure as cs
from commutator_k3 import contract, from_mat2, pair_k2
from qdiss import ansatz as az
from qdiss import muller2405
from qdiss.eom import fullci, wick


def energy(X, h, V, M, N):
    T = from_mat2(X, M)
    g1 = np.einsum("pvqv->pq", T) / (N - 1)
    return (np.einsum("pq,pq", h, g1) + 0.25 * np.einsum("pvqw,pvqw", V, T)).real


def grad_E(h, V, M, N):
    """Hermitian G with dE = Re tr(G^dag dX), found from the h- and V-parts of K_2 by two scalars."""
    P = cs.pairs(M)
    Kh = pair_k2(h, 0 * V, M)
    Kv = pair_k2(0 * h, V, M)
    rng = np.random.default_rng(0)
    A = rng.normal(size=(len(P), len(P))) + 1j * rng.normal(size=(len(P), len(P)))
    dX = A + A.conj().T
    dE_h = energy(dX, h, 0 * V, M, N)
    dE_v = energy(dX, 0 * h, V, M, N)
    ch = dE_h / np.trace(Kh.conj().T @ dX).real
    cv = dE_v / np.trace(Kv.conj().T @ dX).real
    G = ch * Kh + cv * Kv
    # verify on a second random direction
    B = rng.normal(size=(len(P), len(P))) + 1j * rng.normal(size=(len(P), len(P)))
    dX2 = B + B.conj().T
    err = abs(energy(dX2, h, V, M, N) - np.trace(G.conj().T @ dX2).real)
    return G, ch, cv, err


def main():
    for atom, N, spin in (cs.H3, cs.H4):
        tag = "H3 bent" if N == 3 else "H4 trapezoid"
        first = True
        for label, dm, dets, M, h, V in cs.states(atom, N, spin, tag, quench_times=(1.0, 2.0)):
            H = fullci.build_H(h, V, M)
            (g1, g2, g3, d1, d2), rk = cs.ensemble_rdms(dm, dets, M, H)
            l2 = az.cumulant(2, g2, {("G", 1): g1})
            EYE = np.eye(M, dtype=complex)
            X = cs.to_mat2(g2, M)
            if first:
                G, ch, cv, err = grad_E(h, V, M, N)
                print("=== %s: grad_X E = %.6f * K_2(h-part) + %.6f * K_2(V-part) ; linearity check %.1e"
                      % (tag, ch, cv, err))
                print("    expected 1/(2(N-1)) = %.6f for h, 1/2 for V" % (1.0 / (2 * (N - 1))))
                first = False
            tab = {("D", 1): EYE, ("G", 1): g1, ("G", 2): g2, ("G", 3): np.zeros((M,) * 6, dtype=complex)}
            Dm = cs.to_mat2(wick.gamma_dot(2, h, V, tab, M), M)
            Kh, Kv = pair_k2(h, 0 * V, M), pair_k2(0 * h, V, M)
            Ks = Kh + (N - 1) * Kv                     # = 2(N-1) grad E ; exact contraction, conserves E
            flows = {"Hamiltonian pair flow i[grad E, X]": 1j * (G @ X - X @ G),
                     "rescaled K* = K_h + (N-1) K_V": 1j * (Ks @ X - X @ Ks),
                     "drift D = i[K_2, X] alone": Dm}
            for name, g3r in (("lambda_3 = 0", az.separable(3, {("G", 1): g1, ("L", 2): l2})),
                              ("muller2405 + Lackner", muller2405.gamma3(g1, g2, N, True))):
                tab[("G", 3)] = g3r
                flows["D + Phi[%s]" % name] = cs.to_mat2(wick.gamma_dot(2, h, V, tab, M), M)
            Xd = cs.to_mat2(d2, M)
            print("--- %s  rank(rho)~%d  |Xdot_exact|=%.4f  |gamma_1-dot|=%.4f" % (label, rk, np.linalg.norm(Xd), np.linalg.norm(d1)))
            for name, F in flows.items():
                dE = energy(F, h, V, M, N)
                e2 = np.linalg.norm(F - Xd) / np.linalg.norm(Xd)
                c = contract(F, M, N)
                e1 = np.linalg.norm(c - d1) / max(np.linalg.norm(d1), 1e-30)
                U = np.linalg.eigh(0.5 * (X + X.conj().T))[1]
                spec = np.linalg.norm(np.diag(U.conj().T @ F @ U)) / np.linalg.norm(F)
                print("    %-38s rel err rank-2 %.4f   rel err induced rank-1 %.4f   dE/dt %+.1e   spectral frac %.4f"
                      % (name, e2, e1, dE, spec))
            print()


if __name__ == "__main__":
    main()
