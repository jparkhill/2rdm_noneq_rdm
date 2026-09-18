#!/usr/bin/env python3
r"""Route F, episode 2: can the collision part be a commutator, and what does that cost the 1-RDM?

Symbols (as in docs/COMMUTATOR-STRUCTURE.md): X = Mat_2, X[(p,v),(q,w)] = 2 gamma_2[pvqw] over sorted
pairs; d/dt gamma_2 = D + Phi with D the gamma_2-carrying Wick terms (= i[K_2, X], K_2 Hermitian on
pair space, measured in episode 1) and Phi the gamma_3-carrying terms, linear in gamma_3.
x_i, U: eigenvalues and eigenvectors of X. "block-diagonal" = entries with |x_i - x_j| < DEG_TOL.

TWO QUESTIONS.

(i) BOUNDEDNESS. Split Phi = Phi_spec + Phi_off (block-diagonal and off-block parts in the eigenbasis
    of X). Phi_off is representable as i[K_3, X] with the UNIQUE off-block Hermitian
        K_3[ij] = Phi_off[ij] / ( i (x_j - x_i) )        (block-diagonal of K_3 free; set to 0)
    Report |K_3|_2 and |K_3|_F against |K_2|, and the same for the Tikhonov-regularised
        K_3^sigma[ij] = Phi_off[ij] (x_j - x_i) / ( i ((x_j - x_i)^2 + sigma^2) )
    with the fraction of Phi_off it reproduces.

(ii) CONTRACTION. gamma_1 = Tr_2 gamma_2 / (N-1), so any rank-2 flow induces a rank-1 flow
        c(T)[p,q] = sum_v T[p,v,q,v] / (N-1)     (T the 4-index form of a pair-space matrix).
    The exact rank-1 EOM has no gamma_3, so c(Phi[gamma_3]) depends on gamma_3 only through
    Tr_3 gamma_3 = (N-2) gamma_2 and is exactly the rank-1 collision term C[lambda_2]. A commutator
    closure keeps c(D) + c(Phi_off) and DROPS c(Phi_spec). Report |c(Phi_spec)| / |c(Phi)| and
    |c(Phi_spec)| / |gamma_1-dot|: the unavoidable error in the induced 1-RDM flow (equivalently in
    the rank-2 trace identity Tr_2 lambda_2 = -gamma_1 eta_1, which is the contraction relation
    rewritten). Also the spectral fraction of each induced rank-1 flow w.r.t. gamma_1.

Ground truth qdiss.eom.fullci; wick vs fullci checked on every row.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
import commutator_structure as cs
from qdiss import ansatz as az
from qdiss.eom import fullci, wick

DEG_TOL = 1e-6
SIGMAS = (1e-2, 1e-1)


def from_mat2(X, M):
    P = cs.pairs(M)
    T = np.zeros((M,) * 4, dtype=complex)
    for I, (p, v) in enumerate(P):
        for J, (q, w) in enumerate(P):
            t = X[I, J] / 2.0
            T[p, v, q, w] = t
            T[v, p, q, w] = -t
            T[p, v, w, q] = -t
            T[v, p, w, q] = t
    return T


def contract(X, M, N):
    return np.einsum("pvqv->pq", from_mat2(X, M)) / (N - 1)


def spectral_frac(A, Adot, tol=DEG_TOL):
    x, U = np.linalg.eigh(0.5 * (A + A.conj().T))
    Dm = U.conj().T @ (0.5 * (Adot + Adot.conj().T)) @ U
    B = np.where(np.abs(x[:, None] - x[None, :]) < tol, Dm, 0.0)
    n = np.linalg.norm(Adot)
    return np.linalg.norm(B) / n if n > 0 else 0.0


def pair_k2(h, V, M):
    """K_2 on pair space such that D = i[K_2, X]: h on each index plus V on the pair."""
    P = cs.pairs(M)
    K = np.zeros((len(P), len(P)), dtype=complex)
    for I, (p, v) in enumerate(P):
        for J, (q, w) in enumerate(P):
            K[I, J] = (h[p, q] * (v == w) + h[v, w] * (p == q) - h[p, w] * (v == q) - h[v, q] * (p == w)
                       + V[p, v, q, w])
    return K


for atom, N, spin in (cs.H3, cs.H4):
    tag = "H3 bent" if N == 3 else "H4 trapezoid"
    for label, dm, dets, M, h, V in cs.states(atom, N, spin, tag, quench_times=(1.0, 2.0)):
        if "HF determinant" in label:
            continue                       # spectral motion vanishes there by time reversal
        H = fullci.build_H(h, V, M)
        (g1, g2, g3, d1, d2), rk = cs.ensemble_rdms(dm, dets, M, H)
        EYE = np.eye(M, dtype=complex)
        X = cs.to_mat2(g2, M)
        tab = {("D", 1): EYE, ("G", 1): g1, ("G", 2): g2, ("G", 3): np.zeros((M,) * 6, dtype=complex)}
        Dm = cs.to_mat2(wick.gamma_dot(2, h, V, tab, M), M)
        tab[("G", 3)] = g3
        Xd = cs.to_mat2(wick.gamma_dot(2, h, V, tab, M), M)
        Phi = Xd - Dm
        print("=== %s   M=%d N=%d rank(rho)~%d ===" % (label, M, N, rk))
        print("    [check] wick vs fullci rank 2: %.1e ; contraction of exact Xdot vs fullci gamma_1-dot: %.1e"
              % (np.abs(Xd - cs.to_mat2(d2, M)).max(), np.abs(contract(Xd, M, N) - d1).max()))

        # K_2 in closed form, checked against the drift
        K2 = pair_k2(h, V, M)
        sgn = None
        for s in (+1, -1):
            r = np.abs(s * 1j * (K2 @ X - X @ K2) - Dm).max()
            if r < 1e-10:
                sgn = s
        print("    [check] D = %s i[K_2, X] in closed form: %s   |K_2|_2 = %.3f"
              % ("+" if sgn == 1 else "-" if sgn == -1 else "?", "yes" if sgn else "NO",
                 np.linalg.norm(K2, 2)))

        # (i) split Phi in the eigenbasis of X and solve for K_3
        x, U = np.linalg.eigh(0.5 * (X + X.conj().T))
        P = U.conj().T @ Phi @ U
        same = np.abs(x[:, None] - x[None, :]) < DEG_TOL
        Pspec = np.where(same, P, 0.0)
        Poff = P - Pspec
        dx = x[None, :] - x[:, None]
        with np.errstate(divide="ignore", invalid="ignore"):
            K3e = np.where(same, 0.0, Poff / (1j * dx))
        K3 = U @ K3e @ U.conj().T
        Phi_off = U @ Poff @ U.conj().T
        Phi_spec = U @ Pspec @ U.conj().T
        rep = np.abs(1j * (K3 @ X - X @ K3) - Phi_off).max()
        print("    Phi: |Phi|=%.4f  spectral frac=%.4f ; i[K_3,X] reproduces Phi_off to %.1e ; "
              "|K_3|_2=%.3f  |K_3|_F=%.3f  (|K_2|_2=%.3f)  min gap used=%.1e"
              % (np.linalg.norm(Phi), np.linalg.norm(Pspec) / np.linalg.norm(Phi), rep,
                 np.linalg.norm(K3, 2), np.linalg.norm(K3), np.linalg.norm(K2, 2),
                 np.abs(dx[~same]).min()))
        for sg in SIGMAS:
            K3s = U @ np.where(same, 0.0, Poff * dx / (1j * (dx ** 2 + sg ** 2))) @ U.conj().T
            got = 1j * (K3s @ X - X @ K3s)
            print("      sigma=%.0e : |K_3^sigma|_2=%.3f  reproduces %.4f of |Phi_off|"
                  % (sg, np.linalg.norm(K3s, 2),
                     1 - np.linalg.norm(got - Phi_off) / np.linalg.norm(Phi_off)))

        # (ii) contraction to the 1-RDM flow
        cD, cOff, cSpec, cPhi = (contract(A, M, N) for A in (Dm, Phi_off, Phi_spec, Phi))
        print("    induced gamma_1-dot: |d1|=%.4f  |c(D)|=%.4f  |c(Phi)|=%.4f  |c(Phi_off)|=%.4f  "
              "|c(Phi_spec)|=%.4f" % (np.linalg.norm(d1), np.linalg.norm(cD), np.linalg.norm(cPhi),
                                      np.linalg.norm(cOff), np.linalg.norm(cSpec)))
        print("    DROPPED by a commutator closure: |c(Phi_spec)|/|c(Phi)| = %.4f ; "
              "|c(Phi_spec)|/|gamma_1-dot| = %.4f"
              % (np.linalg.norm(cSpec) / np.linalg.norm(cPhi), np.linalg.norm(cSpec) / np.linalg.norm(d1)))
        print("    rank-1 spectral fraction (w.r.t. gamma_1): exact %.4f  c(D) %.4f  c(D)+c(Phi_off) %.4f  "
              "c(Phi_spec) %.4f"
              % (spectral_frac(g1, d1), spectral_frac(g1, cD), spectral_frac(g1, cD + cOff),
                 spectral_frac(g1, cSpec)))
        print()
