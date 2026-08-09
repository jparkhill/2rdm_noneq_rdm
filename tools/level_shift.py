#!/usr/bin/env python3
"""Level-shifted adiabatic closure: does regularisation buy Lipschitz continuity, and at what cost?

THE CLOSURE. The unregularised adiabatic amplitude solves L0[lambda_2] = -D, i.e.
lambda_2 = -L0^+ [D] with L0^+ the Moore-Penrose pseudo-inverse. L0 is the commutator with the
generalised Fock operator f = h + V.gamma_1 acting on the two-particle cumulant; D is the
lambda_2-independent driving. That closure FAILS the Lipschitz screen: L0's eigenvalues
i(eps_p + eps_v - eps_q - eps_w) cross zero as gamma_1 varies, so the pseudo-inverse's norm
diverges and its rank jumps -- the closure is discontinuous in gamma_1, not merely non-Lipschitz.

The level-shifted family replaces the pseudo-inverse by a Tikhonov-regularised solve, in the
Hermitian antisymmetric pair basis where L0 is a real matrix Lm and D has coefficient vector d:

    lambda_2(sigma) = - (Lm^T Lm + sigma^2 I)^{-1} Lm^T d

For every FIXED sigma > 0 this is Lipschitz in gamma_1, because Lm depends polynomially on
gamma_1 and the inverted matrix has condition number bounded by (sigma_max^2 + sigma^2)/sigma^2,
with no rank jump. sigma is part of the closure's definition, not a numerical tolerance. As
sigma -> 0 it recovers the unregularised amplitude, so the MP2 gate becomes a statement about the
limit rather than a property required at finite sigma.

THREE THINGS MEASURED, at fixed molecular Hamiltonian (jittered asymmetric-trapezoid H4/sto-3g,
M = 8 spin-orbitals, N = 4 electrons):

1. LIPSCHITZ CONSTANT, empirically. Perturb gamma_1 by a small traceless Hermitian step and
   measure |delta lambda_2| / |delta gamma_1|. Unbounded as sigma -> 0 is the expected signature;
   a finite plateau at sigma > 0 is what regularisation is supposed to buy.

2. THE PAULI TANGENCY CONDITION, at the Hartree-Fock density where occupations are exactly 0 and
   1 so the condition applies literally. With X[p,q] = sum_{r,t,u} V[r,t,p,u] lambda_2[r,t,q,u],
   the rate of change of occupation p is dn_p/dt = -Im X[p,p], so the PTC requires
       Im X[p,p] <= 0 at empty orbitals (n_p = 0),  Im X[p,p] >= 0 at full orbitals (n_p = 1).
   The exact answer is Im X[p,p] = 0 -- tangency -- so any nonzero value at finite sigma is a
   regularisation artifact and its SIGN decides whether occupations leave [0,1].

3. ACCURACY COST, at the exact FCI density: |lambda_2(sigma) - lambda_2_exact| against the exact
   cumulant, plus the energy. A plateau means Lipschitz continuity is nearly free; monotone
   degradation means we are trading one requirement against the other and should know the rate.
"""
import sys
import numpy as np
sys.path.insert(0, "."); sys.path.insert(0, "tools")
import hylleraas as H
from qdiss import closures
from qdiss.eom import fullci

M, N, h, V, enuc = H.M, H.N, H.h, H.V, H.enuc
g_hf, g_fci, E_HF, E_FCI = H.g_hf, H.g_fci, H.E_HF, H.E_FCI
NB, BASIS, coeff = H.NB, H.BASIS, H.coeff

SIGMAS = (0.0, 1e-6, 1e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0)


def L0_matrix(g1):
    Lm = np.zeros((NB, NB))
    for k, bk in enumerate(BASIS):
        Lm[:, k] = coeff(H.apply_L0(g1, bk))
    return Lm


def lam2_sigma(g1, sigma, Lm=None, d=None):
    if Lm is None:
        Lm = L0_matrix(g1)
    if d is None:
        d = coeff(H.driving(g1))
    if sigma == 0.0:
        x = np.linalg.lstsq(Lm, -d, rcond=1e-8)[0]
    else:
        A = Lm.T @ Lm + sigma ** 2 * np.eye(NB)
        x = -np.linalg.solve(A, Lm.T @ d)
    return sum(xi * b for xi, b in zip(x, BASIS))


def Xmat(l2):
    return np.einsum("rtpu,rtqu->pq", V, l2)


if __name__ == "__main__":
    print("system: jittered asymmetric-trapezoid H4/sto-3g, M=%d N=%d" % (M, N))
    print("        E_HF = %.9f   E_FCI = %.9f" % (E_HF, E_FCI))
    g2ex = fullci.rdm(H.S.psi0, M, 2)
    l2ex = g2ex - closures.wedge(g_fci, g_fci)
    print("        |lambda_2_exact| = %.5f" % np.abs(l2ex).max())
    print()

    # --- 1. empirical Lipschitz constant, measured at the FCI density -------------------------------
    rng = np.random.default_rng(4)
    K = rng.normal(size=(M, M)) + 1j * rng.normal(size=(M, M))
    K = 0.5 * (K + K.conj().T)
    K -= np.eye(M) * np.trace(K).real / M           # traceless: preserves tr gamma_1 = N
    K /= np.abs(K).max()
    EPS = 1e-5

    Lm_f, d_f = L0_matrix(g_fci), coeff(H.driving(g_fci))
    Lm_h, d_h = L0_matrix(g_hf), coeff(H.driving(g_hf))
    print("TABLE 1  empirical Lipschitz constant |d lambda_2| / |d gamma_1| at the FCI density")
    print("         (finite difference, step |d gamma_1| = %.0e traceless Hermitian)" % (EPS,))
    print("   %10s | %14s %14s" % ("sigma", "Lipschitz est", "|lambda_2|"))
    lip = {}
    for s in SIGMAS:
        a = lam2_sigma(g_fci, s, Lm_f, d_f)
        b = lam2_sigma(g_fci + EPS * K, s)
        lip[s] = np.abs(b - a).max() / (EPS * np.abs(K).max())
        print("   %10.0e | %14.3e %14.5f" % (s, lip[s], np.abs(a).max()))
    print()

    # --- 2. Pauli Tangency Condition, at the HF density where occupations are exactly 0 and 1 --------
    n_hf = np.real(np.diag(g_hf))
    occ_idx = [p for p in range(M) if n_hf[p] > 0.5]
    vir_idx = [p for p in range(M) if n_hf[p] < 0.5]
    print("TABLE 2  Pauli Tangency Condition at gamma_HF (occupations exactly 1 on %s, 0 on %s)"
          % (occ_idx, vir_idx))
    print("         dn_p/dt = -Im X[p,p];  PTC needs Im X <= 0 on empty, >= 0 on full; exact = 0")
    print("   %10s | %14s %14s | %s" % ("sigma", "max Im X (empty)", "min Im X (full)", "verdict"))
    for s in SIGMAS:
        l2 = lam2_sigma(g_hf, s, Lm_h, d_h)
        X = Xmat(l2)
        im = np.imag(np.diag(X))
        hi, lo = im[vir_idx].max(), im[occ_idx].min()
        ok = (hi <= 1e-10) and (lo >= -1e-10)
        print("   %10.0e | %+14.3e %+14.3e | %s"
              % (s, hi, lo, "PASS" if ok else "VIOLATED"))
    print()

    # --- 3. accuracy cost, at the FCI density -------------------------------------------------------
    A = g_fci - g_fci @ g_fci
    PAIRS = H.PAIRS


    def mineig2(T):
        Aq = np.array([[2.0 * T[p, v, q, w] for (q, w) in PAIRS] for (p, v) in PAIRS])
        return np.linalg.eigvalsh(0.5 * (Aq + Aq.conj().T)).min()


    print("TABLE 3  accuracy cost of the level shift, at the exact FCI density")
    print("   %10s | %11s %11s %12s %11s" % ("sigma", "|l2-l2ex|", "trace-id", "E-E_FCI", "min eig D"))
    for s in SIGMAS:
        l2 = lam2_sigma(g_fci, s, Lm_f, d_f)
        g2 = closures.wedge(g_fci, g_fci) + l2
        print("   %10.0e | %11.5f %11.5f %+12.6f %11.6f"
              % (s, np.abs(l2 - l2ex).max(),
                 np.abs(np.einsum("pxqx->pq", l2) + A).max(),
                 float(np.real(np.einsum("pq,pq->", h, g_fci)
                               + 0.25 * np.einsum("pvqw,pvqw->", V, g2))) + enuc - E_FCI,
                 mineig2(g2)))
