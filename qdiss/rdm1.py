"""
qdiss.rdm1 -- Dissipative one-electron (1-RDM) equation of motion obtained by
closing the BBGKY-like hierarchy with Mukherjee-Kutzelnigg (generalized)
normal ordering about the *instantaneous, correlated* 1-RDM.

Model
-----
    H   = H_S + H_B + H_SB
    H_S = sum_pq h_pq a^dag_p a_q
    H_B = sum_k omega_k b^dag_k b_k
    H_SB= sum_alpha A_alpha (x) B_alpha ,  A_alpha = sum_pq A^alpha_pq a^dag_p a_q

The coupling is particle-number conserving, so N = tr(gamma) is a constant of
the motion and the stationary state is Fermi-Dirac at the bath temperature with
mu fixed by N.

Key structural point
--------------------
Because A_alpha is quadratic, the TCL2 dissipator is quartic in fermion
operators and d(gamma)/dt involves the 2-RDM.  MK normal ordering about gamma
plus truncation of the 2-cumulant (lambda_2 = 0),

    <a^dag_p a^dag_q a_s a_r> = gamma_pr gamma_qs - gamma_ps gamma_qr + lambda_pqrs
                              ~= gamma_pr gamma_qs - gamma_ps gamma_qr

closes the hierarchy and produces the Pauli-blocking factors (1 - n_p), which
are literally the *hole contractions* <a_p a^dag_p> = 1 - gamma_pp of the
generalized normal ordering.  Those factors are what make the fixed point
Fermi-Dirac rather than Boltzmann, and what make 0 <= gamma <= 1 an invariant
region.
"""
from __future__ import annotations
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


# ----------------------------------------------------------------------------
# Bath: one-sided Fourier transform (real part) of the bath correlation function
# ----------------------------------------------------------------------------
def bose(w, beta):
    w = np.asarray(w, dtype=float)
    out = np.empty_like(w)
    small = np.abs(beta * w) < 1e-12
    out[~small] = 1.0 / np.expm1(beta * w[~small])
    out[small] = 1.0 / (beta * np.where(small, 1.0, 1.0)) - 0.5   # placeholder
    return out


def ohmic_J(w, eta=0.05, wc=5.0):
    """Ohmic spectral density with exponential cutoff, odd in w."""
    w = np.asarray(w, dtype=float)
    return eta * w * np.exp(-np.abs(w) / wc)


def gamma_rate(w, beta, eta=0.05, wc=5.0):
    """
    Gamma(w) = int dt e^{i w t} <B(t) B(0)> = 2 pi J(w) [ n(w) + 1 ].

    Satisfies the KMS / detailed-balance relation  Gamma(-w) = e^{-beta w} Gamma(w)
    exactly, which is the *only* property needed for the Fermi-Dirac fixed point.
    Handles w -> 0 by the classical limit 2 pi J'(0)/beta.
    """
    w = np.atleast_1d(np.asarray(w, dtype=float))
    out = np.empty_like(w)
    small = np.abs(beta * w) < 1e-10
    ws = w[~small]
    # J(w)[n(w)+1] = J(w) / (1 - e^{-beta w})   (odd J makes this positive for both signs)
    out[~small] = 2.0 * np.pi * ohmic_J(ws, eta, wc) / (1.0 - np.exp(-beta * ws))
    out[small] = 2.0 * np.pi * eta / beta            # lim_{w->0} 2pi * eta*w /(beta w)
    return out


# ----------------------------------------------------------------------------
# Fermi-Dirac reference
# ----------------------------------------------------------------------------
def fermi_dirac(eps, beta, mu):
    x = beta * (np.asarray(eps, float) - mu)
    return 0.5 * (1.0 - np.tanh(0.5 * x))          # numerically safe 1/(1+e^x)


def solve_mu(eps, beta, N):
    f = lambda mu: fermi_dirac(eps, beta, mu).sum() - N
    lo, hi = eps.min() - 200.0 / beta, eps.max() + 200.0 / beta
    return brentq(f, lo, hi, xtol=1e-14, rtol=1e-15)


# ----------------------------------------------------------------------------
# The EOM
# ----------------------------------------------------------------------------
class DissipativeRDM1:
    """
    h        : (M,M) hermitian one-body hamiltonian
    couplings: list of (M,M) hermitian matrices A^alpha (system side of H_SB)
    beta     : inverse bath temperature
    blocking : if False, drop the MK hole factors (1-n)  -> naive/linear Redfield
    """

    def __init__(self, h, couplings, beta, eta=0.05, wc=5.0, blocking=True):
        self.h = np.asarray(h, dtype=complex)
        self.M = self.h.shape[0]
        self.beta = float(beta)
        self.blocking = bool(blocking)
        self.eps, self.U = np.linalg.eigh(self.h)          # h = U diag(eps) U^dag
        # rotate couplings into the eigenbasis of h
        self.A = [self.U.conj().T @ np.asarray(A, dtype=complex) @ self.U
                  for A in couplings]
        # W[p, q] = rate for  q -> p   (population transfer)
        dE = self.eps[None, :] - self.eps[:, None]          # dE[p,q] = eps_q - eps_p
        G = gamma_rate(dE.ravel(), beta, eta, wc).reshape(dE.shape)
        self.W = np.zeros((self.M, self.M))
        for A in self.A:
            self.W += np.abs(A) ** 2 * G
        np.fill_diagonal(self.W, 0.0)

    # -- right-hand side, gamma in the eigenbasis of h ------------------------
    def rhs(self, t, y):
        M = self.M
        g = y.view(complex).reshape(M, M) if y.dtype == np.float64 else y.reshape(M, M)
        n = np.real(np.diag(g))
        hole = (1.0 - n) if self.blocking else np.ones(M)

        # coherent part: -i [eps, gamma]
        dg = -1j * (self.eps[:, None] - self.eps[None, :]) * g

        # populations: nonlinear fermionic Pauli master equation
        gain = self.W @ (n) * hole                      # sum_q W[p,q] n_q (1-n_p)
        loss = n * (self.W.T @ hole)                    # sum_q W[q,p] n_p (1-n_q)
        dpop = gain - loss

        # coherences: secular decay at the mean of the total in/out rates
        Gam = (self.W.T @ hole) + (self.W @ n)          # total rate out of / into p
        deco = -0.5 * (Gam[:, None] + Gam[None, :]) * g
        np.fill_diagonal(deco, 0.0)
        dg = dg + deco
        dg[np.diag_indices(M)] = dpop
        return dg

    def propagate(self, g0, t_end, n_out=200, rtol=1e-10, atol=1e-12):
        M = self.M
        f = lambda t, y: self.rhs(t, y.view(np.complex128)).ravel().view(np.float64)
        sol = solve_ivp(f, (0.0, t_end), np.asarray(g0, complex).ravel().view(np.float64),
                        t_eval=np.linspace(0, t_end, n_out), rtol=rtol, atol=atol,
                        method="RK45")
        gs = sol.y.T.copy().view(np.complex128).reshape(-1, M, M)
        return sol.t, gs

    # -- diagnostics ----------------------------------------------------------
    def free_energy(self, g):
        n = np.clip(np.real(np.linalg.eigvalsh(g)), 1e-15, 1 - 1e-15)
        E = np.real(np.trace(np.diag(self.eps) @ g))
        S = -np.sum(n * np.log(n) + (1 - n) * np.log(1 - n))
        return E - S / self.beta
