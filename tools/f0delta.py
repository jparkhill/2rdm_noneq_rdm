#!/usr/bin/env python3
r"""Route F, episode 7: the closure  d/dt X = F0 + delta.

Symbols (docs/COMMUTATOR-STRUCTURE.md). X = Mat_2 (pair space, sorted pairs, X = 2 gamma_2), gamma_1 = Tr_2 gamma_2/(N-1).
F0(gamma_1, gamma_2): the lambda_3 = 0 flow in closed channel form (tools/closed_form_sep.py, verified 1e-15).
P in {D, Q, G}: the pp, hh, ph Gram matrices, affine in (gamma_1, gamma_2) (tools/channel_cones.py), eigenpairs (x_i, u_i).

delta is the Hermitian pair-space matrix minimising

    ||delta||_F^2  +  mu * sum_P sum_ij  w_ij^2 | u_i^dag dP/dt[F0 + delta] u_j |^2 ,     w_ij = sigma / (x_i + x_j + sigma)

subject to the HARD linear constraints
    (b)  Tr_2 (F0 + delta) / (N-1)  =  exact rank-1 EOM at (gamma_1, gamma_2)     (the rank-1 EOM has no gamma_3)
    (c)  Re tr( grad E . (F0 + delta) ) = 0,   grad E = K_h/(2(N-1)) + K_V/2       (energy conservation)

(a)_sigma is the smoothed tangency: w_ij -> 1 in near-null blocks (x_i + x_j << sigma), where the exact flow's
entries vanish, and w_ij -> 0 in the interior. The exact flow satisfies (a) up to O(sigma) (its entries scale as
x_i + x_j), (b) and (c) exactly, so it is a fixed point of the correction up to O(sigma).

Everything linear in delta is precomputed once: a real basis {B_k} of Hermitian pair matrices and, for each B_k,
its induced gamma_2-dot, gamma_1-dot, Qdot, Gdot (all constant linear images). Per evaluation: three
eigendecompositions, rotation of the precomputed images into the cone eigenbases, one least-squares solve in the
null space of the hard constraints.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
import commutator_structure as cs
from commutator_k3 import from_mat2, pair_k2
from closed_form_sep import closed_form
from channel_cones import q_matrix, g_matrix
from qdiss import ansatz as az
from qdiss import muller2405
from qdiss.eom import wick


class F0Delta:
    def __init__(self, h, V, M, N, sigma=1e-3, mu=1e2, cones=("D", "Q", "G"), tau=0.0):
        self.h, self.V, self.M, self.N = h, V, M, N
        self.sigma, self.mu, self.cones, self.tau = sigma, mu, cones, tau
        self.P = cs.pairs(M)
        nP = len(self.P)
        self.nP = nP
        self.EYE = np.eye(M, dtype=complex)
        # real basis of Hermitian nP x nP matrices
        B = []
        for i in range(nP):
            E = np.zeros((nP, nP), dtype=complex); E[i, i] = 1.0; B.append(E)
        for i in range(nP):
            for j in range(i + 1, nP):
                E = np.zeros((nP, nP), dtype=complex); E[i, j] = E[j, i] = 1 / np.sqrt(2); B.append(E)
                E = np.zeros((nP, nP), dtype=complex); E[i, j] = 1j / np.sqrt(2); E[j, i] = -1j / np.sqrt(2); B.append(E)
        self.B = np.array(B)                                   # (n2, nP, nP), orthonormal in Frobenius
        n2 = len(B)
        self.n2 = n2
        # constant linear images of each basis element
        g2d = np.array([from_mat2(Bk, M) for Bk in B])         # (n2, M,M,M,M)
        self.g1d = np.einsum("kpvqv->kpq", g2d) / (N - 1)       # (n2, M, M)
        Q0 = q_matrix(0 * self.EYE, np.zeros((M,) * 4, dtype=complex), M)
        self.QB = np.array([q_matrix(self.g1d[k], g2d[k], M) - Q0 for k in range(n2)])   # (n2, nP, nP)
        self.GB = np.array([g_matrix(self.g1d[k], g2d[k], M) for k in range(n2)])       # (n2, M^2, M^2)
        self.DB = self.B / 2.0
        # hard constraints: contraction (M^2 real: Hermitian g1dot -> real vec) and energy (1 real)
        gradE = pair_k2(h, 0 * V, M) / (2 * (N - 1)) + pair_k2(0 * h, V, M) / 2
        self.gradE = gradE
        Cc = np.array([self._hvec(self.g1d[k]) for k in range(n2)]).T                   # (M^2, n2) real
        Ce = np.array([np.real(np.trace(gradE.conj().T @ Bk)) for Bk in B])[None, :]    # (1, n2)
        self.C = np.vstack([Cc, Ce])
        # null space and pseudo-inverse of the constant constraint matrix
        U, s, Vt = np.linalg.svd(self.C, full_matrices=True)
        r = int((s > 1e-10 * s[0]).sum())
        self.Cpinv = Vt[:r].T @ np.diag(1 / s[:r]) @ U[:, :r].T
        self.Z = Vt[r:].T                                       # (n2, n2 - r)
        self.rankC = r
        self.last = {}

    @staticmethod
    def _hvec(A):
        """Real vector of a Hermitian matrix (diagonal, then Re/Im of the upper triangle)."""
        n = A.shape[0]
        iu = np.triu_indices(n, 1)
        return np.concatenate([np.real(np.diag(A)), np.sqrt(2) * np.real(A[iu]), np.sqrt(2) * np.imag(A[iu])])

    # ------------------------------------------------------------------ pieces
    def state(self, g2):
        g1 = np.einsum("pvqv->pq", g2) / (self.N - 1)
        return g1, cs.to_mat2(g2, self.M)

    def exact_g1dot(self, g1, g2):
        tab = {("D", 1): self.EYE, ("G", 1): g1, ("G", 2): g2}
        return wick.gamma_dot(1, self.h, self.V, tab, self.M)

    def F0(self, g1, g2):
        drift, sandwich, qterm, Phi_ph, _ = closed_form(self.h, self.V, g1, g2, self.M)
        return drift + sandwich + qterm + Phi_ph

    def muller(self, g1, g2):
        tab = {("D", 1): self.EYE, ("G", 1): g1, ("G", 2): g2,
               ("G", 3): muller2405.gamma3(g1, g2, self.N, True)}
        return cs.to_mat2(wick.gamma_dot(2, self.h, self.V, tab, self.M), self.M)

    # ------------------------------------------------------------------ the correction
    def delta(self, g1, g2, F, cones=None):
        """Minimum-norm delta for a given base flow F (pair-space Hermitian rate)."""
        M, N, nP = self.M, self.N, self.nP
        cones = self.cones if cones is None else cones
        X = cs.to_mat2(g2, M)
        Dm = X / 2
        Qm = q_matrix(g1, g2, M)
        Gm = g_matrix(g1, g2, M)
        g2F = from_mat2(F, M)
        g1F = np.einsum("pvqv->pq", g2F) / (N - 1)
        Q0 = q_matrix(0 * g1, 0 * g2, M)
        rateF = {"D": F / 2, "Q": q_matrix(g1F, g2F, M) - Q0, "G": g_matrix(g1F, g2F, M)}
        rateB = {"D": self.DB, "Q": self.QB, "G": self.GB}
        mats = {"D": Dm, "Q": Qm, "G": Gm}
        rows, rhs = [], []
        hard_rows, hard_rhs = [], []
        info = {}
        for c in cones:
            Pm = 0.5 * (mats[c] + mats[c].conj().T)
            x, U = np.linalg.eigh(Pm)
            w = self.sigma / (np.abs(x[:, None]) + np.abs(x[None, :]) + self.sigma)
            RF = U.conj().T @ (0.5 * (rateF[c] + rateF[c].conj().T)) @ U
            RB = np.matmul(np.matmul(U.conj().T[None], rateB[c]), U)
            hard = (np.abs(x)[:, None] < self.tau) & (np.abs(x)[None, :] < self.tau)   # near-null block: exact tangency
            keep = (w > 1e-3) & ~hard if self.mu > 0 else np.zeros_like(hard)   # mu = 0: hard constraints only
            wk = w[keep]
            rows.append(np.sqrt(self.mu) * (RB[:, keep] * wk[None, :]).T)      # (n_keep, n2) complex
            rhs.append(-np.sqrt(self.mu) * (RF[keep] * wk))
            if hard.any():
                hard_rows.append(RB[:, hard].T)
                hard_rhs.append(-RF[hard])
            info["nhard_" + c] = int(hard.sum())
            info["min_eig_" + c] = x.min()
            info["nkeep_" + c] = int(keep.sum())
        A = np.vstack(rows)
        b = np.concatenate(rhs)
        A = np.vstack([A.real, A.imag])
        b = np.concatenate([b.real, b.imag])
        # hard constraints: C c = r_hard
        d1x = self.exact_g1dot(g1, g2)
        r_hard = np.concatenate([self._hvec(d1x - g1F), [-np.real(np.trace(self.gradE.conj().T @ F))]])
        if hard_rows:
            Ch = np.vstack(hard_rows); rh = np.concatenate(hard_rhs)
            Cfull = np.vstack([self.C, Ch.real, Ch.imag])
            rfull = np.concatenate([r_hard, rh.real, rh.imag])
            Uc, sc, Vct = np.linalg.svd(Cfull, full_matrices=True)
            rk = int((sc > 1e-10 * sc[0]).sum())
            Cpinv = Vct[:rk].T @ np.diag(1 / sc[:rk]) @ Uc[:, :rk].T
            Z = Vct[rk:].T
            c_p = Cpinv @ rfull
            info["rank_hard"] = rk
        else:
            Cpinv, Z = self.Cpinv, self.Z
            c_p = Cpinv @ r_hard
        # minimise ||c||^2 + ||A c - b||^2 over c = c_p + Z y
        AZ = A @ Z
        Ap = np.vstack([AZ, Z])                      # ||A(c_p+Zy) - b||^2 + ||c_p + Z y||^2
        bp = np.concatenate([b - A @ c_p, -c_p])
        y, *_ = np.linalg.lstsq(Ap, bp, rcond=None)
        c = c_p + Z @ y
        delta = np.einsum("k,kij->ij", c, self.B)
        self.last = info
        return delta

    def flow(self, g2, mode="f0delta"):
        g1, X = self.state(g2)
        if mode == "f0":
            return self.F0(g1, g2)
        if mode == "muller":
            return self.muller(g1, g2)
        if mode == "f0delta":
            F = self.F0(g1, g2)
            return F + self.delta(g1, g2, F)
        if mode == "f0delta_D":
            F = self.F0(g1, g2)
            return F + self.delta(g1, g2, F, cones=("D",))
        if mode == "f0delta_hard":
            F = self.F0(g1, g2)
            tau0, self.tau = self.tau, 1e-5
            try:
                return F + self.delta(g1, g2, F)
            finally:
                self.tau = tau0
        if mode == "f0delta_hardonly":
            F = self.F0(g1, g2)
            tau0, mu0, self.tau, self.mu = self.tau, self.mu, 1e-5, 0.0
            try:
                return F + self.delta(g1, g2, F)
            finally:
                self.tau, self.mu = tau0, mu0
        if mode == "muller_delta":
            F = self.muller(g1, g2)
            return F + self.delta(g1, g2, F)
        raise ValueError(mode)

    def g2dot(self, g2, mode="f0delta"):
        return from_mat2(self.flow(g2, mode), self.M)
