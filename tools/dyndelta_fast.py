#!/usr/bin/env python3
r"""Episode 10b: the O(M^6) formulation of the dynamical-delta force (docs/DYNDELTA-SPEC.md section 6), torch/CUDA.

The forces on delta are  g_k = sum_P <A_P, Pdot(B_k)>_F  with A_P a state-dependent Hermitian matrix on cone P
and Pdot the FIXED linear map from a pair-space flow F to the rate of cone P. Since {B_k} is an orthonormal basis of
Hermitian pair matrices, the vector g is the coordinate vector of  sum_P Pdot^dag(A_P)  -- the ADJOINT maps applied
once. No basis images are stored or rotated. The hard-constraint projector Z Z^T = I - C^+ C is applied implicitly.

The adjoints are obtained here numerically-exactly by transposing the linear maps ONCE in a small basis? No: they
are written out explicitly below, and VERIFIED against the stored-image implementation of tools/dyndelta_gpu.py at
round-off on a random state before use. Delta is then represented directly as the Hermitian matrix, not as
coefficients (the coefficient basis is only needed for the constant constraint matrix C, which is (M^2+1) x n2 and
applied as a linear functional on Hermitian matrices: rows = contraction to gamma_1 and energy).

Cost per RHS: F0 O(M^6); three eigh O(M^6); adjoints O(M^4); constraints O(M^6). Memory O(M^4).
"""
import sys, os, time
import numpy as np
import torch
sys.path.insert(0, "."); sys.path.insert(0, "tools")
import commutator_structure as cs
from f0delta import F0Delta
from dyndelta_gpu import GPUFlow, up, dev, C

F64 = torch.float64


class FastFlow(GPUFlow):
    """Same equations as GPUFlow; forces via adjoint maps; delta stored as a Hermitian pair matrix."""

    def __init__(self, base, tau, kappa, sigma, kappa_neg=0.0):
        super().__init__(base, tau, kappa, sigma, kappa_neg)
        # release the big precomputed images (not used here)
        self.QB = self.GB = self.DB = None
        M, N, nP = self.M, self.N, self.nP
        # hard constraints as linear functionals on Hermitian pair matrices F:
        #   contraction:  hvec( Tr_2 from_mat2(F) / (N-1) )     (M^2 real numbers)
        #   energy:       Re tr(gradE^dag F)                       (1 real number)
        # Represent the constraint operator by its ADJOINT action too: for the projector we need
        #   P_perp(D) = D - C^+ (C D),  C^+ = C^dag (C C^dag)^{-1}  (C has full row rank M^2+1 on Hermitian matrices)
        # Build C C^dag (small, (M^2+1)^2) by applying C to the adjoint images of unit constraint vectors.
        self.nc = M * M + 1
        E = torch.eye(self.nc, dtype=F64, device=dev)
        adj = torch.stack([self._C_adj(E[i]) for i in range(self.nc)])          # (nc, nP, nP) Hermitian
        CCd = torch.stack([self._C(adj[i]) for i in range(self.nc)])           # (nc, nc)
        self.CCd_inv = torch.linalg.pinv(CCd, rtol=1e-12)
        self._adj = adj

    # ---- constraint functional C: Hermitian pair matrix -> R^{M^2+1}
    def _C(self, F):
        g2F = self.from_mat2(F); g1F = self.g1_of(g2F)
        return torch.cat([self.hvec(g1F), torch.real(torch.trace(self.gradE.conj().T @ F))[None]])

    def _hvec_adj(self, v):
        """adjoint of hvec: real vector -> Hermitian M x M."""
        M = self.M; iu = self.iu
        A = torch.zeros((M, M), dtype=C, device=dev)
        A[torch.arange(M), torch.arange(M)] = v[:M].to(C)
        nu = iu.shape[1]
        re, im = v[M:M + nu], v[M + nu:]
        A[iu[0], iu[1]] = (re + 1j * im).to(C) / np.sqrt(2)
        A[iu[1], iu[0]] = (re - 1j * im).to(C) / np.sqrt(2)
        return A

    def _C_adj(self, v):
        """adjoint of C: R^{M^2+1} -> Hermitian pair matrix."""
        M, N = self.M, self.N
        A1 = self._hvec_adj(v[:-1])                    # M x M Hermitian
        # adjoint of F -> Tr_2 from_mat2(F)/(N-1): the pair lift of A1 times the right factor.
        # <A1, Tr_2 from_mat2(F)>/(N-1) = <T, F> with T = to_mat2( A1 (x) I  antisymmetrised ) / (2 (N-1)) ...
        # do it exactly via the known lift: Tr_2 from_mat2(F)[p,q] = sum_v from_mat2(F)[p,v,q,v]; from_mat2 spreads
        # F/2 over four antisymmetric slots, so the adjoint is 4 * (1/2) * antisym( A1[p,q] d_vw ) restricted to pairs:
        Tfull = torch.einsum("pq,vw->pvqw", A1, self.I2)
        Tfull = self.antisym4(Tfull) * 4.0 / 2.0
        Tpair = self.to_mat2(Tfull) / 2.0 / (N - 1)   # to_mat2 has the factor 2; undo it
        return Tpair + v[-1] * self.gradE

    def P_perp(self, D):
        """Project a Hermitian pair matrix onto the null space of the hard constraints."""
        y = self.CCd_inv @ self._C(D)
        return D - torch.einsum("i,iab->ab", y.to(C), self._adj)

    def C_particular(self, r):
        """Minimum-norm Hermitian pair matrix with C(.) = r."""
        y = self.CCd_inv @ r
        return torch.einsum("i,iab->ab", y.to(C), self._adj)

    # ---- adjoints of the cone-rate maps (Hermitian on cone -> Hermitian pair matrix)
    def Ddot_adj(self, A):
        return A / 2.0

    def _g1_adj(self, B1):
        """adjoint of F -> g1F = Tr_2 from_mat2(F)/(N-1), as a pair matrix (same as the contraction adjoint)."""
        Tfull = self.antisym4(torch.einsum("pq,vw->pvqw", B1, self.I2)) * 2.0
        return self.to_mat2(Tfull) / 2.0 / (self.N - 1)

    def _g2_adj(self, T4):
        """adjoint of F -> from_mat2(F) (4-index) : T4 -> pair matrix = to_mat2(antisym T4)/2 * 4/2 ... derive:
        <T4, from_mat2(F)> = sum over the 4 antisymmetric copies of F/2 -> = 2 <antisym(T4) restricted, F>."""
        return self.to_mat2(self.antisym4(T4)) * 2.0 / 2.0

    def Qdot_adj(self, A):
        """Qdot(F) = to_mat2(T)/2 with T = -d_vw g1F[q,p] + d_vq g1F[w,p] + d_pw g1F[q,v] - d_pq g1F[w,v] + g2F[w,q,v,p].
        Adjoint: A (Hermitian on pairs) -> pair matrix."""
        M = self.M; I = self.I2
        # step 1: adjoint of to_mat2(.)/2 : A -> 4-index tensor (A placed on sorted pairs, antisymmetric copies)
        TA = self.from_mat2(A)          # from_mat2 spreads A/2 antisymmetrically: exactly the adjoint of to_mat2(.)/2 up to factor
        # <A, to_mat2(T)/2> = <A, T[pairs]> = <from_mat2(A)*?, T> ... to_mat2 picks entries (2x); /2 cancels: <A,T_sorted>.
        # <A, T_sorted> = sum_{p<v,q<w} conj(A) T = (1/4) sum_all conj(antisym-extended A) T  -> adjoint tensor = from_mat2(A)/... :
        # from_mat2(A)[p,v,q,w] = A/2 on sorted, antisymmetric elsewhere; sum over all 4 copies gives 4*(A/2)/... use explicit:
        TA = 2.0 * TA                   # now TA = A on sorted pairs, antisymmetric copies; <TA, T>_all = 4 <A, T_sorted>
        TA = TA / 4.0
        # step 2: adjoints of the five terms w.r.t. (g1F, g2F)
        B1 = (- torch.einsum("pvqw,vw->qp", TA, I)     # from -d_vw g1F[q,p]
              + torch.einsum("pvqw,vq->wp", TA, I)     # + d_vq g1F[w,p]
              + torch.einsum("pvqw,pw->qv", TA, I)     # + d_pw g1F[q,v]
              - torch.einsum("pvqw,pq->wv", TA, I))    # - d_pq g1F[w,v]
        # careful: adjoint w.r.t. g1F[a,b] means coefficient of g1F[a,b]; the einsum outputs above index (a,b) directly
        T2 = TA.permute(3, 2, 1, 0)                     # coefficient of g2F[w,q,v,p] -> tensor indexed [w,q,v,p]
        return self._g1_adj(B1.conj().T.conj()) + self._g2_adj(T2)

    def Gdot_adj(self, A):
        """Gdot(F) = (d_qs g1F[p,r] - g2F[p,s,r,q]).reshape(M^2,M^2); adjoint of A (M^2 x M^2 Hermitian)."""
        M = self.M
        TA = A.reshape(M, M, M, M)                      # indexed [p,q,r,s]
        B1 = torch.einsum("pqrs,qs->pr", TA, self.I2)  # coefficient of g1F[p,r]
        T2 = -TA.permute(0, 3, 2, 1)                    # coefficient of g2F[p,s,r,q] -> tensor indexed [p,s,r,q]... permute (p,q,r,s)->(p,s,r,q)
        return self._g1_adj(B1) + self._g2_adj(T2)

    # ---- the RHS
    def rhs(self, g2, delta):
        g1 = self.g1_of(g2)
        F0 = self.F0(g1, g2)
        F = F0 + delta
        D, Q, G = self.cones(g1, g2)
        Dd, Qd, Gd, g1F = self.cone_rates(F)
        force = torch.zeros_like(F); phi = 0.0
        for Pm, Pd, adj in ((D, Dd, self.Ddot_adj), (Q, Qd, self.Qdot_adj), (G, Gd, self.Gdot_adj)):
            x, U = torch.linalg.eigh(0.5 * (Pm + Pm.conj().T))
            sf = torch.sqrt(self.sigma ** 2 / (x ** 2 + self.sigma ** 2))
            Wm = U * sf[None, :].to(C)
            RF = Wm.conj().T @ (0.5 * (Pd + Pd.conj().T)) @ Wm
            phi += float(torch.linalg.norm(RF) ** 2)
            A_P = 2.0 * self.kappa * (Wm @ RF @ Wm.conj().T)             # d/dPdot of kappa ||W^dag Pdot W||^2
            if self.kappa_neg > 0:
                eps = self.eps_bar
                sp = torch.nn.functional.softplus(-x / eps); bp = -2 * eps * sp * torch.sigmoid(-x / eps)
                A_P = A_P + self.kappa_neg * (U * bp[None, :].to(C)) @ U.conj().T
                phi += self.kappa_neg / max(self.kappa, 1e-300) * float((eps ** 2 * sp ** 2).sum())
            force = force + adj(0.5 * (A_P + A_P.conj().T))
        # hard constraints: particular solution for F0, then project the driving force
        g2F0 = self.from_mat2(F0); g1F0 = self.g1_of(g2F0)
        r_hard = torch.cat([self.hvec(self.exact_g1dot(g1, g2) - g1F0),
                            (-torch.real(torch.trace(self.gradE.conj().T @ F0)))[None]])
        delta_hard = self.C_particular(r_hard)
        ddelta = -(delta - delta_hard) - self.P_perp(force)
        self.last_phi = phi
        return self.from_mat2(F), ddelta / self.tau


def validate():
    """Compare FastFlow.rhs against GPUFlow.rhs (stored images) at a random correlated state and random delta."""
    atom, N, spin = cs.H4
    label, dm, dets, M, h, V = next(cs.states(atom, N, spin, "H4"))
    from qdiss.eom import fullci
    H = fullci.build_H(h, V, M); (g1, g2, g3, d1, d2), _ = cs.ensemble_rdms(dm, dets, M, H)
    base = F0Delta(h, V, M, N)
    for kn in (0.0, 1e4):
        slow = GPUFlow(base, 0.5, 100.0, 1e-3, kn); fast = FastFlow(base, 0.5, 100.0, 1e-3, kn)
        rng = np.random.default_rng(1)
        dc = torch.as_tensor(rng.normal(size=base.n2) * 1e-2, device=dev)
        delta = torch.einsum("k,kij->ij", dc.to(C), slow.B)
        # push the state slightly outside the cone for the barrier test
        g2p = g2 - 0.02 * (g2 - np.einsum("pq,vw->pvqw", g1, g1) + np.einsum("pw,vq->pvqw", g1, g1)) if kn > 0 else g2
        g2t = up(g2p)
        gd_s, dcd_s = slow.rhs(g2t, dc)
        gd_f, dd_f = fast.rhs(g2t, delta)
        # convert fast's matrix rate to coefficients for comparison
        dd_f_c = torch.real(torch.einsum("kij,ij->k", slow.B.conj(), dd_f))
        print("kappa_neg=%.0e : |g2dot diff| %.1e   |delta-dot diff| %.1e  (|delta-dot| = %.3f)   phi slow %.3e fast %.3e"
              % (kn, float((gd_s - gd_f).abs().max()), float((dcd_s - dd_f_c).abs().max()), float(dcd_s.norm()), slow.last_phi, fast.last_phi))
        t0 = time.time()
        for _ in range(10): fast.rhs(g2t, delta)
        torch.cuda.synchronize(); tf = (time.time() - t0) / 10
        t0 = time.time()
        for _ in range(10): slow.rhs(g2t, dc)
        torch.cuda.synchronize(); ts = (time.time() - t0) / 10
        print("    time per rhs: stored-images %.3fs   adjoint %.3fs" % (ts, tf))


if __name__ == "__main__":
    validate()
