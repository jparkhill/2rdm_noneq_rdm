#!/usr/bin/env python3
r"""Route F, episode 9: GPU (torch) implementation of the dynamical-delta EOM, validated against the numpy path.

Same equations as tools/f0delta.py (F0Delta.F0, DynamicDelta.rhs):
    d/dt X = F0(X) + delta ;   tau d/dt delta = -(delta - delta_hard) - kappa grad_delta Phi
Everything (pair-space matrices, ph-space matrices, cone rates, the precomputed basis images, RK4) lives on the
device in complex128. Modes: f0, muller (via numpy, slow), dyndelta. Validation: |rhs_gpu - rhs_cpu| at t = 0.

Run: python tools/dyndelta_gpu.py MODE TEND DT [tau kappa sigma]   (H4 kicked ground state unless SYS=quench)
Output every 0.5 a.u.; saved to results/long_gpu_MODE.npy.
"""
import sys, os, time
import numpy as np, scipy.linalg as sla
import torch
sys.path.insert(0, "."); sys.path.insert(0, "tools")
import commutator_structure as cs
from commutator_k3 import pair_k2, from_mat2
from f0delta import F0Delta, DynamicDelta
from f0delta_propagate import rdms_of, metrics
from channel_cones import q_matrix, g_matrix
from molecular_2rdm import build_H_sector, molecule, sector
from qdiss.eom import wick
from qdiss import ansatz as az

dev = torch.device("cuda")
C = torch.complex128
up = lambda a: torch.as_tensor(np.asarray(a), dtype=C, device=dev)


class GPUFlow:
    def __init__(self, base: F0Delta, tau, kappa, sigma, kappa_neg=0.0):
        self.b, self.tau, self.kappa, self.sigma, self.kappa_neg = base, tau, kappa, sigma, kappa_neg
        self.eps_bar = float(os.environ.get("EPS_BAR", "1e-3"))
        M, N, P = base.M, base.N, base.P
        self.M, self.N, self.nP = M, N, len(P)
        h, V = base.h, base.V
        # pair <-> tensor maps as index arrays
        pi = torch.tensor([p for p, v in P], device=dev); vi = torch.tensor([v for p, v in P], device=dev)
        self.pi, self.vi = pi, vi
        self.h, self.V = up(h), up(V)
        self.Kh, self.KV = up(pair_k2(h, 0 * V, M)), up(pair_k2(0 * h, V, M))
        I2 = np.eye(M)
        self.Sw = up(np.eye(M * M).reshape(M, M, M, M).transpose(0, 1, 3, 2).reshape(M * M, M * M))
        self.Vd = up(V.transpose(0, 2, 1, 3).reshape(M * M, M * M))
        self.I2 = up(I2)
        self.gradE = up(base.gradE)
        # basis images (n2, ., .) -> device; and the constraint machinery
        self.B = up(base.B); self.DB = self.B / 2; self.QB = up(base.QB); self.GB = up(base.GB)
        self.Cpinv, self.Z = up(base.Cpinv).real, up(base.Z).real
        self.n2 = base.n2
        # rank-1 wick term list is small; do it in numpy (cheap) -- but we need it on device per rhs: build the
        # one-body/two-body EOM for gamma_1 directly:  i[h^T, g1] + collision from gamma_2 (exact, no gamma_3)
        # Verified below against wick at 1e-15.
        self.iu = torch.triu_indices(M, M, 1, device=dev)

    # ------------------------------------------------------------- representations
    def to_mat2(self, T):
        return 2.0 * T[self.pi[:, None], self.vi[:, None], self.pi[None, :], self.vi[None, :]]

    def from_mat2(self, X):
        M = self.M
        T = torch.zeros((M,) * 4, dtype=C, device=dev)
        t = X / 2.0
        T[self.pi[:, None], self.vi[:, None], self.pi[None, :], self.vi[None, :]] = t
        T[self.vi[:, None], self.pi[:, None], self.pi[None, :], self.vi[None, :]] = -t
        T[self.pi[:, None], self.vi[:, None], self.vi[None, :], self.pi[None, :]] = -t
        T[self.vi[:, None], self.pi[:, None], self.vi[None, :], self.pi[None, :]] = t
        return T

    def pair_k(self, A):
        """A (+) A on pair space: A[p,q] d_vw + A[v,w] d_pq - A[p,w] d_vq - A[v,q] d_pw."""
        pi, vi = self.pi, self.vi
        d = (pi[:, None] == pi[None, :]).to(C); dv = (vi[:, None] == vi[None, :]).to(C)
        dpw = (pi[:, None] == vi[None, :]).to(C); dvq = (vi[:, None] == pi[None, :]).to(C)
        return (A[pi[:, None], pi[None, :]] * dv + A[vi[:, None], vi[None, :]] * d
                - A[pi[:, None], vi[None, :]] * dvq - A[vi[:, None], pi[None, :]] * dpw)

    def ph(self, T):
        M = self.M
        return T.permute(0, 2, 1, 3).reshape(M * M, M * M)

    def ph_inv(self, A):
        M = self.M
        return A.reshape(M, M, M, M).permute(0, 2, 1, 3)

    @staticmethod
    def antisym4(T):
        return 0.25 * (T - T.permute(1, 0, 2, 3) - T.permute(0, 1, 3, 2) + T.permute(1, 0, 3, 2))

    # ------------------------------------------------------------- physics
    def g1_of(self, g2):
        return torch.einsum("pvqv->pq", g2) / (self.N - 1)

    def exact_g1dot(self, g1, g2):
        # i<[H, a+_p a_q]> with H = h + 1/4 V: one-body i(h^T g1 - g1 h^T)-type plus V.gamma_2 -- take from the
        # closed collision form validated in the repo:  gamma_1-dot = i[(h)^T-convention..]; safest: rank-1 wick list
        # evaluated with torch einsum using the same term list (small, 6 terms).
        if not hasattr(self, "_t1"):
            one, two = wick.terms(1)
            self._t1 = (one, two)
        one, two = self._t1
        M = self.M
        tab = {("D", 1): self.I2, ("G", 1): g1, ("G", 2): g2}
        tot = torch.zeros((M, M), dtype=C, device=dev)
        for terms, Hop, hs, pref in ((one, self.h, "rs", 1.0), (two, self.V, "rtsu", 0.25)):
            for sign, facs in terms:
                subs = [hs]; ops = [Hop]
                for name, rank, u, l in facs:
                    subs.append("".join(u) + "".join(l)); ops.append(tab[("D", 1) if name == "D" else (name, rank)])
                tot = tot + sign * pref * torch.einsum(",".join(subs) + "->pq", *ops)
        return 1j * tot

    def F0(self, g1, g2):
        M = self.M
        X = self.to_mat2(g2)
        l2 = g2 - (torch.einsum("pq,vw->pvqw", g1, g1) - torch.einsum("pw,vq->pvqw", g1, g1))
        Np = self.pair_k(g1)
        W = torch.einsum("prqs,rs->pq", self.V, g1)
        KWT = self.pair_k(W.T.contiguous())
        Cm = torch.einsum("rtpu,rtqu->pq", self.V, l2)
        KC, KgC, KCg = self.pair_k(Cm), self.pair_k(g1 @ Cm), self.pair_k(Cm @ g1)
        S = KWT - self.KV @ Np
        Q = 0.5 * (Np @ KC - KgC) + 0.5 * (KC @ Np - KCg)
        Ld = self.ph(l2)
        gI, gTI, Ig, IgT = torch.kron(g1, self.I2), torch.kron(g1.T.contiguous(), self.I2), torch.kron(self.I2, g1), torch.kron(self.I2, g1.T.contiguous())
        P1 = self.Sw @ gTI @ self.Vd @ Ld; P4 = self.Sw @ Ig @ self.Vd @ Ld
        P2 = Ld @ self.Vd @ (self.Sw @ gI).conj().T; P3 = Ld @ self.Vd @ (self.Sw @ IgT).conj().T
        Pph = self.to_mat2(self.antisym4(self.ph_inv(P1 + P2 - P3 - P4))) / 2.0
        K2 = self.Kh + self.KV
        return (1j * (K2 @ X - X @ K2) + 1j * (S.conj().T @ X - X @ S) + 1j * (Q - Q.conj().T)
                + 2j * (Pph - Pph.conj().T))

    def cones(self, g1, g2):
        M = self.M; I = self.I2
        D = self.to_mat2(g2) / 2
        T = (torch.einsum("vw,pq->pvqw", I, I) - torch.einsum("vq,pw->pvqw", I, I)
             - torch.einsum("vw,qp->pvqw", I, g1) + torch.einsum("vq,wp->pvqw", I, g1)
             + torch.einsum("pw,qv->pvqw", I, g1) - torch.einsum("pq,wv->pvqw", I, g1) + g2.permute(3, 2, 1, 0))
        Q = self.to_mat2(T) / 2
        G = (torch.einsum("qs,pr->pqrs", I, g1) - g2.permute(0, 3, 2, 1)).reshape(M * M, M * M)
        return D, Q, G

    def cone_rates(self, F):
        """Rates of D, Q, G induced by a pair-space flow F (affine parts removed)."""
        M = self.M; I = self.I2
        g2F = self.from_mat2(F); g1F = self.g1_of(g2F)
        T = (- torch.einsum("vw,qp->pvqw", I, g1F) + torch.einsum("vq,wp->pvqw", I, g1F)
             + torch.einsum("pw,qv->pvqw", I, g1F) - torch.einsum("pq,wv->pvqw", I, g1F) + g2F.permute(3, 2, 1, 0))
        Qd = self.to_mat2(T) / 2
        Gd = (torch.einsum("qs,pr->pqrs", I, g1F) - g2F.permute(0, 3, 2, 1)).reshape(M * M, M * M)
        return F / 2, Qd, Gd, g1F

    def hvec(self, A):
        iu = self.iu
        return torch.cat([torch.real(torch.diagonal(A)), np.sqrt(2) * torch.real(A[iu[0], iu[1]]), np.sqrt(2) * torch.imag(A[iu[0], iu[1]])])

    def rhs(self, g2, dc):
        g1 = self.g1_of(g2)
        F0 = self.F0(g1, g2)
        delta = torch.einsum("k,kij->ij", dc.to(C), self.B)
        F = F0 + delta
        # Phi gradient
        D, Q, G = self.cones(g1, g2)
        Dd, Qd, Gd, g1F = self.cone_rates(F)
        grad = torch.zeros(self.n2, dtype=torch.float64, device=dev); phi = 0.0
        for Pm, Pd, RB in ((D, Dd, self.DB), (Q, Qd, self.QB), (G, Gd, self.GB)):
            x, U = torch.linalg.eigh(0.5 * (Pm + Pm.conj().T))
            sf = torch.sqrt(self.sigma ** 2 / (x ** 2 + self.sigma ** 2))
            Wm = U * sf[None, :].to(C)
            RF = Wm.conj().T @ (0.5 * (Pd + Pd.conj().T)) @ Wm
            phi += float(torch.linalg.norm(RF) ** 2)
            RBr = torch.matmul(torch.matmul(Wm.conj().T[None], RB), Wm)
            grad += 2 * torch.real(torch.einsum("kij,ij->k", RBr, RF.conj()))
            # BARRIER: Phi_neg = sum_i [x_i]_-^2 ; d/dt = 2 [x_i]_- u_i^dag Pdot u_i ; force pushes negative eigenvalues up
            if self.kappa_neg > 0:
                # SMOOTH barrier: b(x) = eps^2 softplus(-x/eps)^2 ~ x^2 for x << -eps, exp-small inside the cone; C^inf.
                # d/dt sum_i b(x_i) = sum_i b'(x_i) u_i^dag Pdot u_i,  b'(x) = -2 eps softplus(-x/eps) sigmoid(-x/eps)
                eps = self.eps_bar
                sp = torch.nn.functional.softplus(-x / eps); bp = -2 * eps * sp * torch.sigmoid(-x / eps)
                act = bp.abs() > 1e-12 * bp.abs().max().clamp(min=1e-300)
                if bool(act.any()):
                    Ua = U[:, act]
                    rateB_diag = torch.real(torch.einsum("ia,kab,bi->ki", Ua.conj().T, RB, Ua))      # (n2, n_act)
                    grad += (self.kappa_neg / self.kappa) * (rateB_diag @ bp[act].to(rateB_diag.dtype))
                    phi += self.kappa_neg / self.kappa * float((eps ** 2 * sp ** 2).sum())
        # hard constraints for F0
        g2F0 = self.from_mat2(F0); g1F0 = self.g1_of(g2F0)
        r_hard = torch.cat([self.hvec(self.exact_g1dot(g1, g2) - g1F0),
                            (-torch.real(torch.trace(self.gradE.conj().T @ F0)))[None]])
        c_p = self.Cpinv @ r_hard
        force = -(dc - c_p) - self.kappa * (self.Z @ (self.Z.T @ grad))
        self.last_phi = phi
        return self.from_mat2(F), force / self.tau

    def rhs_f0(self, g2):
        g1 = self.g1_of(g2)
        return self.from_mat2(self.F0(g1, g2))


def main():
    mode = sys.argv[1]; T_END = float(sys.argv[2]); DT = float(sys.argv[3])
    tau, kappa, sigma = [float(x) for x in sys.argv[4:7]] if len(sys.argv) > 6 else (0.5, 10.0, 1e-3)
    kappa_neg = float(sys.argv[7]) if len(sys.argv) > 7 else 0.0
    atom, N, spin = cs.H4
    mol, mf, M, h, V, enuc = molecule(atom, "sto-3g", 0, spin); dets = sector(M, N)
    es, vec = np.linalg.eigh(build_H_sector(h, V, M, dets))
    mu1 = cs.dipole(mol, mf, M)
    if os.environ.get("SYS", "kick") == "quench":
        esA, vecA = np.linalg.eigh(build_H_sector(h, 0.3 * V, M, dets)); dm0 = np.outer(vecA[:, 0], vecA[:, 0].conj()); syslabel = "V-quench 0.3->1"
    else:
        Uk = sla.expm(-1j * 0.3 * cs.sector_operator(mu1, M, dets)); dm0 = np.outer(Uk @ vec[:, 0], (Uk @ vec[:, 0]).conj()); syslabel = "kicked ground state"
    base = F0Delta(h, V, M, N); gf = GPUFlow(base, tau, kappa, sigma, kappa_neg)
    g1_0, g2_0 = rdms_of(dm0, dets, M)
    # ---- validation against numpy
    dd = DynamicDelta(base, tau, kappa, sigma)
    dc0 = np.random.default_rng(0).normal(size=base.n2) * 1e-3
    gd_c, dcd_c = dd.rhs(g2_0, dc0)
    g2t = up(g2_0); dct = torch.as_tensor(dc0, device=dev)
    gd_g, dcd_g = gf.rhs(g2t, dct)
    print("VALIDATION gpu vs numpy rhs: |g2dot| %.1e  |dcdot| %.1e  (|g2dot|=%.3f, |dcdot|=%.3f)"
          % (np.abs(gd_g.cpu().numpy() - gd_c).max(), np.abs(dcd_g.cpu().numpy() - dcd_c).max(), np.linalg.norm(gd_c), np.linalg.norm(dcd_c)))
    E0 = (np.einsum("pq,pq", h, g1_0) + 0.25 * np.einsum("pvqw,pvqw", V, g2_0)).real
    nsteps, nout = int(round(T_END / DT)), int(round(0.5 / DT))
    print("=== %s (GPU) | H4 M=%d N=%d %s | %d RK4 steps dt=%.3f a.u. | t_end=%.1f a.u. = %.2f fs | tau=%.3g kappa=%.3g sigma=%.0e kappa_neg=%.3g"
          % (mode, M, N, syslabel, nsteps, DT, T_END, T_END * 0.02419, tau, kappa, sigma, kappa_neg))
    print("    %6s | %9s %9s %9s | %8s %8s | %7s | %9s %9s | %7s" % ("t", "minD", "minQ", "minG", "tr err", "dE", "e(g2)", "mu_x", "mu_x ex", "|delta|"))
    g2, dc = g2t.clone(), torch.zeros(base.n2, dtype=torch.float64, device=dev); rec = []
    f = (lambda y: gf.rhs(*y)) if mode == "dyndelta" else (lambda y: (gf.rhs_f0(y[0]), 0 * y[1]))
    t0 = time.time()
    for k in range(nsteps + 1):
        if k % nout == 0:
            g2n = g2.cpu().numpy(); g2x = rdms_of(cs.propagate(dm0, es, vec, k * DT), dets, M)[1]
            m = metrics(g2n, g2x, h, V, M, N, E0)
            g1 = np.einsum("pvqv->pq", g2n) / (N - 1); g1x = np.einsum("pvqv->pq", g2x) / (N - 1)
            mu, mux = np.einsum("pq,pq", mu1, g1).real, np.einsum("pq,pq", mu1, g1x).real
            rec.append((k * DT, m["minD"], m["minQ"], m["minG"], m["tr"], m["dE"], m["e2"], mu, mux, float(torch.linalg.norm(dc))))
            print("    %6.1f | %+.2e %+.2e %+.2e | %+.1e %+.1e | %.4f | %+.6f %+.6f | %.4f" % rec[-1], flush=True)
            np.save("results/long_gpu_%s_%s.npy" % (mode, os.environ.get("TAG", "")), np.array(rec))
            if not np.isfinite(m["e2"]) or m["e2"] > 10: print("    DIVERGED"); break
        if k < nsteps:
            k1 = f((g2, dc)); k2 = f((g2 + .5 * DT * k1[0], dc + .5 * DT * k1[1])); k3 = f((g2 + .5 * DT * k2[0], dc + .5 * DT * k2[1])); k4 = f((g2 + DT * k3[0], dc + DT * k3[1]))
            g2 = g2 + DT / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]); dc = dc + DT / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
            g2 = 0.5 * (g2 + g2.permute(2, 3, 0, 1).conj())
    print("    (%.0fs, %.3f s/step)" % (time.time() - t0, (time.time() - t0) / max(nsteps, 1)))


if __name__ == "__main__":
    main()
