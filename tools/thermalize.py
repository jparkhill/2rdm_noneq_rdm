#!/usr/bin/env python3
r"""Episode 10c: does  d/dt X = F0 + delta + D2  thermalise, and to what?

SYSTEM: H4 asymmetric trapezoid, sto-3g, M = 8, N = 4. Initial state: the HF determinant (all cumulants zero).
Bath: Ohmic J = 0.05 w exp(-w), kT = 0.5, four one-body coupling operators (x-dipole + 3 random symmetric),
Davies-secular in the eigenbasis of the INSTANTANEOUS f = h + W[gamma_1] (CP at each instant, quadratic).
Closure: F0 (lambda_3 = 0 channel form) + dynamical delta (softplus barrier only, kappa_neg, eps, tau; hard
contraction + energy constraints REPLACED by contraction only -- energy is not conserved with a bath) + D2 (rank-2
Wick contraction of the dissipator with gamma_3 -> separable).

REFERENCES (exact, on the 70-dim sector): (i) Davies-in-f generator on rho with f frozen at... no: f depends on
gamma_1(t); the reference propagates rho with the same instantaneous-f jumps (the same nonlinear generator on rho),
so it is the exact-within-A3 dynamics. (ii) The FCI Gibbs state at kT.
Reported every 2 a.u.: entropy S2 = -tr X ln X (pair space), energy, min eig D/Q/G, |gamma_2 - gamma_2^(i)|,
|gamma_2 - gamma_2^Gibbs|, and the same distances for the reference (i) itself vs Gibbs.
"""
import sys, os, time
import numpy as np, torch
sys.path.insert(0, "."); sys.path.insert(0, "tools")
import commutator_structure as cs
from commutator_k3 import pair_k2
from dyndelta_fast import FastFlow, up, dev, C
from bath_rank2 import davies_onebody, sector_op_general
from bath_contract import rdm_of_rho
from channel_cones import q_matrix, g_matrix
from molecular_2rdm import build_H_sector, molecule, sector
from qdiss.eom import dissipator
from qdiss import ansatz as az

F64 = torch.float64


class Thermal(FastFlow):
    """F0 + barrier-delta + D2; contraction constraint only (no energy constraint)."""

    def __init__(self, base, A1, kT, tau, kappa_neg, eps):
        super().__init__(base, tau, 1e-300, 1e-3, kappa_neg); self.kappa = 0.0; self.eps_bar = eps
        self.A1, self.kT = A1, kT
        M = self.M
        # contraction-only constraint machinery (drop the energy row)
        self.nc = M * M
        E = torch.eye(self.nc, dtype=F64, device=dev)
        adj = torch.stack([self._C_adj_c(E[i]) for i in range(self.nc)])
        CCd = torch.stack([self._C_c(adj[i]) for i in range(self.nc)])
        self.CCd_inv_c = torch.linalg.pinv(CCd, rtol=1e-12); self._adj_c = adj
        # dissipator term list, once
        self.tl2, self.up2, self.lo2 = dissipator.terms(2)

    def _C_c(self, F):
        return self.hvec(self.g1_of(self.from_mat2(F)))

    def _C_adj_c(self, v):
        A1 = self._hvec_adj(v)
        T = self.antisym4(torch.einsum("pq,vw->pvqw", A1, self.I2)) * 2.0
        return self.to_mat2(T) / 2.0 / (self.N - 1)

    def _torch_terms(self, rank):
        """Precompile the term list into (coef, einsum string, operand keys) with deltas substituted."""
        key = "_tt%d" % rank
        if hasattr(self, key): return getattr(self, key)
        tl, upi, loi = dissipator.terms(rank); out_sub = "".join(upi) + "".join(loi); comp = []
        for coef, deltas, g in tl:
            rep = dissipator._sub(deltas); r = lambda i: rep.get(i, i)
            subs, keys = ["".join(r(i) for i in "abcd")], ["R"]
            if g is not None:
                _, m, cre, ann = g; subs.append("".join(r(i) for i in cre) + "".join(r(i) for i in ann)); keys.append(("G", m))
            for o in out_sub:
                if r(o) != o: subs.append(o + r(o)); keys.append("I")
            comp.append((coef, ",".join(subs) + "->" + out_sub, keys))
        setattr(self, key, comp); return comp

    def _gdot(self, rank, R, tab):
        tot = 0
        for coef, expr, keys in self._torch_terms(rank):
            ops = [R if k == "R" else self.I2 if k == "I" else tab[k] for k in keys]
            tot = tot + coef * torch.einsum(expr, *ops)
        return tot

    def D2(self, g1, g2):
        """Closed rank-2 dissipator as a pair-space matrix; jumps in numpy (M x M eigenproblem), contractions on device."""
        M = self.M
        g1n = g1.cpu().numpy()
        f = self.b.h + np.einsum("prqs,rs->pq", self.b.V, g1n)
        jumps, R = davies_onebody(f, self.A1, self.kT, M)
        Rt = up(R)
        l2 = g2 - (torch.einsum("pq,vw->pvqw", g1, g1) - torch.einsum("pw,vq->pvqw", g1, g1))
        # separable gamma_3 = A_3(g1 g1 g1) + A_3(g1 l2): use the MKNO term list once, on device
        if not hasattr(self, "_sep3"):
            from qdiss.symbolic import mkno
            self._sep3 = [(sg, ",".join("".join(u) + "".join(l) for n, r_, u, l in fct) + "->pvxqwy", [(n, r_) for n, r_, u, l in fct])
                          for sg, fct in mkno.gamma_to_cumulants(("p", "v", "x"), ("q", "w", "y")) if not any(n == "L" and r_ == 3 for n, r_, u, l in fct)]
        tabc = {("G", 1): g1, ("L", 2): l2}
        g3 = 0
        for sg, expr, keys in self._sep3:
            g3 = g3 + sg * torch.einsum(expr, *[tabc[k] for k in keys])
        d2 = self._gdot(2, Rt, {("G", 1): g1, ("G", 2): g2, ("G", 3): g3})
        d1 = self._gdot(1, Rt, {("G", 1): g1, ("G", 2): g2})
        return self.to_mat2(d2), d1

    def rhs(self, g2, delta):
        g1 = self.g1_of(g2)
        F0 = self.F0(g1, g2)
        D2m, d1_diss = self.D2(g1, g2)
        Fbase = F0 + D2m
        F = Fbase + delta
        D, Q, G = self.cones(g1, g2)
        Dd, Qd, Gd, g1F = self.cone_rates(F)
        force = torch.zeros_like(F); phi = 0.0
        eps = self.eps_bar
        for Pm, Pd, adj in ((D, Dd, self.Ddot_adj), (Q, Qd, self.Qdot_adj), (G, Gd, self.Gdot_adj)):
            x, U = torch.linalg.eigh(0.5 * (Pm + Pm.conj().T))
            sp = torch.nn.functional.softplus(-x / eps); bp = -2 * eps * sp * torch.sigmoid(-x / eps)
            A_P = self.kappa_neg * (U * bp[None, :].to(C)) @ U.conj().T
            phi += float((eps ** 2 * sp ** 2).sum())
            force = force + adj(0.5 * (A_P + A_P.conj().T))
        # contraction constraint: Tr_2 (Fbase + delta)/(N-1) = exact rank-1 EOM (Hamiltonian, no gamma_3) + exact rank-1 dissipator
        target = self.exact_g1dot(g1, g2) + d1_diss
        r = self.hvec(target - self.g1_of(self.from_mat2(Fbase)))
        delta_hard = torch.einsum("i,iab->ab", (self.CCd_inv_c @ r).to(C), self._adj_c)
        y = self.CCd_inv_c @ self._C_c(force)
        force_perp = force - torch.einsum("i,iab->ab", y.to(C), self._adj_c)
        self.last_phi = phi
        return self.from_mat2(F), (-(delta - delta_hard) - force_perp) / self.tau


def main():
    kT = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
    T_END = float(sys.argv[2]) if len(sys.argv) > 2 else 300.0
    DT = float(sys.argv[3]) if len(sys.argv) > 3 else 0.01
    kneg, eps, tau = 1e4, 1e-3, 0.5
    atom, N, spin = cs.H4
    mol, mf, M, h, V, enuc = molecule(atom, "sto-3g", 0, spin); dets = sector(M, N)
    Hs = build_H_sector(h, V, M, dets); es, Ue = np.linalg.eigh(Hs)
    rng = np.random.default_rng(3)
    A1 = [cs.dipole(mol, mf, M)]
    for _ in range(3):
        B = rng.normal(size=(M, M)); A1.append((0.5 * (B + B.T) / np.sqrt(M)).astype(complex))
    A_sec = [cs.sector_operator(A, M, dets) for A in A1]
    # Gibbs
    w = np.exp(-(es - es[0]) / kT); w /= w.sum(); rhoG = (Ue * w) @ Ue.conj().T
    g1G, g2G = rdm_of_rho(rhoG, dets, M, 1), rdm_of_rho(rhoG, dets, M, 2)
    EG = (np.einsum("pq,pq", h, g1G) + 0.25 * np.einsum("pvqw,pvqw", V, g2G)).real
    # initial: HF determinant
    hf = sum(1 << p for p in range(N)); i0 = dets.index(hf)
    rho = np.zeros((len(dets),) * 2, dtype=complex); rho[i0, i0] = 1.0
    g1_0, g2_0 = rdm_of_rho(rho, dets, M, 1), rdm_of_rho(rho, dets, M, 2)
    # closure
    base = type("B", (), {})(); base.h, base.V, base.M, base.N, base.P = h, V, M, N, cs.pairs(M); base.nP = len(base.P)
    base.gradE = pair_k2(h, 0 * V, M) / (2 * (N - 1)) + pair_k2(0 * h, V, M) / 2
    base.B = np.zeros((1, base.nP, base.nP)); base.QB = base.GB = np.zeros((1, 1, 1)); base.Cpinv = np.zeros((1, 1)); base.Z = np.zeros((1, 1)); base.n2 = 1
    th = Thermal(base, A1, kT, tau, kneg, eps)
    # reference (i): rho under -i[H,rho] + D_f(gamma_1(rho))[rho], RK4 same dt
    # sector matrices of the elementary one-body operators a+_p a_q, once (M^2 of them, 70x70)
    Epq = np.array([[sector_op_general(np.eye(M)[:, [p]] @ np.eye(M)[[q], :] + 0j, M, dets) for q in range(M)] for p in range(M)])
    def ref_rhs(r):
        g1r = np.einsum("pqij,ji->pq", Epq, r)
        f = h + np.einsum("prqs,rs->pq", V, g1r)
        jumps, _ = davies_onebody(f, A1, kT, M)
        out = -1j * (Hs @ r - r @ Hs)
        for g, l, wv in jumps:
            L = np.einsum("pq,pqij->ij", l, Epq); LdL = L.conj().T @ L
            out += g * (L @ r @ L.conj().T - 0.5 * (LdL @ r + r @ LdL))
        return out
    def S2(g2n):
        x = np.linalg.eigvalsh(cs.to_mat2(g2n, M) / 2); x = x[x > 1e-14]; return -np.sum(x * np.log(x))
    def Erg(g2n):
        g1n = np.einsum("pvqv->pq", g2n) / (N - 1); return (np.einsum("pq,pq", h, g1n) + 0.25 * np.einsum("pvqw,pvqw", V, g2n)).real
    def mins(g2n):
        g1n = np.einsum("pvqv->pq", g2n) / (N - 1)
        return [np.linalg.eigvalsh(0.5 * (P + P.conj().T)).min() for P in (cs.to_mat2(g2n, M) / 2, q_matrix(g1n, g2n, M), g_matrix(g1n, g2n, M))]
    nsteps, nout = int(round(T_END / DT)), int(round(2.0 / DT))
    print("=== thermalisation from the HF determinant | H4 M=%d N=%d | kT=%.2f Ohmic eta=0.05 w_c=1, 4 couplings, Davies in eig(f) | closure F0 + barrier-delta(kneg=%.0e eps=%.0e tau=%.2g) + D2 | RK4 dt=%.3f, t_end=%.0f a.u. = %.1f fs"
          % (M, N, kT, kneg, eps, tau, DT, T_END, T_END * 0.02419), flush=True)
    print("Gibbs: E=%.5f  S2(pair)=%.4f  min eig D/Q/G %s" % (EG, S2(g2G), ["%.1e" % v for v in mins(g2G)]), flush=True)
    print("    %6s | %8s %7s | %8s %8s %8s | %9s %9s | %8s %7s | %9s | %7s" % ("t", "E", "S2", "minD", "minQ", "minG", "|g2-ref|", "|g2-Gibbs|", "E_ref", "S2_ref", "|ref-Gibbs|", "|delta|"), flush=True)
    g2 = up(g2_0); delta = torch.zeros((base.nP, base.nP), dtype=C, device=dev); r = rho.copy()
    f = lambda y: th.rhs(*y)
    def rk4r(y, dt):
        k1 = ref_rhs(y); k2 = ref_rhs(y + .5 * dt * k1); k3 = ref_rhs(y + .5 * dt * k2); k4 = ref_rhs(y + dt * k3); return y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    # reference trajectory once (cached on disk keyed by kT, T_END, DT), coarser dt is fine for the 70-dim ODE
    cache = "results/thermalize_ref_kT%.2f_T%.0f.npz" % (kT, T_END)
    if os.path.exists(cache):
        ref_g2 = np.load(cache)["ref_g2"]; print("reference loaded from %s" % cache, flush=True)
    else:
        t0 = time.time(); DTR = 0.05; nr = int(round(T_END / DTR)); nor = int(round(2.0 / DTR)); ref_g2 = []
        for k in range(nr + 1):
            if k % nor == 0: ref_g2.append(rdm_of_rho(r, dets, M, 2))
            if k < nr: r = rk4r(r, DTR)
        ref_g2 = np.array(ref_g2); np.savez(cache, ref_g2=ref_g2); print("reference computed in %.0fs (dt=%.2f)" % (time.time() - t0, DTR), flush=True)
    t0 = time.time(); rec = []
    for k in range(nsteps + 1):
        if k % nout == 0:
            g2n = g2.cpu().numpy(); g2r = ref_g2[k // nout]
            row = (k * DT, Erg(g2n), S2(g2n), *mins(g2n), np.linalg.norm(g2n - g2r) / np.linalg.norm(g2G), np.linalg.norm(g2n - g2G) / np.linalg.norm(g2G),
                   Erg(g2r), S2(g2r), np.linalg.norm(g2r - g2G) / np.linalg.norm(g2G), float(torch.linalg.norm(delta)))
            rec.append(row)
            print("    %6.1f | %8.5f %7.4f | %+.1e %+.1e %+.1e | %9.2e %9.2e | %8.5f %7.4f | %9.2e | %7.4f" % row, flush=True)
            np.save("results/thermalize_kT%.2f.npy" % kT, np.array(rec))
            if not np.isfinite(row[1]) or row[6] > 10: print("    DIVERGED"); break
        if k < nsteps:
            k1 = f((g2, delta)); k2 = f((g2 + .5 * DT * k1[0], delta + .5 * DT * k1[1]))
            k3 = f((g2 + .5 * DT * k2[0], delta + .5 * DT * k2[1])); k4 = f((g2 + DT * k3[0], delta + DT * k3[1]))
            g2 = g2 + DT / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]); delta = delta + DT / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
            g2 = 0.5 * (g2 + g2.permute(2, 3, 0, 1).conj()); delta = 0.5 * (delta + delta.conj().T)
    print("    (%.0fs)" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
