#!/usr/bin/env python3
r"""Episode 11: dipole spectra from the THERMAL fixed point, bath on, closure vs exact-within-A3 reference.

SYSTEM: H4 asymmetric trapezoid, sto-3g, M=8, N=4. Bath as in tools/thermalize.py (kT, Ohmic eta=0.05 w_c=1, four
one-body couplings, Davies-secular in eig(f), f = h + W[gamma_1]).

PHASES
  1. Fixed points. Closure: loaded from results/thermal_fixedpoint_kT.npz (tools/thermalize.py, 300 a.u.).
     Reference rho: propagate the FCI Gibbs state under -i[H,.] + D_f until |rhs| < 1e-9 (same basin; the closure's
     own reference was computed from the HF determinant and lands at the same state).
  2. Kick. The reference rho is kicked at the density-operator level, rho -> U rho U^dag, U = exp(-i 0.3 mu_x) (this IS
     a wavefunction-level kick applied to each ensemble member). The closure's gamma_2 is kicked by the exact
     one-body rotation of its RDMs, gamma_2 -> (u (x) u) gamma_2 (u (x) u)^dag with u = exp(-i 0.3 mu_x) on orbitals --
     the exact image of the same unitary on any 2-RDM, so no representability is lost.
  3. Propagate both with the bath ON for T a.u.; record <mu_x>(t). Spectrum: |FT| of the mean-subtracted dipole with
     a Hann window (bath damping makes the lines Lorentzian; window length sets only the floor).
  Closures: dyndelta (F0 + barrier-delta + D2) and f0_D2 (F0 + D2, no delta). Also FCI Gibbs kicked under the
  same generator, as a second reference (how much the 8% fixed-point error moves the spectrum).
usage: thermal_spectra.py kT T_END DT
"""
import sys, os, time
import numpy as np, scipy.linalg as sla, torch
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "."); sys.path.insert(0, "tools")
import commutator_structure as cs
from commutator_k3 import pair_k2
from thermalize import Thermal, up, dev, C
from bath_rank2 import davies_onebody, sector_op_general
from bath_contract import rdm_of_rho
from molecular_2rdm import build_H_sector, molecule, sector


def main():
    kT, T_END, DT = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
    atom, N, spin = cs.H4
    mol, mf, M, h, V, enuc = molecule(atom, "sto-3g", 0, spin); dets = sector(M, N)
    Hs = build_H_sector(h, V, M, dets); es, Ue = np.linalg.eigh(Hs)
    rng = np.random.default_rng(3)
    A1 = [cs.dipole(mol, mf, M)]
    for _ in range(3):
        B = rng.normal(size=(M, M)); A1.append((0.5 * (B + B.T) / np.sqrt(M)).astype(complex))
    mu1 = A1[0]
    Epq = np.array([[sector_op_general(np.eye(M)[:, [p]] @ np.eye(M)[[q], :] + 0j, M, dets) for q in range(M)] for p in range(M)])
    Mu_sec = np.einsum("pq,pqij->ij", mu1, Epq)

    # ---- exact rho-level generator on the GPU: L_w = sum_pq l_pq E_pq is a contraction of the fixed (M^2,n,n) tensor
    Epq_t = up(Epq).reshape(M * M, len(dets), len(dets)); Hs_t = up(Hs); V_t = up(V); h_t = up(h)
    def ref_rhs(r):
        g1r = torch.einsum("kij,ji->k", Epq_t, r).reshape(M, M)
        f = (h_t + torch.einsum("prqs,rs->pq", V_t, g1r)).cpu().numpy()
        jumps, _ = davies_onebody(f, A1, kT, M)
        ls = up(np.array([l for g, l, wv in jumps])).reshape(len(jumps), M * M)          # (nj, M^2)
        gs = torch.as_tensor(np.array([g for g, l, wv in jumps]), device=dev)
        L = torch.einsum("jk,kab->jab", ls, Epq_t)                                         # (nj, n, n)
        Ld = L.conj().transpose(1, 2)
        LdL = torch.einsum("j,jab,jbc->ac", gs.to(C), Ld, L)
        jump = torch.einsum("j,jab,bc,jcd->ad", gs.to(C), L, r, Ld)
        return -1j * (Hs_t @ r - r @ Hs_t) + jump - 0.5 * (LdL @ r + r @ LdL)

    def rk4(f, y, dt):
        k1 = f(y); k2 = f(y + .5 * dt * k1); k3 = f(y + .5 * dt * k2); k4 = f(y + dt * k3); return y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

    # ---- phase 1: reference fixed point from Gibbs
    w = np.exp(-(es - es[0]) / kT); w /= w.sum(); rhoG = (Ue * w) @ Ue.conj().T
    fp_cache = "results/thermal_rho_fixedpoint_kT%.2f.npy" % kT
    if os.path.exists(fp_cache):
        rho_fp = np.load(fp_cache)
    else:
        r = up(rhoG); t0 = time.time()
        for k in range(20000):
            r = rk4(ref_rhs, r, 0.05)
            if k % 400 == 0:
                nr = float(torch.linalg.norm(ref_rhs(r))); print("  ref fixed point: t=%.0f |rhs|=%.2e" % (k * 0.05, nr), flush=True)
                if nr < 1e-9: break
        rho_fp = r.cpu().numpy(); np.save(fp_cache, rho_fp); print("  (%.0fs)" % (time.time() - t0), flush=True)
    fpz = np.load("results/thermal_fixedpoint_kT%.2f.npz" % kT)
    g2_cl, delta_cl = fpz["g2"], fpz["delta"]
    g2_fp = rdm_of_rho(rho_fp, dets, M, 2); g2G = rdm_of_rho(rhoG, dets, M, 2)
    print("fixed points: |closure - ref|/|Gibbs| = %.4f   |ref - Gibbs|/|Gibbs| = %.4f" % (np.linalg.norm(g2_cl - g2_fp) / np.linalg.norm(g2G), np.linalg.norm(g2_fp - g2G) / np.linalg.norm(g2G)))

    # ---- phase 2: kick
    Uk = sla.expm(-1j * 0.3 * Mu_sec); u = sla.expm(-1j * 0.3 * mu1)
    def kick_g2(g2):   # gamma_2[p,v,q,w] = <a+_p a+_v a_w a_q>, U a+_p U^dag = sum_a a+_a u[a,p]  => gamma_2 -> u* u* gamma_2 u u
        return np.einsum("pa,vb,abcd,qc,wd->pvqw", u.conj(), u.conj(), g2, u, u)
    # sanity: kick_g2 on the Gibbs gamma_2 equals the RDM of the kicked Gibbs rho
    chk = np.abs(kick_g2(g2G) - rdm_of_rho(Uk @ rhoG @ Uk.conj().T, dets, M, 2)).max()
    print("kick consistency (RDM of kicked rho vs rotated RDM): %.1e" % chk, flush=True)
    assert chk < 1e-10, "kick convention wrong"
    starts = {"ref: exact-within-A3 fixed point": ("rho", Uk @ rho_fp @ Uk.conj().T),
              "ref: FCI Gibbs": ("rho", Uk @ rhoG @ Uk.conj().T),
              "closure dyndelta+D2 from its fixed point": ("cl", kick_g2(g2_cl), True),
              "closure F0+D2 (no delta) from same start": ("cl", kick_g2(g2_cl), False)}
    base = type("B", (), {})(); base.h, base.V, base.M, base.N, base.P = h, V, M, N, cs.pairs(M); base.nP = len(base.P)
    base.gradE = pair_k2(h, 0 * V, M) / (2 * (N - 1)) + pair_k2(0 * h, V, M) / 2
    base.B = np.zeros((1, base.nP, base.nP)); base.QB = base.GB = np.zeros((1, 1, 1)); base.Cpinv = np.zeros((1, 1)); base.Z = np.zeros((1, 1)); base.n2 = 1
    nsteps, nout = int(round(T_END / DT)), int(round(0.5 / DT))
    print("=== kicked thermal dynamics, bath ON | H4 M=%d N=%d kT=%.2f | kick exp(-i 0.3 mu_x) | RK4 dt=%.3f, T=%.0f a.u. = %.1f fs | Hann resolution ~ %.3f Ha"
          % (M, N, kT, DT, T_END, T_END * 0.02419, 4 * np.pi / T_END), flush=True)
    tr_path = "results/thermal_spectra_kT%.2f.npy" % kT
    traces = np.load(tr_path, allow_pickle=True).item() if os.path.exists(tr_path) else {}
    for name, spec in starts.items():
        if name in traces and len(traces[name]) == int(round(T_END / 0.5)) + 1 and spec[0] == "rho":
            print("  %-45s loaded" % name, flush=True); continue
        t0 = time.time(); mu_t = []
        if spec[0] == "rho":
            r = up(spec[1]); Mu_t = up(Mu_sec)
            DTR = 0.05; nr = int(round(T_END / DTR)); nor = int(round(0.5 / DTR))
            for k in range(nr + 1):
                if k % nor == 0: mu_t.append(float(torch.real(torch.trace(Mu_t @ r))))
                if k < nr: r = rk4(ref_rhs, r, DTR)
        else:
            th = Thermal(base, A1, kT, 0.5, 1e4 if spec[2] else 0.0, 1e-3)
            g2 = up(spec[1]); delta = torch.zeros((base.nP, base.nP), dtype=C, device=dev)
            f = lambda y: th.rhs(*y); mins = []
            for k in range(nsteps + 1):
                if k % nout == 0:
                    g2n = g2.cpu().numpy(); g1n = np.einsum("pvqv->pq", g2n) / (N - 1)
                    mu_t.append(np.einsum("pq,pq", mu1, g1n).real)
                    mins.append(np.linalg.eigvalsh(cs.to_mat2(g2n, M) / 2).min())
                    if not np.isfinite(mu_t[-1]): print("    DIVERGED"); break
                if k < nsteps:
                    k1 = f((g2, delta)); k2 = f((g2 + .5 * DT * k1[0], delta + .5 * DT * k1[1]))
                    k3 = f((g2 + .5 * DT * k2[0], delta + .5 * DT * k2[1])); k4 = f((g2 + DT * k3[0], delta + DT * k3[1]))
                    g2 = g2 + DT / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]); delta = delta + DT / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
                    g2 = 0.5 * (g2 + g2.permute(2, 3, 0, 1).conj()); delta = 0.5 * (delta + delta.conj().T)
            print("    %s: min eig D over run %.2e" % (name, min(mins)), flush=True)
        traces[name] = np.array(mu_t)
        print("  %-45s done in %.0fs, %d samples, <mu> range [%.4f, %.4f]" % (name, time.time() - t0, len(mu_t), min(mu_t), max(mu_t)), flush=True)
        np.save("results/thermal_spectra_kT%.2f.npy" % kT, traces, allow_pickle=True)
    # ---- spectra
    t = np.arange(len(next(iter(traces.values())))) * 0.5
    fig, axes = plt.subplots(2, 1, figsize=(11, 9))
    exc = es - es[0]; tdm = np.abs(Ue[:, 0].conj() @ Mu_sec @ Ue)
    for name, mu_t in traces.items():
        n = min(len(mu_t), len(t)); x = mu_t[:n] - mu_t[:n].mean(); win = np.hanning(n)
        F = np.abs(np.fft.rfft(x * win, n=16 * n)) * 0.5; om = 2 * np.pi * np.fft.rfftfreq(16 * n, d=0.5)
        lw = 2.2 if name.startswith("ref: exact") else 1.2
        axes[0].plot(om, F, lw=lw, label=name); axes[1].plot(t[:n] * 0.02419, mu_t[:n], lw=lw, label=name)
    for e, f_ in zip(exc, tdm):
        if f_ > 1e-3 and 0 < e < 1.8: axes[0].axvline(e, color="k", ls=":", lw=0.7, alpha=0.5)
    axes[0].set_xlim(0, 1.8); axes[0].set_xlabel("frequency (Ha)"); axes[0].set_ylabel("|FT <mu_x>| (a.u.)")
    axes[0].set_title("Kicked THERMAL state, bath on, kT=%.2f, %.1f fs: closure vs exact-within-A3 (dotted: FCI x-allowed excitations from the ground state)" % (kT, T_END * 0.02419), fontsize=10)
    axes[0].legend(fontsize=8); axes[1].set_xlabel("t (fs)"); axes[1].set_ylabel("<mu_x>(t)"); axes[1].legend(fontsize=8)
    plt.tight_layout(); plt.savefig("results/thermal_spectra_kT%.2f.png" % kT, dpi=130)
    print("saved results/thermal_spectra_kT%.2f.png" % kT)


if __name__ == "__main__":
    main()
