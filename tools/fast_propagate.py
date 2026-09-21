#!/usr/bin/env python3
r"""Long propagation with the O(M^6) dynamical delta (tools/dyndelta_fast.py).

SYSTEM: H4 asymmetric trapezoid, sto-3g, M = 8, N = 4; initial state exp(-i 0.3 mu_x)|Psi_0> (x-dipole-kicked FCI
ground state) unless SYS=quench (ground state of h + 0.3V under h + V). RK4, dt and t_end from argv (a.u.).
Closure: d/dt X = F0 + delta, tau d/dt delta = -(delta - delta_hard) - P_perp[ kappa grad Phi + kappa_neg grad Phi_neg ],
Phi the smooth tangency potential (sigma), Phi_neg the softplus cone barrier (eps). kappa = 0 switches tangency off,
kappa_neg = 0 switches the barrier off. Exact FCI propagation alongside; output every 0.5 a.u. to results/fast_TAG.txt/.npy.

usage: python tools/fast_propagate.py TAG T_END DT tau kappa sigma kappa_neg [eps]
"""
import sys, os, time
import numpy as np, scipy.linalg as sla, torch
sys.path.insert(0, "."); sys.path.insert(0, "tools")
import commutator_structure as cs
from f0delta import F0Delta
from dyndelta_fast import FastFlow, up, dev, C
from f0delta_propagate import rdms_of, metrics
from molecular_2rdm import build_H_sector, molecule, sector


def main():
    tag, T_END, DT = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
    tau, kappa, sigma, kneg = [float(x) for x in sys.argv[4:8]]
    eps = float(sys.argv[8]) if len(sys.argv) > 8 else 1e-3
    atom, N, spin = cs.H4
    mol, mf, M, h, V, enuc = molecule(atom, "sto-3g", 0, spin); dets = sector(M, N)
    es, vec = np.linalg.eigh(build_H_sector(h, V, M, dets)); mu1 = cs.dipole(mol, mf, M)
    if os.environ.get("SYS", "kick") == "quench":
        esA, vecA = np.linalg.eigh(build_H_sector(h, 0.3 * V, M, dets)); dm0 = np.outer(vecA[:, 0], vecA[:, 0].conj()); lab = "V-quench 0.3->1"
    else:
        Uk = sla.expm(-1j * 0.3 * cs.sector_operator(mu1, M, dets)); dm0 = np.outer(Uk @ vec[:, 0], (Uk @ vec[:, 0]).conj()); lab = "x-kicked ground state"
    base = F0Delta.__new__(F0Delta)
    from commutator_k3 import pair_k2
    base.h, base.V, base.M, base.N, base.P = h, V, M, N, cs.pairs(M); base.nP = len(base.P)
    base.gradE = pair_k2(h, 0 * V, M) / (2 * (N - 1)) + pair_k2(0 * h, V, M) / 2
    base.B = np.zeros((1, base.nP, base.nP)); base.QB = base.GB = np.zeros((1, 1, 1)); base.Cpinv = np.zeros((1, 1)); base.Z = np.zeros((1, 1)); base.n2 = 1
    ff = FastFlow(base, tau, max(kappa, 1e-300), sigma, kneg); ff.eps_bar = eps
    if kappa == 0:
        ff.kappa = 0.0
    g1_0, g2_0 = rdms_of(dm0, dets, M)
    E0 = (np.einsum("pq,pq", h, g1_0) + 0.25 * np.einsum("pvqw,pvqw", V, g2_0)).real
    nsteps, nout = int(round(T_END / DT)), int(round(0.5 / DT))
    print("=== %s | H4 M=%d N=%d %s | %d RK4 steps dt=%.4f a.u. | t_end=%.1f a.u. = %.2f fs | tau=%.3g kappa=%.3g sigma=%.0e kappa_neg=%.3g eps=%.0e"
          % (tag, M, N, lab, nsteps, DT, T_END, T_END * 0.02419, tau, kappa, sigma, kneg, eps), flush=True)
    print("    %6s | %9s %9s %9s | %8s %8s | %7s | %9s %9s | %7s" % ("t", "minD", "minQ", "minG", "tr err", "dE", "e(g2)", "mu_x", "mu_x ex", "|delta|"), flush=True)
    g2 = up(g2_0); delta = torch.zeros((base.nP, base.nP), dtype=C, device=dev); rec = []
    f = lambda y: ff.rhs(*y)
    t0 = time.time()
    for k in range(nsteps + 1):
        if k % nout == 0:
            g2n = g2.cpu().numpy(); g2x = rdms_of(cs.propagate(dm0, es, vec, k * DT), dets, M)[1]
            m = metrics(g2n, g2x, h, V, M, N, E0)
            g1 = np.einsum("pvqv->pq", g2n) / (N - 1); g1x = np.einsum("pvqv->pq", g2x) / (N - 1)
            rec.append((k * DT, m["minD"], m["minQ"], m["minG"], m["tr"], m["dE"], m["e2"],
                        np.einsum("pq,pq", mu1, g1).real, np.einsum("pq,pq", mu1, g1x).real, float(torch.linalg.norm(delta))))
            print("    %6.1f | %+.2e %+.2e %+.2e | %+.1e %+.1e | %.4f | %+.6f %+.6f | %.4f" % rec[-1], flush=True)
            np.save("results/fast_%s.npy" % tag, np.array(rec))
            if not np.isfinite(m["e2"]) or m["e2"] > 10:
                print("    DIVERGED", flush=True); break
        if k < nsteps:
            k1 = f((g2, delta)); k2 = f((g2 + .5 * DT * k1[0], delta + .5 * DT * k1[1]))
            k3 = f((g2 + .5 * DT * k2[0], delta + .5 * DT * k2[1])); k4 = f((g2 + DT * k3[0], delta + DT * k3[1]))
            g2 = g2 + DT / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]); delta = delta + DT / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
            g2 = 0.5 * (g2 + g2.permute(2, 3, 0, 1).conj()); delta = 0.5 * (delta + delta.conj().T)
    print("    (%.0fs, %.3f s/step)" % (time.time() - t0, (time.time() - t0) / max(nsteps, 1)), flush=True)


if __name__ == "__main__":
    main()
