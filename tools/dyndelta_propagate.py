#!/usr/bin/env python3
r"""Route F, episode 8: dynamical (overdamped, Car-Parrinello-like) delta.  See DynamicDelta in tools/f0delta.py.

SYSTEM: H4 asymmetric trapezoid (0.80/1.20 parallel sides), sto-3g, M = 8 spin-orbitals, N = 4 electrons.
TIME:   t_end = 2.0 a.u. = 48.4 as = 0.048 fs (1 a.u. = 24.19 as), RK4 dt = 0.01 a.u.; delta initialised at 0.
TRAJECTORIES: A dipole-kicked ground state (pure); C interaction quench 0.3V -> V (pure).
Reported vs exact FCI propagation: min eig of D, Q, G, tr error, energy drift, relative gamma_1/gamma_2/lambda_2 errors,
and the potential Phi.
"""
import sys, time
import numpy as np, scipy.linalg as sla
sys.path.insert(0, "."); sys.path.insert(0, "tools")
import commutator_structure as cs
from f0delta import F0Delta, DynamicDelta
from f0delta_propagate import rdms_of, metrics
from molecular_2rdm import build_H_sector, molecule, sector

DT, T_END, T_OUT = float(__import__("os").environ.get("DT","0.01")), float(__import__("os").environ.get("TEND","2.0")), 0.25


def run(label, dm0, dets, M, N, h, V, es, vec, tau, kappa, sigma):
    base = F0Delta(h, V, M, N)
    dd = DynamicDelta(base, tau=tau, kappa=kappa, sigma=sigma)
    g1_0, g2_0 = rdms_of(dm0, dets, M)
    E0 = (np.einsum("pq,pq", h, g1_0) + 0.25 * np.einsum("pvqw,pvqw", V, g2_0)).real
    nsteps, nout = int(round(T_END / DT)), int(round(T_OUT / DT))
    print("=== %s | tau=%.3g kappa=%.3g sigma=%.0e | H4 M=%d N=%d, %d steps of dt=%.3f a.u., t_end=%.1f a.u.=%.3f fs"
          % (label, tau, kappa, sigma, M, N, nsteps, DT, T_END, T_END * 0.02419))
    hdr = "    %5s | %9s %9s %9s | %9s %9s | %8s %8s | %7s %7s %7s | %8s %8s"
    print(hdr % ("t", "minD", "minQ", "minG", "ex minD", "ex minG", "tr err", "dE", "e(g1)", "e(g2)", "e(l2)", "Phi", "|delta|"))
    g2, dc = g2_0.copy(), np.zeros(base.n2)
    t0 = time.time()

    def f(y):
        g, d = y
        gd, dd_ = dd.rhs(g, d)
        return (gd, dd_)
    for k in range(nsteps + 1):
        if k % nout == 0:
            g2x = rdms_of(cs.propagate(dm0, es, vec, k * DT), dets, M)[1]
            m = metrics(g2, g2x, h, V, M, N, E0)
            dd.rhs(g2, dc)
            print(hdr % ("%.2f" % (k * DT), "%+.1e" % m["minD"], "%+.1e" % m["minQ"], "%+.1e" % m["minG"], "%+.1e" % m["xminD"], "%+.1e" % m["xminG"],
                         "%+.1e" % m["tr"], "%+.1e" % m["dE"], "%.4f" % m["e1"], "%.4f" % m["e2"], "%.4f" % m["el2"], "%.1e" % dd.last_phi, "%.3f" % np.linalg.norm(dc)))
            if not np.isfinite(m["e2"]) or m["e2"] > 10:
                print("    DIVERGED"); break
        if k < nsteps:
            k1 = f((g2, dc)); k2 = f((g2 + 0.5 * DT * k1[0], dc + 0.5 * DT * k1[1]))
            k3 = f((g2 + 0.5 * DT * k2[0], dc + 0.5 * DT * k2[1])); k4 = f((g2 + DT * k3[0], dc + DT * k3[1]))
            g2 = g2 + DT / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]); dc = dc + DT / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
            g2 = 0.5 * (g2 + g2.transpose(2, 3, 0, 1).conj())
    print("    (%.0fs)" % (time.time() - t0))


def main():
    atom, N, spin = cs.H4
    mol, mf, M, h, V, enuc = molecule(atom, "sto-3g", 0, spin); dets = sector(M, N)
    es, vec = np.linalg.eigh(build_H_sector(h, V, M, dets))
    Uk = sla.expm(-1j * 0.3 * cs.sector_operator(cs.dipole(mol, mf, M), M, dets))
    dmA = np.outer(Uk @ vec[:, 0], (Uk @ vec[:, 0]).conj())
    esA, vecA = np.linalg.eigh(build_H_sector(h, 0.3 * V, M, dets)); dmC = np.outer(vecA[:, 0], vecA[:, 0].conj())
    tau, kappa, sigma = [float(x) for x in sys.argv[1:4]] if len(sys.argv) > 3 else (0.05, 1e2, 1e-3)
    which = sys.argv[4] if len(sys.argv) > 4 else "AC"
    if "A" in which: run("A kicked ground state (pure)", dmA, dets, M, N, h, V, es, vec, tau, kappa, sigma)
    if "C" in which: run("C V-quench 0.3->1 (pure)", dmC, dets, M, N, h, V, es, vec, tau, kappa, sigma)


if __name__ == "__main__":
    main()
