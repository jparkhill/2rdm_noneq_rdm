#!/usr/bin/env python3
r"""Route F, episode 9: long-time propagation and the dipole spectrum.

SYSTEM: H4 asymmetric trapezoid, sto-3g, M = 8, N = 4. Kicked ground state exp(-i 0.3 mu_x)|0>.
Closures: f0 (lambda_3 = 0 channel form), muller (muller2405 + Lackner), dyndelta (tau, kappa, sigma from argv).
Integrator RK4, dt from env DT (a.u.), T_END from env; output every 0.5 a.u. Metrics: min eig D/Q/G, tr err, dE,
dipole <mu_x>(t) vs exact, e(gamma_2) vs exact. The dipole trace is saved for Fourier analysis.
"""
import sys, os, time
import numpy as np, scipy.linalg as sla
sys.path.insert(0, "."); sys.path.insert(0, "tools")
import commutator_structure as cs
from f0delta import F0Delta, DynamicDelta
from f0delta_propagate import rdms_of, metrics, rk4
from molecular_2rdm import build_H_sector, molecule, sector

DT = float(os.environ.get("DT", "0.02")); T_END = float(os.environ.get("TEND", "50")); T_OUT = 0.5


def main():
    mode = sys.argv[1]
    atom, N, spin = cs.H4
    mol, mf, M, h, V, enuc = molecule(atom, "sto-3g", 0, spin); dets = sector(M, N)
    es, vec = np.linalg.eigh(build_H_sector(h, V, M, dets))
    mu1 = cs.dipole(mol, mf, M)
    Uk = sla.expm(-1j * 0.3 * cs.sector_operator(mu1, M, dets))
    dm0 = np.outer(Uk @ vec[:, 0], (Uk @ vec[:, 0]).conj())
    base = F0Delta(h, V, M, N)
    dd = DynamicDelta(base, *[float(x) for x in sys.argv[2:5]]) if mode == "dyndelta" else None
    g1_0, g2_0 = rdms_of(dm0, dets, M)
    E0 = (np.einsum("pq,pq", h, g1_0) + 0.25 * np.einsum("pvqw,pvqw", V, g2_0)).real
    nsteps, nout = int(round(T_END / DT)), int(round(T_OUT / DT))
    print("=== %s | H4 M=%d N=%d kicked ground state | %d RK4 steps dt=%.3f a.u. | t_end=%.1f a.u. = %.2f fs" % (mode, M, N, nsteps, DT, T_END, T_END * 0.02419))
    print("    %6s | %9s %9s %9s | %8s %8s | %7s | %9s %9s" % ("t", "minD", "minQ", "minG", "tr err", "dE", "e(g2)", "mu_x", "mu_x ex"))
    g2, dc = g2_0.copy(), np.zeros(base.n2); rec = []
    t0 = time.time()
    for k in range(nsteps + 1):
        if k % nout == 0:
            g2x = rdms_of(cs.propagate(dm0, es, vec, k * DT), dets, M)[1]
            m = metrics(g2, g2x, h, V, M, N, E0)
            g1 = np.einsum("pvqv->pq", g2) / (N - 1); g1x = np.einsum("pvqv->pq", g2x) / (N - 1)
            mu, mux = np.einsum("pq,pq", mu1, g1).real, np.einsum("pq,pq", mu1, g1x).real
            rec.append((k * DT, m["minD"], m["minQ"], m["minG"], m["tr"], m["dE"], m["e2"], mu, mux))
            print("    %6.1f | %+.2e %+.2e %+.2e | %+.1e %+.1e | %.4f | %+.6f %+.6f" % rec[-1], flush=True)
            if not np.isfinite(m["e2"]) or m["e2"] > 10: print("    DIVERGED"); break
        if k < nsteps:
            if dd is None:
                g2 = rk4(lambda y: base.g2dot(y, mode), g2, DT)
            else:
                f = lambda y: dd.rhs(*y)
                k1 = f((g2, dc)); k2 = f((g2 + .5 * DT * k1[0], dc + .5 * DT * k1[1])); k3 = f((g2 + .5 * DT * k2[0], dc + .5 * DT * k2[1])); k4 = f((g2 + DT * k3[0], dc + DT * k3[1]))
                g2 = g2 + DT / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]); dc = dc + DT / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
            g2 = 0.5 * (g2 + g2.transpose(2, 3, 0, 1).conj())
    print("    (%.0fs)" % (time.time() - t0))
    np.save("results/long_%s.npy" % mode, np.array(rec))


if __name__ == "__main__":
    main()
