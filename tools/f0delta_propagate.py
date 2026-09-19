#!/usr/bin/env python3
r"""Route F, episode 7, step 3: the first propagation. H4 trapezoid, M = 8, N = 4.

Closures propagated on gamma_2 (gamma_1 by contraction), RK4 with fixed dt, from the same initial gamma_2 as
the exact FCI propagation (the initial state is a wavefunction/ensemble; nothing is kicked at the RDM level):
    f0delta    d/dt X = F0 + delta      (tools/f0delta.py, sigma = 1e-3, mu = 1e2, cones D, Q, G)
    f0delta_D  same with the D cone only
    f0         F0 alone (lambda_3 = 0, closed channel form)
    muller     muller2405 + Lackner reconstruction of gamma_3
Trajectories:
    A  dipole-kicked ground state (pure), kick exp(-i 0.3 mu_x)
    B  dipole-kicked thermal state kT = 0.5 (mixed)                    <- required mixed-state row
    C  interaction quench: ground state of h + 0.3 V propagated under h + V (pure)
Metrics at each output time, closure vs exact: min eigenvalue of D, Q, G (exact clearance alongside),
tr X - N(N-1), E(t) - E(0), relative errors of gamma_1, gamma_2, lambda_2 against the exact trajectory.
A dt-halving check is run on trajectory A for f0delta.
"""
import sys
import time
import numpy as np
import scipy.linalg as sla
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
import commutator_structure as cs
from commutator_k3 import from_mat2
from channel_cones import q_matrix, g_matrix
from f0delta import F0Delta
from molecular_2rdm import build_H_sector, molecule, sector
from qdiss import ansatz as az
from qdiss.eom import fullci

T_END = {"A": 2.0, "B": 2.0, "C": 2.0}
DT = 0.02
T_OUT = 0.25


def rdms_of(dm, dets, M):
    w, U = np.linalg.eigh(0.5 * (dm + dm.conj().T))
    g1 = g2 = 0
    for i in np.where(w > 1e-12)[0]:
        psi = np.zeros(1 << M, dtype=complex); psi[dets] = U[:, i]
        g1 = g1 + w[i] * fullci.rdm(psi, M, 1)
        g2 = g2 + w[i] * fullci.rdm(psi, M, 2)
    return g1, g2


def metrics(g2, g2x, h, V, M, N, E0):
    g1 = np.einsum("pvqv->pq", g2) / (N - 1)
    g1x = np.einsum("pvqv->pq", g2x) / (N - 1)
    l2, l2x = az.cumulant(2, g2, {("G", 1): g1}), az.cumulant(2, g2x, {("G", 1): g1x})
    mins = []
    for gg1, gg2 in ((g1, g2), (g1x, g2x)):
        for P in (cs.to_mat2(gg2, M) / 2, q_matrix(gg1, gg2, M), g_matrix(gg1, gg2, M)):
            mins.append(np.linalg.eigvalsh(0.5 * (P + P.conj().T)).min())
    E = (np.einsum("pq,pq", h, g1) + 0.25 * np.einsum("pvqw,pvqw", V, g2)).real
    return dict(minD=mins[0], minQ=mins[1], minG=mins[2], xminD=mins[3], xminQ=mins[4], xminG=mins[5],
                tr=np.einsum("pvpv->", g2).real - N * (N - 1), dE=E - E0,
                e1=np.linalg.norm(g1 - g1x) / np.linalg.norm(g1x), e2=np.linalg.norm(g2 - g2x) / np.linalg.norm(g2x),
                el2=np.linalg.norm(l2 - l2x) / np.linalg.norm(l2x))


def rk4(f, y, dt):
    k1 = f(y); k2 = f(y + 0.5 * dt * k1); k3 = f(y + 0.5 * dt * k2); k4 = f(y + dt * k3)
    return y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def run(label, dm0, dets, M, N, h, V, es, vec, modes, dt=DT, t_end=2.0, out=None):
    cl = F0Delta(h, V, M, N)
    g1_0, g2_0 = rdms_of(dm0, dets, M)
    E0 = (np.einsum("pq,pq", h, g1_0) + 0.25 * np.einsum("pvqw,pvqw", V, g2_0)).real
    nsteps = int(round(t_end / dt)); nout = int(round(T_OUT / dt))
    exact = {}
    for k in range(0, nsteps + 1, nout):
        exact[k] = rdms_of(cs.propagate(dm0, es, vec, k * dt), dets, M)[1]
    print("=== trajectory %s   dt=%.3f  t_end=%.1f   E0=%.6f   exact min eig at t=0: D %.1e Q %.1e G %.1e"
          % (label, dt, t_end, E0, *[metrics(g2_0, g2_0, h, V, M, N, E0)[k] for k in ("xminD", "xminQ", "xminG")]))
    header = "    %-10s %5s | %9s %9s %9s | %9s %9s %9s | %8s %8s | %7s %7s %7s"
    print(header % ("closure", "t", "minD", "minQ", "minG", "ex minD", "ex minQ", "ex minG", "tr err", "dE", "e(g1)", "e(g2)", "e(l2)"))
    results = {}
    for mode in modes:
        g2 = g2_0.copy()
        rows = []
        t0 = time.time()
        for k in range(nsteps + 1):
            if k % nout == 0:
                m = metrics(g2, exact[k], h, V, M, N, E0); m["t"] = k * dt
                rows.append(m)
                print(header % (mode, "%.2f" % (k * dt), "%+.2e" % m["minD"], "%+.2e" % m["minQ"], "%+.2e" % m["minG"],
                                "%+.2e" % m["xminD"], "%+.2e" % m["xminQ"], "%+.2e" % m["xminG"],
                                "%+.1e" % m["tr"], "%+.1e" % m["dE"], "%.4f" % m["e1"], "%.4f" % m["e2"], "%.4f" % m["el2"]))
                if not np.isfinite(m["e2"]) or m["e2"] > 10:
                    print("    %s DIVERGED" % mode); break
            if k < nsteps:
                g2 = rk4(lambda y: cl.g2dot(y, mode), g2, dt)
                g2 = 0.5 * (g2 + g2.transpose(2, 3, 0, 1).conj())        # keep Hermitian (round-off)
        print("    (%s: %.0fs)" % (mode, time.time() - t0))
        results[mode] = rows
    if out is not None:
        np.save(out, results, allow_pickle=True)
    return results


def main():
    atom, N, spin = cs.H4
    mol, mf, M, h, V, enuc = molecule(atom, "sto-3g", 0, spin)
    dets = sector(M, N)
    es, vec = np.linalg.eigh(build_H_sector(h, V, M, dets))
    Uk = sla.expm(-1j * 0.3 * cs.sector_operator(cs.dipole(mol, mf, M), M, dets))
    modes = ("f0", "muller", "f0delta_D", "f0delta")
    # A kicked pure
    dmA = np.outer(Uk @ vec[:, 0], (Uk @ vec[:, 0]).conj())
    # dt check on A, f0delta only, short
    print("### dt check (trajectory A, f0delta, t_end 0.5)")
    ra = run("A dt=0.02", dmA, dets, M, N, h, V, es, vec, ("f0delta",), dt=0.02, t_end=0.5)
    rb = run("A dt=0.01", dmA, dets, M, N, h, V, es, vec, ("f0delta",), dt=0.01, t_end=0.5)
    print("    dt-halving: e(g2) at t=0.5: %.6f (dt=0.02) vs %.6f (dt=0.01)" % (ra["f0delta"][-1]["e2"], rb["f0delta"][-1]["e2"]))
    print()
    run("A kicked ground state (pure)", dmA, dets, M, N, h, V, es, vec, modes, out="results/traj_A.npy")
    # B kicked thermal
    w = np.exp(-(es - es[0]) / 0.5); w /= w.sum()
    dmB = Uk @ (vec * w) @ vec.conj().T @ Uk.conj().T
    run("B kicked thermal kT=0.5 (mixed)", dmB, dets, M, N, h, V, es, vec, modes, out="results/traj_B.npy")
    # C quench
    esA, vecA = np.linalg.eigh(build_H_sector(h, 0.3 * V, M, dets))
    dmC = np.outer(vecA[:, 0], vecA[:, 0].conj())
    run("C V-quench 0.3 -> 1 (pure)", dmC, dets, M, N, h, V, es, vec, modes, out="results/traj_C.npy")


if __name__ == "__main__":
    main()
