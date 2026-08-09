#!/usr/bin/env python3
"""Does the positivity clearance track MIXEDNESS or ENERGY? A thermal family confounds them.

In tools/mixed_benchmark.py the exact 2-RDM positivity clearance grew five orders of magnitude
from the pure FCI ground state to a beta=2 thermal state, and the closure ranking inverted. But
beta controls mixedness and mean energy simultaneously, so that experiment cannot say which is
responsible. It matters: if the effect tracks MIXEDNESS it transfers to nonequilibrium states,
which are mixed but not thermal; if it tracks ENERGY it only says hot states have small cumulants.

Four families of density operator rho = sum_k w_k |Psi_k><Psi_k| over the exact eigenstates of H
in the N-electron sector, plotted against a common mixedness coordinate, the participation ratio

    effective rank  =  1 / sum_k w_k^2

which equals 1 for a pure state and equals K for an equal mixture of K states.

    thermal      w_k proportional to exp(-beta (E_k - E_0)).  Mixedness and energy rise together.
    cold-K       equal weights on the K LOWEST eigenstates. Mixedness rises, energy stays low.
    hot-K        equal weights on K eigenstates starting mid-spectrum. Same effective rank as
                 cold-K at the same K, but much higher mean energy. THIS IS THE DECISIVE PAIR:
                 if clearance depends on mixedness alone, cold-K and hot-K agree at matched K.
    rank2-far    equal mixture of the ground state and one excited state k. Effective rank is
                 fixed at 2 while the mean energy is swept, isolating energy at constant mixedness.

Reported per state: mean energy, the EXACT clearance (smallest eigenvalue of the 2-RDM pair
matrix D[I,J] = 2 gamma_2[p,v,q,w] over unordered pairs), the size of the exact cumulant, and the
RELATIVE cumulant error |lambda_2 - lambda_2_exact| / |lambda_2_exact| of two closures whose
ranking inverted: blocked (occupation-only, trace-exact) and pt (coherent amplitude, sigma=1e-2).
"""
import sys
import numpy as np
sys.path.insert(0, "."); sys.path.insert(0, "tools")
import hylleraas as H
import level_shift as LS
from molecular_2rdm import build_H_sector, sector
from qdiss import closures
from qdiss.eom import fullci

M, N, h, V, enuc = H.M, H.N, H.h, H.V, H.enuc
SIGMA, PAIRS = 1.0e-2, H.PAIRS
dets = sector(M, N)
es, vec = np.linalg.eigh(build_H_sector(h, V, M, dets))
NS = len(dets)


def rdms(w):
    keep = np.where(w > 1e-13)[0]
    g1 = np.zeros((M, M), dtype=complex)
    g2 = np.zeros((M,) * 4, dtype=complex)
    for k in keep:
        psi = np.zeros(1 << M, dtype=complex); psi[dets] = vec[:, k]
        g1 += w[k] * fullci.rdm(psi, M, 1)
        g2 += w[k] * fullci.rdm(psi, M, 2)
    return g1, g2, float((w * es).sum()) + enuc


def clearance(g2):
    A = np.array([[2.0 * g2[p, v, q, w] for (q, w) in PAIRS] for (p, v) in PAIRS])
    return np.linalg.eigvalsh(0.5 * (A + A.conj().T)).min()


def row(tag, w):
    g1, g2ex, E = rdms(w)
    l2ex = g2ex - closures.wedge(g1, g1)
    nex = np.abs(l2ex).max()
    lb = closures.lambda2(g1, "blocked")
    lp = LS.lam2_sigma(g1, SIGMA)
    print("   %-14s %8.2f %+11.6f %+11.3e %10.5f %10.4f %10.4f"
          % (tag, 1.0 / (w ** 2).sum(), E, clearance(g2ex), nex,
             np.abs(lb - l2ex).max() / nex, np.abs(lp - l2ex).max() / nex))


print("mixedness vs energy, jittered asymmetric-trapezoid H4/sto-3g, M=%d N=%d, %d eigenstates"
      % (M, N, NS))
print("   %-14s %8s %11s %11s %10s %10s %10s"
      % ("state", "eff rank", "E_mean", "clearance", "|l2_exact|", "blocked", "pt"))
print("   %s" % ("-" * 80))
for beta in (50.0, 20.0, 10.0, 5.0, 2.0, 1.0):
    w = np.exp(-beta * (es - es[0])); w /= w.sum()
    row("thermal b=%.0f" % beta, w)
print("   %s" % ("-" * 80))
for K in (2, 4, 8, 16, 32, 70):
    w = np.zeros(NS); w[:K] = 1.0 / K
    row("cold-K K=%d" % K, w)
print("   %s" % ("-" * 80))
for K in (2, 4, 8, 16, 32):
    off = 35
    w = np.zeros(NS); w[off:off + K] = 1.0 / K
    row("hot-K K=%d" % K, w)
print("   %s" % ("-" * 80))
for k in (1, 5, 20, 40, 69):
    w = np.zeros(NS); w[0] = 0.5; w[k] = 0.5
    row("rank2 k=%d" % k, w)
