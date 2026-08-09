#!/usr/bin/env python3
"""Cumulant reconstructions benchmarked on MIXED states, not only the FCI ground state.

WHY. Every benchmark in this project up to now was run at the pure FCI ground state, which is the
hardest possible case for N-representability and the least representative of the nonequilibrium
target. The positivity clearance of the exact 2-RDM grows five orders of magnitude between the
pure ground state and a beta=2 thermal state, so a positivity error that is fatal at zero
temperature can sit well inside the margin where the theory is meant to be used. This file adds
the mixed rows.

THE STATE. rho(beta) = sum_k w_k |Psi_k><Psi_k| with w_k proportional to exp(-beta (E_k - E_0))
over the exact eigenstates of H in the N-electron sector. Reduced density matrices are LINEAR in
rho, so gamma_1(beta) = sum_k w_k gamma_1^(k) and likewise gamma_2, and the exact reference
energy is tr(rho H) = sum_k w_k E_k. beta = infinity recovers the pure ground state.

THE APPROXIMATION BEING TESTED. Given the exact gamma_1 of rho(beta), each closure reconstructs
lambda_2, and the reconstructed two-particle reduced density matrix is
gamma_2 = gamma_1 ^ gamma_1 + lambda_2 with (a^b)[pv,qw] = a[p,q]b[v,w] - a[p,w]b[v,q]. No
propagation and no optimisation: this isolates reconstruction error alone. Closures:

    zero           lambda_2 := 0. The mean-field/TDHF cumulant, included as the null control.
    blocked        lambda_2 = A^f + f^A with A = gamma_1 - gamma_1^2, f fixed by demanding the
                   rank-2 partial-trace identity exactly. Retired for other reasons (pole surface,
                   non-extensivity, no V-dependence) but kept as the trace-exact reference point.
    muller         Mueller/Buijse-Baerends shape, exchange weight sqrt(n_p n_v).
    pt(sigma)      the reference-free adiabatic amplitude, lambda_2 = -(Lm^T Lm + sigma^2 I)^{-1}
                   Lm^T d, with Lm the matrix of the generalised-Fock superoperator L0 and d the
                   coefficient vector of the driving D, both from the MKNO rank-2 EOM. sigma =
                   1e-2 is a FIXED part of the closure definition, making it Lipschitz in gamma_1.
    trace_plus_pt  pt(sigma) plus the ker(L0) component that best carries the partial-trace
                   identity.

METRICS. Cumulant error against the exact lambda_2 of the same rho; the partial-trace identity
residual |sum_x lambda_2[px,qx] + (gamma_1 - gamma_1^2)|; energy against tr(rho H); trace errors
in both reduced density matrices; and positivity, reported ALONGSIDE the exact clearance at the
same beta, since that is the scale a violation must be judged against.
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
SIGMA = 1.0e-2
PAIRS = H.PAIRS
BETAS = (None, 50.0, 10.0, 5.0, 2.0)          # None = pure ground state


def pair_matrix(g2):
    A = np.array([[2.0 * g2[p, v, q, w] for (q, w) in PAIRS] for (p, v) in PAIRS])
    return 0.5 * (A + A.conj().T)


dets = sector(M, N)
Hs = build_H_sector(h, V, M, dets)
es, vec = np.linalg.eigh(Hs)


def thermal(beta):
    if beta is None:
        w = np.zeros(len(dets)); w[0] = 1.0
    else:
        w = np.exp(-beta * (es - es[0])); w /= w.sum()
    keep = np.where(w > 1e-13)[0]
    g1 = np.zeros((M, M), dtype=complex)
    g2 = np.zeros((M,) * 4, dtype=complex)
    for k in keep:
        psi = np.zeros(1 << M, dtype=complex); psi[dets] = vec[:, k]
        g1 += w[k] * fullci.rdm(psi, M, 1)
        g2 += w[k] * fullci.rdm(psi, M, 2)
    return g1, g2, float((w * es).sum()) + enuc, len(keep)


def reconstruct(name, g1):
    if name == "zero":
        return np.zeros((M,) * 4, dtype=complex)
    if name in ("blocked", "muller"):
        return closures.lambda2(g1, name)
    Lm = LS.L0_matrix(g1)
    d = H.coeff(H.driving(g1))
    l2 = LS.lam2_sigma(g1, SIGMA, Lm, d)
    if name == "pt":
        return l2
    if name == "trace_plus_pt":
        return H.trace_plus_pt(g1, l2, verbose=False)
    raise ValueError(name)


NAMES = ("zero", "blocked", "muller", "pt", "trace_plus_pt")
print("mixed-state reconstruction benchmark, jittered asymmetric-trapezoid H4/sto-3g, M=%d N=%d"
      % (M, N))
print("closures reconstruct lambda_2 from the exact gamma_1 of rho(beta); pt uses sigma = %.0e"
      % SIGMA)
for beta in BETAS:
    g1, g2ex, Eex, nk = thermal(beta)
    l2ex = g2ex - closures.wedge(g1, g1)
    A = g1 - g1 @ g1
    clear = np.linalg.eigvalsh(pair_matrix(g2ex)).min()
    occ = np.linalg.eigvalsh(0.5 * (g1 + g1.conj().T))
    print()
    print("=== beta = %-8s  rank(rho) = %-3d  E_exact = %+.9f  |lambda_2_exact| = %.5f"
          % ("inf" if beta is None else "%.0f" % beta, nk, Eex, np.abs(l2ex).max()))
    print("    exact occupations %.6f .. %.6f;  EXACT 2-RDM clearance min eig D = %+.3e"
          % (occ.min(), occ.max(), clear))
    print("    %-14s %10s %10s %12s %10s %10s %11s"
          % ("closure", "|l2-l2ex|", "trace-id", "E-E_exact", "tr g1 err", "tr g2 err",
             "min eig D"))
    for nm in NAMES:
        l2 = reconstruct(nm, g1)
        g2 = closures.wedge(g1, g1) + l2
        E = float(np.real(np.einsum("pq,pq->", h, g1)
                          + 0.25 * np.einsum("pvqw,pvqw->", V, g2))) + enuc
        print("    %-14s %10.5f %10.5f %+12.6f %10.2e %10.2e %+11.3e"
              % (nm, np.abs(l2 - l2ex).max(),
                 np.abs(np.einsum("pxqx->pq", l2) + A).max(),
                 E - Eex, abs(np.trace(g1).real - N),
                 abs(np.einsum("pvpv->", g2).real - N * (N - 1)),
                 np.linalg.eigvalsh(pair_matrix(g2)).min()))
