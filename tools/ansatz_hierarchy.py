#!/usr/bin/env python3
"""The hierarchy of cumulant reconstructions, scored at ranks 2 AND 3 by identical code.

NOTHING IS PROPAGATED HERE. This is a static scoring of reconstructions at fixed states, so that
the pre-investment screens in CLAUDE.md are applied before any dynamics is run. Poles, trace
errors along a trajectory and positivity errors are the separate downstream measurement.

THE OBJECT RECONSTRUCTED. lambda_k, the k-cumulant of the MKNO generalized normal ordering, with
    gamma_k[p_1..p_k, q_1..q_k] = <a+_p1 .. a+_pk a_qk .. a_q1>
    Tr_k T[p_1..p_{k-1}, q_1..q_{k-1}] = sum_y T[p_1..p_{k-1} y, q_1..q_{k-1} y]
    eta_1 = I - gamma_1,  E = gamma_1 - I/2
    anti(E, X) = E contracted into each of X's 2(k-1) slots and summed
At rank 2 the inputs are gamma_1 alone; at rank 3 they are gamma_1 and lambda_2. Definitions of
each level are in qdiss.ansatz.

THE LEVELS, defined identically at both ranks:
    null          lambda_k = 0
    scaled        lambda_k = 0 with the separable groups rescaled to restore Tr gamma_k exactly.
                  At k = 3 the entry reported is Mueller et al. arXiv:2405.05093 Eq. 14 without
                  the Lackner correction, which tools/ansatz_place_muller.py showed IS a
                  two-coefficient member of this level (fit residual 3.5e-16, lambda_3 content 0).
    trace         the minimum-Frobenius-norm lambda_k with the exact partial trace. Equals
                  closures.Corr at k = 2 and Lackner Eq. 43 / gpu.trace_inverse at k = 3, both
                  verified to 0.000% norm gap in tools/ansatz_validate.py.
    scaled+trace  `scaled` with its trace channel replaced by `trace`. At k = 3 this is
                  Mueller Eq. 14 WITH the Lackner Eq. 44/45 correction.
    muller_rdmft  k = 2 only: the RDMFT Mueller / Buijse-Baerends exchange shape used as a
                  cumulant, c = n_p n_v - sqrt(n_p n_v). A different Mueller from the above.

COLUMNS.
    trace res     max |Tr_k lambda_k^rec - exact target|. Zero is the rigorous constraint.
    rel err       ||lambda_k^rec - lambda_k^exact||_F / ||lambda_k^exact||_F. 1.000 means no
                  better than lambda_k = 0.
    captured      ||P lambda_k^exact||^2 / ||lambda_k^exact||^2 with P the orthogonal projector
                  onto the trace row space: the fraction of the exact cumulant that the
                  single-trace identity can possibly determine. A property of the state, not of
                  the ansatz, so it is printed once per rank.
    lipschitz     max over 200 random directions of ||d lambda_k|| / ||d input||, at step 1e-4.
                  Screen 1 in CLAUDE.md rejects any closure whose value here blows up.

STATES. The pure FCI ground state, and a mixed state rho = sum_n w_n |n><n| with Boltzmann
weights at kT = 0.5 hartree over the sector eigenbasis -- CLAUDE.md requires a mixed-state row
because positivity and reconstruction quality are least representative at the pure ground state.
Ground truth is qdiss.eom.fullci in the N-electron sector.
"""
import sys
import numpy as np
import torch
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from molecular_2rdm import ORB_B, build_H_sector, molecule, sector
from qdiss import ansatz as az, closures, gpu
from qdiss.eom import fullci

H4_TRAP = "H 0.00 0.00 0.0; H 0.80 0.00 0.0; H 1.05 0.72 0.0; H -0.15 0.72 0.0"
H3_BENT = "H 0.00 0.00 0.0; H 1.00 0.00 0.0; H 0.35 0.90 0.0"
SYS = ((H3_BENT, "sto-3g", 3, 0, 1, "H3 bent"), (H4_TRAP, "sto-3g", 4, 0, 0, "H4 trapezoid"))
KT = 0.5
dev = gpu.device()


def fro(x):
    return float(np.sqrt(np.abs(np.vdot(x, x)).real))


def rdms(dm, M, N, dets):
    """gamma_1, gamma_2, gamma_3 of a sector density matrix dm (dim x dim over the determinants)."""
    w, U = np.linalg.eigh(dm)
    g = [np.zeros((M,) * (2 * k), dtype=complex) for k in (1, 2, 3)]
    for i, wi in enumerate(w):
        if abs(wi) < 1e-13:
            continue
        psi = np.zeros(1 << M, dtype=complex)
        psi[dets] = U[:, i]
        for j, k in enumerate((1, 2, 3)):
            g[j] += wi * fullci.rdm(psi, M, k)
    return g


def scaled_rank2(g1, N):
    """gamma_2 = a (gamma_1 wedge gamma_1) with a fixed by Tr gamma_2 = N(N-1); lambda_2 = the
    difference from the unscaled separable part."""
    sep = az.wedge(g1, g1)
    tr = np.einsum("pvpv->", sep).real
    a = N * (N - 1) / tr
    return (a - 1.0) * sep


def rank2_set(g1, N, M):
    S = az.trace_target(2, g1, None)
    out = {"null": np.zeros((M,) * 4, dtype=complex),
           "scaled": scaled_rank2(g1, N),
           "trace": az.min_norm_trace(S, M, 2),
           "muller_rdmft": closures.lambda2(g1, "muller_rdmft")}
    out["scaled+trace"] = az.retrace(out["scaled"], S, M, 2)
    out["muller_rdmft+trace"] = az.retrace(out["muller_rdmft"], S, M, 2)
    return out


def rank3_set(g1, l2, N, M):
    g2 = l2 + az.wedge(g1, g1)
    g1t, g2t = gpu.up(g1, dev), gpu.up(g2, dev)
    tab = {("G", 1): g1, ("L", 2): l2}
    sep = az.separable(3, tab)
    S = az.trace_target(3, g1, l2)
    raw = gpu.down(gpu.muller2405_gamma3(g1t, g2t, N, M, False, ORB_B)) - sep
    cor = gpu.down(gpu.muller2405_gamma3(g1t, g2t, N, M, True, ORB_B)) - sep
    return {"null": np.zeros((M,) * 6, dtype=complex),
            "scaled": raw,
            "trace": az.min_norm_trace(S, M, 3),
            "scaled+trace": cor}


for atom, basis, N, chg, spin, label in (SYS if __name__ == "__main__" else ()):
    mol, mf, M, h, V, enuc = molecule(atom, basis, chg, spin)
    dets = sector(M, N)
    es, vec = np.linalg.eigh(build_H_sector(h, V, M, dets))
    print("=== %s   M=%d N=%d   sector %d   E_FCI %.9f ==="
          % (label, M, N, len(dets), es[0].real + enuc))

    states = []
    dm = np.zeros((len(dets),) * 2, dtype=complex)
    dm[0, 0] = 1.0
    states.append(("pure ground", vec @ dm @ vec.conj().T))
    w = np.exp(-(es - es[0]).real / KT)
    w /= w.sum()
    states.append(("thermal kT=%.1f (S=%.2f)" % (KT, -float((w * np.log(w + 1e-300)).sum())),
                   vec @ np.diag(w.astype(complex)) @ vec.conj().T))

    for sname, dmat in states:
        g1, g2, g3 = rdms(dmat, M, N, dets)
        l2 = az.cumulant(2, g2, {("G", 1): g1})
        l3 = az.cumulant(3, g3, {("G", 1): g1, ("L", 2): l2})
        occ = np.linalg.eigvalsh(0.5 * (g1 + g1.conj().T))
        print()
        print("  --- state: %s   occupations [%.4f .. %.4f] ---" % (sname, occ.min(), occ.max()))
        for k, exact, sets in ((2, l2, rank2_set(g1, N, M)), (3, l3, rank3_set(g1, l2, N, M))):
            P = az.project_trace(exact, M, k)
            cap = (fro(P) / fro(exact)) ** 2
            print("    rank %d   |lambda_%d| = %.5f   trace-identity captures %.1f%% of its norm^2"
                  % (k, k, fro(exact), 100.0 * cap))
            # The moment error is the fair column for Mueller, whose target is gamma_k, not
            # lambda_k: gamma_k^rec = separable(lower cumulants) + lambda_k^rec, and the
            # separable part is identical for every ansatz here, so the two columns differ only
            # by the normalisation ||lambda_k|| vs ||gamma_k||.
            tab = {("G", 1): g1} if k == 2 else {("G", 1): g1, ("L", 2): l2}
            sep = az.separable(k, tab)
            gex = sep + exact
            print("    %-20s %11s %9s %9s %11s"
                  % ("ansatz", "trace res", "rel err", "gam err", "|lambda_k|"))
            for nm, L in sets.items():
                print("    %-20s %11.2e %9.4f %9.4f %11.5f"
                      % (nm, az.trace_residual(L, az.trace_target(k, g1, l2 if k == 3 else None), k),
                         fro(L - exact) / fro(exact),
                         fro(sep + L - gex) / fro(gex), fro(L)))
    print()
