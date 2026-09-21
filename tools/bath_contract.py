#!/usr/bin/env python3
r"""Episode 10, CHECK 1-2: contraction of the bath dissipator to ranks 1 and 2, and what it closes on.

Two generators on the Fock-space sector (H4, M=8, N=4), both Hermiticity- and trace-preserving:

  REDFIELD (A1 Born, A2 Markov, coupling evolved under H_S exactly, NO secular approx):
      D_R[rho] = sum_alpha ( [A~_alpha rho, A_alpha] + [A_alpha, rho A~_alpha^dag] ),
      A~ = sum_ab |a><a| A |b><b| Gamma(E_b - E_a),  Gamma(w) = gamma(w)/2 + i*Delta(w); Lamb shift Delta set to 0 here.
      A is ONE-BODY but A~ is MANY-BODY (evolved under the interacting H_S), so D_R is not quadratic in a, a^dag.
  DAVIES (A5 secular on top): the CP form of tools/bath_reference.py. Also many-body.

  A3-CLOSED REDFIELD: A~ evolved under the one-body f = h + W[gamma_1] instead of H_S. Then A~ is one-body and
      D is QUADRATIC: d/dt gamma_k|_D is linear in gamma_k and gamma_{k+1}, nothing higher (BATH-DERIVATION.md s.4).

WHAT IS MEASURED. For each generator, the exact dissipative rates d/dt gamma_1 and d/dt gamma_2 (from the Fock-space
D[rho] and the ladder RDM machinery) at a correlated state. Then a LINEAR test of rank closure: perturb the state
in a direction that changes gamma_3 but leaves gamma_1, gamma_2 fixed (a null direction of the (1,2)-contraction
maps in the space of N-electron density operators) and see whether d/dt gamma_2|_D changes. If it does, the rank-2
dissipator depends on gamma_3 (as derived for the quadratic case); if the dependence extends to gamma_4-only
directions, the closure is worse than quadratic (expected for Redfield/Davies with the exact A~).
Ground truth: rates from the Fock-space generator; nothing fitted.
"""
import sys, itertools
import numpy as np
sys.path.insert(0, "."); sys.path.insert(0, "tools")
import commutator_structure as cs
from bath_reference import DaviesBath, gamma_rate
from molecular_2rdm import build_H_sector, molecule, sector
from qdiss.eom import fullci


def rdm_of_rho(rho, dets, M, k):
    """k-RDM of a sector density matrix by eigen-ensemble."""
    w, U = np.linalg.eigh(0.5 * (rho + rho.conj().T)); out = 0
    for i in np.where(np.abs(w) > 1e-13)[0]:
        psi = np.zeros(1 << M, dtype=complex); psi[dets] = U[:, i]
        out = out + w[i] * fullci.rdm(psi, M, k)
    return out


class Redfield:
    def __init__(self, Hs, A_list, kT, one_body_f=None, M=None, dets=None):
        """one_body_f: if given (M x M), evolve A under sum f a+a (A3) instead of H_S."""
        es, U = np.linalg.eigh(Hs)
        self.es, self.U = es, U
        if one_body_f is None:
            ev, Uv = es, U
        else:
            Hf = cs.sector_operator(one_body_f, M, dets)
            ev, Uv = np.linalg.eigh(Hf)
        gaps = ev[None, :] - ev[:, None]
        G = np.vectorize(lambda w: 0.5 * gamma_rate(w, kT))(gaps)          # Gamma(w) with Lamb shift dropped
        self.pairs = []
        for A in A_list:
            Ae = Uv.conj().T @ A @ Uv
            At = Uv @ (G * Ae) @ Uv.conj().T
            self.pairs.append((A, At))

    def D(self, rho):
        out = np.zeros_like(rho)
        for A, At in self.pairs:
            out += (At @ rho @ A - A @ At @ rho) + (A @ rho @ At.conj().T - rho @ At.conj().T @ A)
        return out


def main():
    kT = 0.5
    atom, N, spin = cs.H4
    mol, mf, M, h, V, enuc = molecule(atom, "sto-3g", 0, spin); dets = sector(M, N)
    Hs = build_H_sector(h, V, M, dets)
    rng = np.random.default_rng(3)
    A_list = [cs.sector_operator(cs.dipole(mol, mf, M), M, dets)]
    for _ in range(3):
        B = rng.normal(size=(M, M)); B = 0.5 * (B + B.T) / np.sqrt(M); A_list.append(cs.sector_operator(B.astype(complex), M, dets))
    es, U = np.linalg.eigh(Hs)
    # correlated mixed state: kicked thermal
    w = np.exp(-(es - es[0]) / kT); w /= w.sum()
    import scipy.linalg as sla
    Uk = sla.expm(-1j * 0.3 * A_list[0]); rho = Uk @ (U * w) @ U.conj().T @ Uk.conj().T
    g1 = rdm_of_rho(rho, dets, M, 1)
    f = h + np.einsum("prqs,rs->pq", V, g1)
    gens = {"Redfield (exact A~)": Redfield(Hs, A_list, kT),
            "Redfield A3 (A~ under f)": Redfield(Hs, A_list, kT, one_body_f=f, M=M, dets=dets)}
    dav = DaviesBath(Hs, A_list, kT)
    gens["Davies (secular, CP)"] = type("W", (), {"D": lambda self, r, dav=dav, U=U: U @ dav.D(U.conj().T @ r @ U) @ U.conj().T})()
    # null directions: perturbations of rho that leave gamma_1 AND gamma_2 unchanged (numerically: kernel of the
    # linear map rho -> (gamma_1, gamma_2) on Hermitian traceless sector matrices)
    n = len(dets)
    basis = []
    for i in range(n):
        for j in range(i, n):
            E = np.zeros((n, n), dtype=complex)
            if i == j: E[i, i] = 1
            else: E[i, j] = E[j, i] = 1 / np.sqrt(2)
            basis.append(E)
            if i != j:
                E = np.zeros((n, n), dtype=complex); E[i, j] = 1j / np.sqrt(2); E[j, i] = -1j / np.sqrt(2); basis.append(E)
    import os
    cache = "results/bath_contract_Cmap.npy"
    if os.path.exists(cache):
        Cmap = np.load(cache); print("loaded cached contraction map %s" % str(Cmap.shape))
    else:
        print("building the (gamma_1, gamma_2) contraction map on %d Hermitian basis matrices ..." % len(basis))
        cols = []
        for E in basis:
            g1e = rdm_of_rho(E, dets, M, 1).reshape(-1); g2e = rdm_of_rho(E, dets, M, 2).reshape(-1)
            cols.append(np.concatenate([g1e.real, g1e.imag, g2e.real, g2e.imag, [np.trace(E).real]]))
        Cmap = np.array(cols).T; np.save(cache, Cmap)
    _, s, Vt = np.linalg.svd(Cmap, full_matrices=True)
    r = int((s > 1e-10 * s[0]).sum()); null = Vt[r:]
    print("rank of rho -> (tr, gamma_1, gamma_2): %d ; null directions (change gamma_3+ only): %d" % (r, null.shape[0]))
    null3 = np.zeros((0, null.shape[1]))   # gamma_4-only directions: none exist at M=8, N=4 (measured: rank 400 of 400)
    eps = 1e-4
    for name, gen in gens.items():
        D0 = gen.D(rho)
        d1 = rdm_of_rho(D0, dets, M, 1); d2 = rdm_of_rho(D0, dets, M, 2)
        print("=== %s : |d/dt gamma_1|_D = %.4f  |d/dt gamma_2|_D = %.4f  tr D = %.1e  herm %.1e"
              % (name, np.linalg.norm(d1), np.linalg.norm(d2), abs(np.trace(D0)), np.abs(D0 - D0.conj().T).max()))
        for lab, nl, kmax in (("gamma_3-only directions", null, 12), ("gamma_4-only directions", null3, 12)):
            if nl.shape[0] == 0:
                print("    (no %s available at this M, N)" % lab); continue
            resp1, resp2 = [], []
            for v in nl[rng.choice(nl.shape[0], size=min(kmax, nl.shape[0]), replace=False)]:
                dE = sum(c * E for c, E in zip(v, basis)); dE = dE / np.linalg.norm(dE)
                Dp = gen.D(rho + eps * dE)
                resp1.append(np.linalg.norm(rdm_of_rho(Dp - D0, dets, M, 1)) / eps)
                resp2.append(np.linalg.norm(rdm_of_rho(Dp - D0, dets, M, 2)) / eps)
            print("    response of d/dt gamma_1 | gamma_2 to %-24s : max %.2e | %.2e   (rate norms %.3f | %.3f)"
                  % (lab, max(resp1), max(resp2), np.linalg.norm(d1), np.linalg.norm(d2)))


if __name__ == "__main__":
    main()
