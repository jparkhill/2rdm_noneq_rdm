#!/usr/bin/env python3
r"""Episode 10b, step 3: the rank-1 and rank-2 dissipators by Wick contraction, validated against the Fock space.

Generator: Davies in the ONE-BODY eigenbasis of f = h + W[gamma_1] (CP and quadratic; docs/BATH-DERIVATION.md,
"row 4"). Jump operators L_w = sum_{eps_b - eps_a = w} A_ab a+_a a_b (in the eigenbasis of f; a+_a creates the
f-eigenorbital a), rates gamma(w) = 2 pi J(|w|)[n_B + theta(w)], Ohmic J as in tools/bath_reference.py.
For the Wick contraction all jumps are collected into R[p,q,r,s] = sum_w gamma(w) conj(l^w_pq) l^w_rs (spin-orbital basis).

CHECK A: at a correlated mixed state (kicked thermal, H4 M=8 N=4), d/dt gamma_1 and gamma_2 from the Wick term list
with EXACT gamma_2, gamma_3 must equal the RDMs of D[rho] computed on the 70-dim sector, at round-off.
CHECK B: rank-1 with gamma_2 -> gamma_1 ^ gamma_1 reproduces qdiss.rdm1's Pauli master equation (same h, no V, so f = h).
CHECK C: rank-2 with gamma_3 -> separable_3(gamma_1, lambda_2): the closed dissipator D2. Its contraction to
gamma_1 vs the exact rank-1 dissipator (should agree up to the Tr_3 lambda_3 defect), its trace, and its cone leaks.
"""
import sys
import numpy as np, scipy.linalg as sla
sys.path.insert(0, "."); sys.path.insert(0, "tools")
import commutator_structure as cs
from bath_reference import gamma_rate, TOL
from bath_contract import rdm_of_rho
from molecular_2rdm import build_H_sector, molecule, sector
from qdiss.eom import fullci, dissipator
from qdiss import ansatz as az


def davies_onebody(f, A_list, kT, M):
    """Jump operators (as M x M one-body matrices l^w in the spin-orbital basis) and rates, secular in eig(f)."""
    eps, U = np.linalg.eigh(f)
    jumps = []
    gaps = eps[None, :] - eps[:, None]                     # gaps[a,b] = eps_b - eps_a
    for A in A_list:
        Ae = U.conj().T @ A @ U
        flat = gaps.reshape(-1); order = np.argsort(flat); groups = []; cur = [order[0]]
        for i in order[1:]:
            if abs(flat[i] - flat[cur[0]]) < TOL: cur.append(i)
            else: groups.append(cur); cur = [i]
        groups.append(cur)
        for g in groups:
            w = float(np.mean(flat[g])); le = np.zeros((M, M), dtype=complex)
            for i in g:
                a, b = divmod(i, M); le[a, b] = Ae[a, b]
            if np.abs(le).max() < 1e-14: continue
            jumps.append((gamma_rate(w, kT), U @ le @ U.conj().T, w))      # back to spin-orbital basis
    R = np.zeros((M,) * 4, dtype=complex)
    for g, l, w in jumps:
        R += g * np.einsum("ab,cd->abcd", l.conj(), l)
    return jumps, R


def fock_generator(jumps, M, dets):
    """D[rho] on the N-sector for the same jumps."""
    Ls = [(g, cs.sector_operator_nonherm(l, M, dets) if hasattr(cs, "sector_operator_nonherm") else None) for g, l, w in jumps]
    return Ls


def sector_op_general(l, M, dets):
    idx = {d: i for i, d in enumerate(dets)}
    O = np.zeros((len(dets),) * 2, dtype=complex)
    for j, d in enumerate(dets):
        b = np.zeros(1 << M, dtype=complex); b[d] = 1.0
        col = np.zeros(1 << M, dtype=complex)
        for p in range(M):
            for q in range(M):
                if abs(l[p, q]) > 1e-14: col += l[p, q] * fullci.ad(p, fullci.aa(q, b))
        for d2, i in idx.items():
            if col[d2] != 0: O[i, j] = col[d2]
    return O


def main():
    kT = 0.5
    atom, N, spin = cs.H4
    mol, mf, M, h, V, enuc = molecule(atom, "sto-3g", 0, spin); dets = sector(M, N)
    Hs = build_H_sector(h, V, M, dets); es, U = np.linalg.eigh(Hs)
    rng = np.random.default_rng(3)
    A1 = [cs.dipole(mol, mf, M)]
    for _ in range(3):
        B = rng.normal(size=(M, M)); A1.append((0.5 * (B + B.T) / np.sqrt(M)).astype(complex))
    A_sec = [cs.sector_operator(A, M, dets) for A in A1]
    w = np.exp(-(es - es[0]) / kT); w /= w.sum()
    Uk = sla.expm(-1j * 0.3 * A_sec[0]); rho = Uk @ (U * w) @ U.conj().T @ Uk.conj().T
    g1, g2, g3 = (rdm_of_rho(rho, dets, M, k) for k in (1, 2, 3))
    f = h + np.einsum("prqs,rs->pq", V, g1)
    jumps, R = davies_onebody(f, A1, kT, M)
    print("Davies in eig(f): %d jump operators, all rates >= 0: %s" % (len(jumps), all(g >= 0 for g, _, _ in jumps)))
    # Fock-space D[rho]
    Lsec = [(g, sector_op_general(l, M, dets)) for g, l, wv in jumps]
    Drho = np.zeros_like(rho)
    for g, L in Lsec:
        LdL = L.conj().T @ L
        Drho += g * (L @ rho @ L.conj().T - 0.5 * (LdL @ rho + rho @ LdL))
    d1x, d2x = rdm_of_rho(Drho, dets, M, 1), rdm_of_rho(Drho, dets, M, 2)
    tab = {("G", 1): g1, ("G", 2): g2, ("G", 3): g3}
    d1w = dissipator.gamma_dot(1, R, tab, M); d2w = dissipator.gamma_dot(2, R, tab, M)
    print("CHECK A  rank 1: |wick - fock| = %.1e  (|d1| = %.4f)   rank 2: %.1e  (|d2| = %.4f)   tr D[rho] = %.1e"
          % (np.abs(d1w - d1x).max(), np.linalg.norm(d1x), np.abs(d2w - d2x).max(), np.linalg.norm(d2x), abs(np.trace(Drho))))
    # CHECK B: V = 0, gamma_2 = g1 ^ g1, compare with qdiss.rdm1 Pauli equation (blocking=True) at a diagonal-in-eig(h) gamma_1
    import qdiss.rdm1 as r1
    cls = [c for c in vars(r1).values() if isinstance(c, type) and hasattr(c, "rhs")][0]
    eps_h, Uh = np.linalg.eigh(h); nocc = 1 / (np.exp((eps_h - eps_h[N - 1] - 0.05) / 0.4) + 1)
    g1d = Uh @ np.diag(nocc) @ Uh.conj().T
    jumps0, R0 = davies_onebody(h, A1, kT, M)
    d1_w0 = dissipator.gamma_dot(1, R0, {("G", 1): g1d, ("G", 2): az.wedge(g1d, g1d)}, M)
    bath1 = cls(h, A1, 1 / kT, eta=0.05, wc=1.0)
    d1_paper = Uh @ bath1.rhs(0.0, (Uh.conj().T @ g1d @ Uh).reshape(-1)).reshape(M, M) @ Uh.conj().T
    d1_paper -= Uh @ (-1j * (eps_h[:, None] - eps_h[None, :]) * (Uh.conj().T @ g1d @ Uh)) @ Uh.conj().T   # remove coherent part
    print("CHECK B  rank-1 dissipator at lambda_2 = 0 vs qdiss.rdm1 Pauli equation (V=0, f=h): populations |diff| = %.1e, coherences |diff| = %.1e  (|d1| = %.4f)"
          % (np.abs(np.diag(Uh.conj().T @ (d1_w0 - d1_paper) @ Uh)).max(),
             np.abs((Uh.conj().T @ (d1_w0 - d1_paper) @ Uh) - np.diag(np.diag(Uh.conj().T @ (d1_w0 - d1_paper) @ Uh))).max(), np.linalg.norm(d1_paper)))
    # CHECK C: the closed rank-2 dissipator
    l2 = az.cumulant(2, g2, {("G", 1): g1})
    g3sep = az.separable(3, {("G", 1): g1, ("L", 2): l2})
    d2c = dissipator.gamma_dot(2, R, {("G", 1): g1, ("G", 2): g2, ("G", 3): g3sep}, M)
    from channel_cones import q_matrix, g_matrix
    from commutator_k3 import from_mat2
    def leaks(g1, g2, d2):
        d1c = np.einsum("pvqv->pq", d2) / (N - 1); out = []
        for Pm, Pd in ((cs.to_mat2(g2, M) / 2, cs.to_mat2(d2, M) / 2), (q_matrix(g1, g2, M), q_matrix(d1c, d2, M) - q_matrix(0 * g1, 0 * g2, M)), (g_matrix(g1, g2, M), g_matrix(d1c, d2, M))):
            x, Uu = np.linalg.eigh(0.5 * (Pm + Pm.conj().T)); out.append(x.min())
        return out
    print("CHECK C  closed D2 (gamma_3 -> separable): |D2 - exact| / |exact| = %.4f ; contraction vs exact rank-1: %.1e (rel %.3f) ; tr D2 = %.1e"
          % (np.linalg.norm(d2c - d2x) / np.linalg.norm(d2x), np.abs(np.einsum("pvqv->pq", d2c) / (N - 1) - d1x).max(),
             np.linalg.norm(np.einsum("pvqv->pq", d2c) / (N - 1) - d1x) / np.linalg.norm(d1x), abs(np.einsum("pvpv->", d2c))))
    print("         min eig of D, Q, G at this state: %s  (dissipator adds no cone kernel information here; leaks tested along propagation)" % ["%.2e" % v for v in leaks(g1, g2, d2c)])
    np.save("results/bath_rank2_R.npy", dict(R=R, f=f, A1=A1), allow_pickle=True)


if __name__ == "__main__":
    main()
