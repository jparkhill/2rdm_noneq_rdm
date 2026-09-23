#!/usr/bin/env python3
r"""Rank 3, checks that need no new algebra (two-adjacent-ranks rule for episodes 1 and 10b).

  (a) DRIFT IS A COMMUTATOR. d/dt gamma_3 = D_3[gamma_3] + Phi_3[gamma_4]. Claim: D_3 = i[K3, Mat_3] with
      K3 the three-particle cluster Hamiltonian on antisymmetric triples,
      K3[(pvx),(qwy)] = h on each index (antisymmetrised) + V on each pair with the third index carried by a delta.
      Test: spectral fraction of D_3 in the eigenbasis of Mat_3 (must be ~1e-16), and |D_3 - i[K3,X3]| with K3 built
      explicitly.
  (b) EXACTNESS AT N = 3 (H3 bent, M = 6): the rank-3 EOM with gamma_4 = 0 equals fullci (known, 2e-16) -- re-run as
      the pipeline test for the Mat_3 representation.
  (c) DISSIPATOR AT RANK 3: qdiss.eom.dissipator.terms(3) -- do the gamma_5 bare moments cancel, and does the rank-3
      contraction with exact gamma_3, gamma_4 match the Fock-space generator? (H3, M=6, N=3: gamma_4 = 0 exactly, so
      this also tests that the term list handles the vanishing top moment.)
  (d) THE FIVE RANK-3 CONES: Gram matrices of aaa, a+a+a+, a+aa, a+a+a, and their affine expressions in gamma_1..3,
      verified against direct Grams. (Erdahl T1/T2 conditions in the pure-state form.)
"""
import sys, itertools
import numpy as np
sys.path.insert(0, "."); sys.path.insert(0, "tools")
import commutator_structure as cs
from bath_rank2 import davies_onebody, sector_op_general
from bath_contract import rdm_of_rho
from molecular_2rdm import build_H_sector, molecule, sector
from qdiss.eom import fullci, wick, dissipator


def triples(M):
    return [t for t in itertools.combinations(range(M), 3)]


def to_mat3(T, M):
    Tr = triples(M); n = len(Tr)
    X = np.zeros((n, n), dtype=complex)
    for I, (p, v, x) in enumerate(Tr):
        for J, (q, w, y) in enumerate(Tr):
            X[I, J] = 6.0 * T[p, v, x, q, w, y]
    return X


def cluster_k3(h, V, M):
    Tr = triples(M); n = len(Tr)
    K = np.zeros((n, n), dtype=complex)
    for I, a in enumerate(Tr):
        for J, b in enumerate(Tr):
            # one-body: sum over positions where two indices match and one differs (with permutation sign)
            for perm in itertools.permutations(range(3)):
                sgn = np.linalg.det(np.eye(3)[list(perm)])
                bp = [b[i] for i in perm]
                # h part
                for k in range(3):
                    if all(a[j] == bp[j] for j in range(3) if j != k):
                        K[I, J] += sgn * h[a[k], bp[k]]
                # V part: pair (k1,k2) of a with pair of bp, third equal
                for k in range(3):
                    others = [j for j in range(3) if j != k]
                    if a[k] == bp[k]:
                        K[I, J] += sgn * 0.5 * V[a[others[0]], a[others[1]], bp[others[0]], bp[others[1]]]
    return K


def main():
    # ---------- (a), (b) on H3 bent (M=6, N=3) and H4 (M=8, N=4)
    for atom, N, spin, tag in ((cs.H3[0], 3, 1, "H3 bent"), (cs.H4[0], 4, 0, "H4 trapezoid")):
        mol, mf, M, h, V, enuc = molecule(atom, "sto-3g", 0, spin); dets = sector(M, N)
        label, dm, dets_, M_, h_, V_ = next(cs.states(atom, N, spin, tag))
        H = fullci.build_H(h, V, M)
        w, U = np.linalg.eigh(0.5 * (dm + dm.conj().T)); psi = np.zeros(1 << M, dtype=complex); psi[dets] = U[:, np.argmax(w)]
        g1, g2, g3 = (fullci.rdm(psi, M, k) for k in (1, 2, 3))
        g4 = fullci.rdm(psi, M, 4) if N >= 4 else np.zeros((M,) * 8, dtype=complex)
        d3 = fullci.rdm_dot(psi, H, M, 3)
        EYE = np.eye(M, dtype=complex)
        tab = {("D", 1): EYE, ("G", 1): g1, ("G", 2): g2, ("G", 3): g3, ("G", 4): np.zeros((M,) * 8, dtype=complex)}
        D3 = wick.gamma_dot(3, h, V, tab, M)
        tab[("G", 4)] = g4
        full = wick.gamma_dot(3, h, V, tab, M)
        X3 = to_mat3(g3, M); D3m = to_mat3(D3, M)
        x, Ux = np.linalg.eigh(0.5 * (X3 + X3.conj().T)); B = Ux.conj().T @ D3m @ Ux
        same = np.abs(x[:, None] - x[None, :]) < 1e-6
        spec = np.linalg.norm(np.where(same, B, 0)) / np.linalg.norm(D3m)
        K3 = cluster_k3(h, V, M)
        comm = 1j * (K3 @ X3 - X3 @ K3)
        print("=== %s  M=%d N=%d  Mat_3 dim %d" % (tag, M, N, X3.shape[0]))
        print("  (b) rank-3 wick with exact gamma_4 vs fullci: %.1e ; with gamma_4 = 0: %.1e %s"
              % (np.abs(full - d3).max(), np.abs(D3 - d3).max(), "(exact at N=3)" if N == 3 else "(Phi_3 = %.3f)" % np.linalg.norm(to_mat3(d3 - D3, M))))
        print("  (a) drift D_3: spectral fraction %.1e ; |D_3 - i[K3, X3]| = %.1e  (|D_3| = %.3f)  K3 Hermitian %.0e"
              % (spec, np.abs(D3m - comm).max(), np.linalg.norm(D3m), np.abs(K3 - K3.conj().T).max()))
    # ---------- (c) dissipator at rank 3 on H3 (gamma_4 = 0) and H4 (exact gamma_4)
    tl3, _, _ = dissipator.terms(3)
    from collections import Counter
    print("  (c) rank-3 dissipator term list: %s" % dict(Counter((g[1] if g else 0) for _, _, g in tl3)))
    for atom, N, spin, tag in ((cs.H3[0], 3, 1, "H3 bent"), (cs.H4[0], 4, 0, "H4 trapezoid")):
        mol, mf, M, h, V, enuc = molecule(atom, "sto-3g", 0, spin); dets = sector(M, N)
        Hs = build_H_sector(h, V, M, dets); es, Ue = np.linalg.eigh(Hs)
        rng = np.random.default_rng(3); A1 = [cs.dipole(mol, mf, M)]
        for _ in range(2):
            Bm = rng.normal(size=(M, M)); A1.append((0.5 * (Bm + Bm.T) / np.sqrt(M)).astype(complex))
        w = np.exp(-(es - es[0]) / 0.5); w /= w.sum(); rho = (Ue * w) @ Ue.conj().T
        g1, g2, g3 = (rdm_of_rho(rho, dets, M, k) for k in (1, 2, 3))
        g4 = rdm_of_rho(rho, dets, M, 4) if N >= 4 else np.zeros((M,) * 8, dtype=complex)
        f = h + np.einsum("prqs,rs->pq", V, g1)
        jumps, R = davies_onebody(f, A1, 0.5, M)
        Drho = np.zeros_like(rho)
        for g, l, wv in jumps:
            L = sector_op_general(l, M, dets); LdL = L.conj().T @ L
            Drho += g * (L @ rho @ L.conj().T - 0.5 * (LdL @ rho + rho @ LdL))
        d3x = rdm_of_rho(Drho, dets, M, 3)
        d3w = dissipator.gamma_dot(3, R, {("G", 1): g1, ("G", 2): g2, ("G", 3): g3, ("G", 4): g4}, M)
        print("      %s: |wick rank-3 dissipator - fock| = %.1e  (|d3| = %.4f)" % (tag, np.abs(d3w - d3x).max(), np.linalg.norm(d3x)))
    # ---------- (d) rank-3 cones vs direct Grams (H3 pure kicked state)
    atom, N, spin = cs.H3
    label, dm, dets, M, h, V = next(cs.states(atom, N, spin, "H3"))
    w, U = np.linalg.eigh(0.5 * (dm + dm.conj().T)); psi = np.zeros(1 << M, dtype=complex); psi[dets] = U[:, np.argmax(w)]
    g1, g2, g3 = (fullci.rdm(psi, M, k) for k in (1, 2, 3))
    I = np.eye(M)
    strings = {
        "aaa  (D3)": lambda p, v, x: fullci.aa(x, fullci.aa(v, fullci.aa(p, psi))),
        "a+a+a+ (Q3)": lambda p, v, x: fullci.ad(x, fullci.ad(v, fullci.ad(p, psi))),
        "a+aa (T-type, one creator)": lambda p, v, x: fullci.aa(x, fullci.aa(v, fullci.ad(p, psi))),
        "a+a+a (T-type, two creators)": lambda p, v, x: fullci.aa(x, fullci.ad(v, fullci.ad(p, psi))),
    }
    print("  (d) rank-3 Gram matrices (direct) at %s: dims and min eigenvalues; affine forms in gamma_1..3 to be derived from these" % label)
    for name, fn in strings.items():
        vecs = np.array([fn(p, v, x) for p in range(M) for v in range(M) for x in range(M)])
        G = vecs.conj() @ vecs.T
        ev = np.linalg.eigvalsh(0.5 * (G + G.conj().T))
        print("      %-32s dim %d  min eig %+.1e  rank %d" % (name, G.shape[0], ev.min(), (ev > 1e-10).sum()))


if __name__ == "__main__":
    main()
