#!/usr/bin/env python3
r"""Route F, episode 4c: identify Phi[separable_3] EXACTLY as pair-space matrix products.

Symbols (docs/COMMUTATOR-STRUCTURE.md). Pair-space matrices over sorted pairs (p<v):
    X   = Mat_2 = 2 gamma_2,   Xd = Mat(g1 ^ g1) = 2 (g1 ^ g1),   L = Mat(lambda_2) = X - Xd
    N   = g1 (+) g1   (N[(pv),(qw)] = g1[pq] d_vw + g1[vw] d_pq - g1[pw] d_vq - g1[vq] d_pw)
    K_V[(pv),(qw)] = V[pvqw],  K_W = W (+) W (W[pq] = sum_rs V[prqs] g1[rs]),  K_h = h (+) h
    K_C = C (+) C, C[pq] = sum_rtu V[rtpu] l2[rtqu]  (non-Hermitian; kept as is)
Phi_ggg = Phi[A_3(g1 g1 g1)], Phi_gl = Phi[A_3(g1 l2)].

Each Wick term of Phi_ggg contracts V with three g1's; in pair space every such contraction is a
product of the form A . B . C with B in {K_V, K_W, K_h} (where V's indices go) and A, C in
{I, N, Xd} (how the g1's attach to the ket and bra pairs). So Phi_ggg must lie in the real span of
{P + P^dag, i(P - P^dag)} over those monomials. Same for Phi_gl with B extended by K_C and A, C by L.

The test is a least-squares fit with coefficients common across all states, and the criterion is
RESIDUAL AT ROUND-OFF, not "small". A residual of 1e-2 means the class is incomplete.
"""
import sys
import itertools
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
import commutator_structure as cs
from commutator_k3 import pair_k2
from qdiss import ansatz as az
from qdiss.eom import fullci, wick


def mats(h, V, g1, l2, M):
    I = np.eye(len(cs.pairs(M)), dtype=complex)
    N = pair_k2(g1, 0 * V, M)
    Xd = cs.to_mat2(az.wedge(g1, g1), M)
    L = cs.to_mat2(l2, M)
    W = np.einsum("prqs,rs->pq", V, g1)
    C = np.einsum("rtpu,rtqu->pq", V, l2)
    return {"I": I, "N": N, "Xd": Xd, "L": L,
            "K_V": pair_k2(0 * h, V, M), "K_W": pair_k2(W, 0 * V, M), "K_WT": pair_k2(W.T, 0 * V, M),
            "K_h": pair_k2(h, 0 * V, M),
            "K_C": pair_k2(C, 0 * V, M), "K_CT": pair_k2(C.T, 0 * V, M)}


def monomials(m, sides, cores):
    out = {}
    for a, b, c in itertools.product(sides, cores, sides):
        out["%s.%s.%s" % (a, b, c)] = m[a] @ m[b] @ m[c]
    return out


def fit(recs, target_key, sides, cores, label):
    rows, rhs = [], []
    names = None
    for r in recs:
        mons = monomials(r["m"], sides, cores)
        names = list(mons.keys())
        cols = []
        for P in mons.values():
            cols.append((P + P.conj().T).reshape(-1))
            cols.append((1j * (P - P.conj().T)).reshape(-1))
        A = np.array(cols).T
        b = r[target_key].reshape(-1)
        rows.append(np.vstack([A.real, A.imag]))
        rhs.append(np.concatenate([b.real, b.imag]))
    A = np.vstack(rows)
    b = np.concatenate(rhs)
    c, *_ = np.linalg.lstsq(A, b, rcond=1e-12)
    res = np.linalg.norm(A @ c - b) / np.linalg.norm(b)
    per = [np.linalg.norm(Ai @ c - bi) / np.linalg.norm(bi) for Ai, bi in zip(rows, rhs)]
    sv = np.linalg.svd(A, compute_uv=False)
    rank = int((sv > 1e-10 * sv[0]).sum())
    print("  %-10s over %d monomials (%d real unknowns, rank %d): common residual %.2e ; per-state %s"
          % (label, len(names), len(c), rank, res, " ".join("%.1e" % p for p in per)))
    if res < 1e-6:
        # report the sparse identity: prefer a minimum-norm representative and print nonzeros
        terms = []
        for k, n in enumerate(names):
            ch, ca = c[2 * k], c[2 * k + 1]
            if abs(ch) > 1e-8 or abs(ca) > 1e-8:
                terms.append("%+.4f (P+P^dag) %+.4f i(P-P^dag)  [P = %s]" % (ch, ca, n))
        print("      identity (min-norm representative; the monomials are not independent, rank %d of %d):" % (rank, len(c)))
        for t in terms:
            print("        " + t)
    return c, names


def main():
    recs = []
    for atom, N, spin in (cs.H3, cs.H4):
        tag = "H3" if N == 3 else "H4"
        for label, dm, dets, M, h, V in cs.states(atom, N, spin, tag, quench_times=(1.0,)):
            if "HF determinant" in label or "t=0" in label:
                continue
            H = fullci.build_H(h, V, M)
            (g1, g2, g3, d1, d2), rk = cs.ensemble_rdms(dm, dets, M, H)
            l2 = az.cumulant(2, g2, {("G", 1): g1})
            EYE = np.eye(M, dtype=complex)
            tab = {("D", 1): EYE, ("G", 1): g1, ("G", 2): g2}

            def phi_of(g3x):
                tab[("G", 3)] = np.zeros((M,) * 6, dtype=complex)
                Dm = cs.to_mat2(wick.gamma_dot(2, h, V, tab, M), M)
                tab[("G", 3)] = g3x
                return cs.to_mat2(wick.gamma_dot(2, h, V, tab, M), M) - Dm

            Phi_ggg = phi_of(az.separable(3, {("G", 1): g1, ("L", 2): 0 * l2}))
            Phi_sep = phi_of(az.separable(3, {("G", 1): g1, ("L", 2): l2}))
            recs.append(dict(label=label, Phi_ggg=Phi_ggg, Phi_gl=Phi_sep - Phi_ggg, Phi=phi_of(g3),
                             m=mats(h, V, g1, l2, M)))
            print("state: %s" % label)
    print()
    print("Phi_ggg:")
    fit(recs, "Phi_ggg", ("I", "N", "Xd"), ("K_V", "K_W"), "V, W")
    fit(recs, "Phi_ggg", ("I", "N", "Xd"), ("K_V", "K_WT"), "V, W^T")
    fit(recs, "Phi_ggg", ("I", "N", "Xd"), ("K_V", "K_W", "K_WT"), "V, W, W^T")
    print("Phi_gl:")
    fit(recs, "Phi_gl", ("I", "N", "Xd", "L"), ("K_V", "K_WT"), "V, W^T")
    fit(recs, "Phi_gl", ("I", "N", "Xd", "L"), ("K_V", "K_WT", "K_C", "K_CT"), "V, W^T, C, C^T")
    print("Phi (exact, includes lambda_3) for reference:")
    fit(recs, "Phi", ("I", "N", "Xd", "L"), ("K_V", "K_WT", "K_C", "K_CT"), "V, W^T, C, C^T")


if __name__ == "__main__":
    main()
