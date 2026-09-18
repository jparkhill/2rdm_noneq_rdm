#!/usr/bin/env python3
r"""Route F, episode 4b: is Phi of the SANDWICH form i(S^dag X - X S) with S a polynomial in pair-space
occupation matrices and the interaction?

Symbols (docs/COMMUTATOR-STRUCTURE.md). X = Mat_2. N = g1 (+) g1 on pair space
(N[(pv),(qw)] = g1[p,q] d_vw + g1[v,w] d_pq - g1[p,w] d_vq - g1[v,q] d_pw), K_V[(pv),(qw)] = V[pvqw],
K_W = W (+) W with W[p,q] = sum_rs V[p,r,q,s] g1[r,s], K_h = h (+) h, L = Mat(lambda_2) (same
convention as X, so X = L + Mat(g1 ^ g1)), K_C = C (+) C with C[p,q] = sum_rtu V[r,t,p,u] l2[r,t,q,u]
(non-Hermitian, kept whole).

HYPOTHESIS (from expanding Phi[A_3(g1 g1 g1)] by hand). Each Wick term of Phi contracts V with three
gamma_1 factors, producing pair-space PRODUCTS such as N K_V X, X K_V N, K_W X: phase-space-filling
(Bethe-Goldstone) factors. These are not Hermitian; the flow is of the XR form

    Phi  =  i ( S^dag X - X S ) ,      S = sum_a s_a S_a ,  s_a complex

whose Hermitian part of S is a commutator with X (what a K can carry) and whose anti-Hermitian part is
the spectral motion. We fit S over the monomial basis
    {K_V, N K_V, K_V N, N K_V N, K_W, N K_W, K_W N, K_h, N K_h, K_h N, L K_V, K_V L, K_C, N K_C, K_C N}
with coefficients COMMON across all states. Targets: Phi_ggg, Phi_gl, Phi_sep, Phi. Residual at 1e-3 or
below means the structure is derived; a large residual means the monomial set is incomplete.

Then, for the fitted S: split S = K - i Gamma/2 (K, Gamma Hermitian), and report the norm of Gamma
(the dissipator-like part) against K, and the spectral fraction of the fitted flow against the target's.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
import commutator_structure as cs
from commutator_k3 import pair_k2, contract
from phi_anatomy import pairmat, spectral_split, herm
from qdiss import ansatz as az
from qdiss.eom import fullci, wick


def s_basis(h, V, g1, l2, M):
    N = pair_k2(g1, 0 * V, M)
    KV = pair_k2(0 * h, V, M)
    Kh = pair_k2(h, 0 * V, M)
    W = np.einsum("prqs,rs->pq", V, g1)
    KW = pair_k2(W, 0 * V, M)
    C = np.einsum("rtpu,rtqu->pq", V, l2)
    KC = pair_k2(C, 0 * V, M)
    L = cs.to_mat2(l2, M)
    b = {"K_V": KV, "N K_V": N @ KV, "K_V N": KV @ N, "N K_V N": N @ KV @ N,
         "K_W": KW, "N K_W": N @ KW, "K_W N": KW @ N,
         "K_h": Kh, "N K_h": N @ Kh, "K_h N": Kh @ N,
         "L K_V": L @ KV, "K_V L": KV @ L,
         "K_C": KC, "N K_C": N @ KC, "K_C N": KC @ N}
    return b


def sandwich(S, X):
    return 1j * (S.conj().T @ X - X @ S)


def main():
    recs = []
    for atom, N, spin in (cs.H3, cs.H4):
        tag = "H3" if N == 3 else "H4"
        for label, dm, dets, M, h, V in cs.states(atom, N, spin, tag, quench_times=(1.0, 2.0)):
            if "HF determinant" in label or "t=0" in label:
                continue
            H = fullci.build_H(h, V, M)
            (g1, g2, g3, d1, d2), rk = cs.ensemble_rdms(dm, dets, M, H)
            l2 = az.cumulant(2, g2, {("G", 1): g1})
            EYE = np.eye(M, dtype=complex)
            X = cs.to_mat2(g2, M)
            tab = {("D", 1): EYE, ("G", 1): g1, ("G", 2): g2}

            def phi_of(g3x):
                tab[("G", 3)] = np.zeros((M,) * 6, dtype=complex)
                Dm = cs.to_mat2(wick.gamma_dot(2, h, V, tab, M), M)
                tab[("G", 3)] = g3x
                return cs.to_mat2(wick.gamma_dot(2, h, V, tab, M), M) - Dm

            Phi = phi_of(g3)
            Phi_ggg = phi_of(az.separable(3, {("G", 1): g1, ("L", 2): 0 * l2}))
            Phi_sep = phi_of(az.separable(3, {("G", 1): g1, ("L", 2): l2}))
            recs.append(dict(label=label, X=X, M=M, N=N, Phi=Phi, Phi_ggg=Phi_ggg, Phi_gl=Phi_sep - Phi_ggg,
                             Phi_sep=Phi_sep, B=s_basis(h, V, g1, l2, M)))
    names = list(recs[0]["B"].keys())
    print("sandwich fits  Phi_piece ~ i(S^dag X - X S),  S = sum_a s_a S_a (complex s_a common over %d states)" % len(recs))
    print("basis: " + ", ".join(names))
    for piece in ("Phi_ggg", "Phi_gl", "Phi_sep", "Phi"):
        rows, rhs = [], []
        for r in recs:
            X = r["X"]
            # real-linear in (Re s, Im s): i(S^dag X - X S) with S = (a + ib) S_a
            cols = []
            for S in r["B"].values():
                cols.append(sandwich(S, X).reshape(-1))
                cols.append(sandwich(1j * S, X).reshape(-1))
            A = np.array(cols).T
            b = r[piece].reshape(-1)
            rows.append(np.vstack([A.real, A.imag]))
            rhs.append(np.concatenate([b.real, b.imag]))
        A = np.vstack(rows)
        b = np.concatenate(rhs)
        c, *_ = np.linalg.lstsq(A, b, rcond=None)
        s = c[0::2] + 1j * c[1::2]
        sv = np.linalg.svd(A, compute_uv=False)
        per = [np.linalg.norm(Ai @ c - bi) / np.linalg.norm(bi) for Ai, bi in zip(rows, rhs)]
        print("  target %-8s common residual %.2e  per-state: %s  cond %.1e"
              % (piece, np.linalg.norm(A @ c - b) / np.linalg.norm(b), " ".join("%.1e" % p for p in per), sv[0] / sv[-1]))
        keep = [(n, v) for n, v in zip(names, s) if abs(v) > 1e-6]
        print("      S = " + "  ".join("%+.4f%+.4fi [%s]" % (v.real, v.imag, n) for n, v in keep))
        # K / Gamma split of the fitted S and its spectral content, first and last state
        for r in (recs[0], recs[-1]):
            S = sum(v * r["B"][n] for n, v in zip(names, s))
            K = herm(S)
            G = 1j * (S - S.conj().T)
            F = sandwich(S, r["X"])
            sp_t = spectral_split(r["X"], r[piece])[0]
            sp_f = spectral_split(r["X"], F)[0]
            print("      %-38s |K|_2=%.3f |Gamma|_2=%.3f   spectral frac: target %.4f fit %.4f"
                  % (r["label"], np.linalg.norm(K, 2), np.linalg.norm(G, 2),
                     np.linalg.norm(sp_t) / np.linalg.norm(r[piece]), np.linalg.norm(sp_f) / max(np.linalg.norm(F), 1e-30)))
    print()


if __name__ == "__main__":
    main()
