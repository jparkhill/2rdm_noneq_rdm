#!/usr/bin/env python3
r"""Route F, episode 4: anatomy of the collision part Phi, and polynomial effective generators K.

Symbols (docs/COMMUTATOR-STRUCTURE.md). X = Mat_2, X[(p,v),(q,w)] = 2 gamma_2[pvqw] over sorted pairs.
Exact rank-2 flow  Xdot = D + Phi,  D = i[K_2, X],  K_2 = K_h + K_V  (K_h = h on each pair index,
K_V[(pv),(qw)] = V[p,v,q,w]).  Phi is linear in gamma_3.  c(T) = Tr_2 T / (N-1).
gamma_3 = A_3(g1 g1 g1) + A_3(g1 l2) + lambda_3   (MKNO; the two separable groups are antisymmetrised
products, `qdiss.ansatz.separable`), so

    Phi = Phi_ggg + Phi_gl + Phi_l3 ,   Phi_ggg = Phi[A_3(g1 g1 g1)],  Phi_gl = Phi[A_3(g1 l2)],  Phi_l3 = Phi[lambda_3].

STEP 1 (anatomy). For each piece: norm, spectral fraction in the eigenbasis of X (the part NO
commutator can represent), and its contraction c(.) against the exact identity c(Phi) = (N-2) c(i[K_V,X]).

STEP 2 (polynomial K). Hermitian pair-space structures built from h, V, gamma_1, lambda_2:
    W+W      : W (+) W,  W[p,q] = sum_rs V[p,r,q,s] g1[r,s]        (mean field on each index)
    Vg       : sum_s V[p,v,q,s] g1[s,w], antisymmetrised in (q,w), Hermitised   (one leg dressed)
    Vgg      : sum_rs V[p,v,r,s] g1[r,q] g1[s,w], Hermitised                    (pp/hh ladder dressing)
    gVg      : sum_rs g1[p,r] V[r,v,q,s] g1[s,w], antisymmetrised, Hermitised   (ph dressing)
    Ch+Ch, Ca+Ca : (C + C^dag)/2 and i(C - C^dag)/2 on each index, C[p,q] = sum_rtu V[r,t,p,u] l2[r,t,q,u]
               (the rank-1 collision term is the anti-Hermitian part of C)
    Vl_pp    : sum_rs V[p,v,r,s] l2[r,s,q,w], Hermitised                        (lambda_2 in the pp channel)
    Vl_ph    : sum_rs V[p,r,q,s] l2[s,v,r,w], antisymmetrised, Hermitised        (lambda_2 in the ph channel)
    plus K_h and K_V themselves.
Each target Phi_piece is fitted as sum_a c_a i[K_a, X] with coefficients COMMON across all states
(per-state fits are reported only as the floor). Residual fraction = what no member of the span gives.

STEP 3 (conditions) for the fitted K_fit of the total Phi and of Phi_sep = Phi_ggg + Phi_gl:
    rank-2 flow error |i[K_2 + K_fit, X] - Xdot_exact| / |Xdot_exact|
    (C) contraction error |c(i[K_2 + K_fit, X]) - gamma_1-dot| / |gamma_1-dot|
    dE/dt of the closure flow (physical E), against exact 0.
Reference rows: D, lambda_3 = 0, and the best commutator flow D + Phi_off (error = spectral fraction).

Ground truth qdiss.eom.fullci; wick vs fullci checked per state.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
import commutator_structure as cs
from commutator_k3 import contract, pair_k2
from hamiltonian_pair_flow import energy
from qdiss import ansatz as az
from qdiss.eom import fullci, wick

DEG_TOL = 1e-6


def pairmat(T, M):
    P = cs.pairs(M)
    K = np.zeros((len(P), len(P)), dtype=complex)
    for I, (p, v) in enumerate(P):
        for J, (q, w) in enumerate(P):
            K[I, J] = T[p, v, q, w]
    return K


def herm(K):
    return 0.5 * (K + K.conj().T)


def asym_bra(T):
    return T - T.transpose(0, 1, 3, 2)


def asym_ket(T):
    return T - T.transpose(1, 0, 2, 3)


def k_structures(h, V, g1, l2, M):
    W = np.einsum("prqs,rs->pq", V, g1)
    C = np.einsum("rtpu,rtqu->pq", V, l2)
    I = np.eye(M, dtype=complex)
    out = {}
    out["K_h"] = pair_k2(h, 0 * V, M)
    out["K_V"] = pair_k2(0 * h, V, M)
    out["W+W"] = pair_k2(W, 0 * V, M)
    out["Ch+Ch"] = pair_k2(0.5 * (C + C.conj().T), 0 * V, M)          # Hermitian part of C
    out["Ca+Ca"] = pair_k2(0.5j * (C - C.conj().T), 0 * V, M)         # i x anti-Hermitian part (the rank-1 collision term)
    out["Vg"] = herm(pairmat(asym_bra(np.einsum("pvqs,sw->pvqw", V, g1)), M))
    out["Vgg"] = herm(pairmat(np.einsum("pvrs,rq,sw->pvqw", V, g1, g1), M))
    out["gVg"] = herm(pairmat(asym_bra(asym_ket(np.einsum("pr,rvqs,sw->pvqw", g1, V, g1))), M))
    out["Vl_pp"] = herm(pairmat(np.einsum("pvrs,rsqw->pvqw", V, l2), M))
    out["Vl_ph"] = herm(pairmat(asym_bra(asym_ket(np.einsum("prqs,svrw->pvqw", V, l2))), M))
    for k, K in out.items():
        assert np.abs(K - K.conj().T).max() < 1e-12, k
    return out


def spectral_split(X, F):
    x, U = np.linalg.eigh(herm(X))
    B = U.conj().T @ F @ U
    same = np.abs(x[:, None] - x[None, :]) < DEG_TOL
    spec = U @ np.where(same, B, 0.0) @ U.conj().T
    return spec, F - spec


def main():
    records = []
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
            Xd = cs.to_mat2(d2, M)
            tab = {("D", 1): EYE, ("G", 1): g1, ("G", 2): g2}

            def phi_of(g3x):
                tab[("G", 3)] = np.zeros((M,) * 6, dtype=complex)
                Dm = cs.to_mat2(wick.gamma_dot(2, h, V, tab, M), M)
                tab[("G", 3)] = g3x
                return cs.to_mat2(wick.gamma_dot(2, h, V, tab, M), M) - Dm, Dm

            Phi, Dm = phi_of(g3)
            assert np.abs(Dm + Phi - Xd).max() < 1e-12
            sep_ggg = az.separable(3, {("G", 1): g1, ("L", 2): 0 * l2})
            sep = az.separable(3, {("G", 1): g1, ("L", 2): l2})
            Phi_ggg, _ = phi_of(sep_ggg)
            Phi_sep, _ = phi_of(sep)
            Phi_gl = Phi_sep - Phi_ggg
            Phi_l3 = Phi - Phi_sep
            Ks = k_structures(h, V, g1, l2, M)
            KV = Ks["K_V"]
            cKV = contract(1j * (KV @ X - X @ KV), M, N)
            print("=== %s   rank(rho)~%d  |Xdot|=%.4f  |D|=%.4f" % (label, rk, np.linalg.norm(Xd), np.linalg.norm(Dm)))
            print("  STEP 1  piece        |piece|  spectral  |c(piece)|  |c(piece) - (N-2)c(i[K_V,X])|")
            for nm, F in (("Phi", Phi), ("Phi_ggg", Phi_ggg), ("Phi_gl", Phi_gl), ("Phi_sep", Phi_sep), ("Phi_l3", Phi_l3)):
                sp, off = spectral_split(X, F)
                cF = contract(F, M, N)
                print("          %-10s  %.4f   %.4f    %.4f      %.2e"
                      % (nm, np.linalg.norm(F), np.linalg.norm(sp) / max(np.linalg.norm(F), 1e-30),
                         np.linalg.norm(cF), np.linalg.norm(cF - (N - 2) * cKV)))
            records.append(dict(label=label, M=M, N=N, h=h, V=V, X=X, Xd=Xd, Dm=Dm, d1=d1,
                                Phi=Phi, Phi_ggg=Phi_ggg, Phi_gl=Phi_gl, Phi_sep=Phi_sep, Phi_l3=Phi_l3, Ks=Ks))
            print()

    # ------------------------------------------------------------------ STEP 2: common fits
    names = list(records[0]["Ks"].keys())
    print("STEP 2  common-coefficient fits of  Phi_piece ~ sum_a c_a i[K_a, X]  over %d states (%s)"
          % (len(records), ", ".join(names)))
    fits = {}
    for piece in ("Phi_ggg", "Phi_gl", "Phi_sep", "Phi", "Phi_off"):
        rows, rhs, floors = [], [], []
        per_state = []
        for r in records:
            X = r["X"]
            target = r[piece] if piece != "Phi_off" else spectral_split(X, r["Phi"])[1]
            cols = [(1j * (K @ X - X @ K)).reshape(-1) for K in r["Ks"].values()]
            A = np.array(cols).T
            b = target.reshape(-1)
            rows.append(A)
            rhs.append(b)
            c_s, *_ = np.linalg.lstsq(A, b, rcond=None)
            per_state.append(np.linalg.norm(A @ c_s - b) / np.linalg.norm(b))
        A = np.vstack(rows)
        b = np.concatenate(rhs)
        c, *_ = np.linalg.lstsq(A, b, rcond=None)
        sv = np.linalg.svd(A, compute_uv=False)
        res_total = np.linalg.norm(A @ c - b) / np.linalg.norm(b)
        per = [np.linalg.norm(Ai @ c - bi) / np.linalg.norm(bi) for Ai, bi in zip(rows, rhs)]
        fits[piece] = c
        print("  target %-8s  common residual %.4f   per-state with common c: %s   per-state floor: %s   cond %.1e"
              % (piece, res_total, " ".join("%.3f" % p for p in per), " ".join("%.3f" % p for p in per_state), sv[0] / sv[-1]))
        print("      coefficients: " + "  ".join("%s=%+.3f%+.3fi" % (n, cc.real, cc.imag) for n, cc in zip(names, c)))
    print()

    # ------------------------------------------------------------------ STEP 3: conditions
    print("STEP 3  closure flows  Xdot = i[K_2 + K_fit, X]   (K_fit from the common fit of the named target)")
    print("  %-40s %-14s rank-2 err   (C) contraction err   dE/dt" % ("state", "closure"))
    for r in records:
        X, Xd, Dm, d1, M, N, h, V = r["X"], r["Xd"], r["Dm"], r["d1"], r["M"], r["N"], r["h"], r["V"]
        K2 = r["Ks"]["K_h"] + r["Ks"]["K_V"]
        flows = {"D alone": Dm, "D + Phi_sep (lambda_3=0)": Dm + r["Phi_sep"],
                 "D + Phi_off (best commutator)": Dm + spectral_split(X, r["Phi"])[1]}
        for piece in ("Phi_sep", "Phi", "Phi_off"):
            Kfit = sum(cc * K for cc, K in zip(fits[piece], r["Ks"].values()))
            Kfit = herm(Kfit) if np.abs(fits[piece].imag).max() < 1e-8 else Kfit
            flows["K_fit[%s]" % piece] = 1j * ((K2 + Kfit) @ X - X @ (K2 + Kfit))
        for nm, F in flows.items():
            print("  %-40s %-30s %.4f        %.4f          %+.1e"
                  % (r["label"] if nm == "D alone" else "", nm, np.linalg.norm(F - Xd) / np.linalg.norm(Xd),
                     np.linalg.norm(contract(F, M, N) - d1) / np.linalg.norm(d1), energy(F, h, V, M, N)))
        print()


if __name__ == "__main__":
    main()
