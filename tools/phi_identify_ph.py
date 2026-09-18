#!/usr/bin/env python3
r"""Route F, episode 5a: identify Phi[A_3(g1 l2)] EXACTLY as pp-space and ph-space matrix products.

Symbols (docs/COMMUTATOR-STRUCTURE.md). gamma_2[p,v,q,w] = <a+_p a+_v a_w a_q> (ket p,v ; bra q,w).
pp representation: X, Xd, L = Mat(lambda_2), N = g1(+)g1, K_V, K_WT as before.

ph representations of a 4-index T[p,v,q,w]: a particle-hole matrix pairs one KET index with one BRA index
in each of row and column. Two inequivalent pairings ("direct" and "exchange"):
    ph_d(T)[(p,q),(v,w)] = T[p,v,q,w]          (row: ket p with bra q ; column: ket v with bra w)
    ph_x(T)[(p,w),(v,q)] = T[p,v,q,w]          (row: ket p with bra w ; column: ket v with bra q)
Each row/column pair is ordered (ket, bra), so a one-body matrix acts on the ket slot as A (x) I and on
the bra slot as I (x) A; both g1 and g1^T are allowed (the convention gamma_1[p,q] = <a+_p a_q> makes the
transpose non-trivial for complex states). Dressings: {I, g1 (x) I, g1^T (x) I, I (x) g1, I (x) g1^T}.

The Wick patterns of Phi_gl (tools output, 72 terms) are of three kinds: V sharing two indices of the same
kind (both ket or both bra) with lambda_2 -> pp products K_V . L, L . K_V with N dressings; V closing two
indices on g1 -> W-type, K_WT . L, L . K_WT; and 32 terms where V shares one ket and one bra index with
lambda_2 -> ph products  A . ph(V) . B . ph(lambda_2) . C  and  A . ph(lambda_2) . B . ph(V) . C  with
either pairing for V and for lambda_2 and dressings A, B, C from the list above.

Fit with coefficients COMMON across states; criterion is residual at ROUND-OFF. Then a greedy sparse
re-fit extracts a readable identity.
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


def ph_d(T, M):
    return T.transpose(0, 2, 1, 3).reshape(M * M, M * M)      # [(p,q),(v,w)]


def ph_d_inv(A, M):
    return A.reshape(M, M, M, M).transpose(0, 2, 1, 3)


def ph_x(T, M):
    return T.transpose(0, 3, 1, 2).reshape(M * M, M * M)      # [(p,w),(v,q)]


def ph_x_inv(A, M):
    return A.reshape(M, M, M, M).transpose(0, 2, 3, 1)


def antisym4(T):
    return 0.25 * (T - T.transpose(1, 0, 2, 3) - T.transpose(0, 1, 3, 2) + T.transpose(1, 0, 3, 2))


def build(h, V, g1, l2, M):
    I2 = np.eye(M, dtype=complex)
    pp = {"I": np.eye(len(cs.pairs(M)), dtype=complex), "N": pair_k2(g1, 0 * V, M),
          "Xd": cs.to_mat2(az.wedge(g1, g1), M), "L": cs.to_mat2(l2, M),
          "K_V": pair_k2(0 * h, V, M), "K_WT": pair_k2(np.einsum("prqs,rs->pq", V, g1).T, 0 * V, M)}
    C = np.einsum("rtpu,rtqu->pq", V, l2)          # rank-1 collision matrix, a partial trace of V.lambda_2
    for cn, Cm in (("C", C), ("CT", C.T)):
        pp["K_" + cn] = pair_k2(Cm, 0 * V, M)
        for gn, Gm in (("g", g1), ("gT", g1.T)):
            pp["K_%s%s" % (cn, gn)] = pair_k2(Cm @ Gm, 0 * V, M)
            pp["K_%s%s" % (gn, cn)] = pair_k2(Gm @ Cm, 0 * V, M)
    base = {"I": np.eye(M * M, dtype=complex), "gI": np.kron(g1, I2), "gTI": np.kron(g1.T, I2),
            "Ig": np.kron(I2, g1), "IgT": np.kron(I2, g1.T)}
    Sw = np.eye(M * M, dtype=complex).reshape(M, M, M, M).transpose(0, 1, 3, 2).reshape(M * M, M * M)  # (a,b)->(b,a)
    dress = dict(base)
    for k, A in base.items():
        dress["Sw." + k] = Sw @ A
    dress["_mid"] = list(base.keys())
    # ph_x(T) = -ph_d(T) exactly (antisymmetry in the bra pair), so only the direct pairing is needed
    ph = {"Vd": ph_d(V, M), "Ld": ph_d(l2, M)}
    return pp, dress, ph


def columns(pp, dress, ph, M, use_pp=True, use_ph=True):
    cols, names = [], []

    def add(name, P_pair):
        cols.append((P_pair + P_pair.conj().T).reshape(-1))
        cols.append((1j * (P_pair - P_pair.conj().T)).reshape(-1))
        names.append(name)

    if use_pp:
        for a, b, c in itertools.product(("I", "N", "Xd", "L"), ("K_V", "K_WT"), ("I", "N", "Xd", "L")):
            add("pp:%s.%s.%s" % (a, b, c), pp[a] @ pp[b] @ pp[c])
        ccores = [k for k in pp if k.startswith("K_C") or k.startswith("K_g") or k.startswith("K_gT")]
        for a, b, c in itertools.product(("I", "N"), ccores, ("I", "N")):
            add("pp:%s.%s.%s" % (a, b, c), pp[a] @ pp[b] @ pp[c])
    if use_ph:
        outer = [k for k in dress.keys() if k != "_mid"]
        for a, c in itertools.product(outer, repeat=2):
            for b in dress["_mid"]:
                for order in (("Vd", "Ld"), ("Ld", "Vd")):
                    Pm = dress[a] @ ph[order[0]] @ dress[b] @ ph[order[1]] @ dress[c].conj().T
                    T = antisym4(ph_d_inv(Pm, M))
                    add("ph:%s.%s.%s.%s.%s^dag" % (a, order[0], b, order[1], c), cs.to_mat2(T, M) / 2.0)
    return np.array(cols).T, names


def fit(recs, key, use_pp, use_ph, label, sparse=False):
    rows, rhs = [], []
    for r in recs:
        A, names = columns(r["pp"], r["dress"], r["ph"], r["M"], use_pp, use_ph)
        b = r[key].reshape(-1)
        rows.append(np.vstack([A.real, A.imag]))
        rhs.append(np.concatenate([b.real, b.imag]))
    A = np.vstack(rows)
    b = np.concatenate(rhs)
    scale = np.linalg.norm(A, axis=0)
    scale[scale == 0] = 1.0
    c, *_ = np.linalg.lstsq(A / scale, b, rcond=1e-11)
    c = c / scale
    sv = np.linalg.svd(A / scale, compute_uv=False)
    res = np.linalg.norm(A @ c - b) / np.linalg.norm(b)
    per = [np.linalg.norm(Ai @ c - bi) / np.linalg.norm(bi) for Ai, bi in zip(rows, rhs)]
    print("  %-8s %-9s %4d monomials, rank %4d: common residual %.2e ; per-state %s"
          % (key, label, len(names), int((sv > 1e-10 * sv[0]).sum()), res, " ".join("%.1e" % p for p in per)))
    if sparse and res < 1e-8:
        # greedy orthogonal matching pursuit for a readable representative
        sel, resid = [], b.copy()
        for _ in range(60):
            corr = np.abs(A.T @ resid) / (np.linalg.norm(A, axis=0) + 1e-300)
            j = int(np.argmax(corr))
            if j in sel:
                break
            sel.append(j)
            cs_, *_ = np.linalg.lstsq(A[:, sel], b, rcond=None)
            resid = b - A[:, sel] @ cs_
            if np.linalg.norm(resid) / np.linalg.norm(b) < 1e-10:
                break
        print("      sparse representative with %d columns, residual %.1e:" % (len(sel), np.linalg.norm(resid) / np.linalg.norm(b)))
        for j, v in sorted(zip(sel, cs_), key=lambda t: t[0]):
            kind = "(P+P^dag)" if j % 2 == 0 else "i(P-P^dag)"
            print("        %+.4f %-10s P = %s" % (v, kind, names[j // 2]))
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
            pp, dress, ph = build(h, V, g1, l2, M)
            recs.append(dict(label=label, M=M, Phi_gl=Phi_sep - Phi_ggg, Phi_l3=phi_of(g3) - Phi_sep,
                             Phi=phi_of(g3), pp=pp, dress=dress, ph=ph))
            print("state: %s" % label)
    print()
    fit(recs, "Phi_gl", True, False, "pp only")
    fit(recs, "Phi_gl", False, True, "ph only")
    fit(recs, "Phi_gl", True, True, "pp + ph", sparse=True)
    print()
    fit(recs, "Phi_l3", True, True, "pp + ph")
    fit(recs, "Phi", True, True, "pp + ph")


if __name__ == "__main__":
    main()
