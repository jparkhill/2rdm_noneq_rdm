#!/usr/bin/env python3
r"""Route F, episode 5b: a READABLE identity for Phi[A_3(g1 l2)] from the exact class of tools/phi_identify_ph.py.

The class (pp products + ph products, see phi_identify_ph.py) contains Phi_gl at 4e-14 but is highly
redundant. Here an L1-regularised (IRLS) fit with coefficients common across states selects a sparse
representative, which is then re-fitted by least squares on the selected columns and verified on a
held-out state (H3 quench t=2) that took no part in the selection.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
import commutator_structure as cs
from phi_identify_ph import build, columns
from qdiss import ansatz as az
from qdiss.eom import fullci, wick


def phi_gl_of(h, V, M, g1, g2, l2):
    EYE = np.eye(M, dtype=complex)
    tab = {("D", 1): EYE, ("G", 1): g1, ("G", 2): g2, ("G", 3): np.zeros((M,) * 6, dtype=complex)}
    Dm = cs.to_mat2(wick.gamma_dot(2, h, V, tab, M), M)
    tab[("G", 3)] = az.separable(3, {("G", 1): g1, ("L", 2): l2})
    Ps = cs.to_mat2(wick.gamma_dot(2, h, V, tab, M), M) - Dm
    tab[("G", 3)] = az.separable(3, {("G", 1): g1, ("L", 2): 0 * l2})
    Pg = cs.to_mat2(wick.gamma_dot(2, h, V, tab, M), M) - Dm
    return Ps - Pg


def system(states):
    rows, rhs, names = [], [], None
    for st in states:
        A, names = columns(st["pp"], st["dress"], st["ph"], st["M"], True, True)
        b = st["Phi_gl"].reshape(-1)
        rows.append(np.vstack([A.real, A.imag]))
        rhs.append(np.concatenate([b.real, b.imag]))
    return np.vstack(rows), np.concatenate(rhs), names


def main():
    train, test = [], []
    for atom, N, spin in (cs.H3, cs.H4):
        tag = "H3" if N == 3 else "H4"
        for label, dm, dets, M, h, V in cs.states(atom, N, spin, tag, quench_times=(1.0, 2.0)):
            if "HF determinant" in label or "t=0" in label:
                continue
            H = fullci.build_H(h, V, M)
            (g1, g2, g3, d1, d2), rk = cs.ensemble_rdms(dm, dets, M, H)
            l2 = az.cumulant(2, g2, {("G", 1): g1})
            pp, dress, ph = build(h, V, g1, l2, M)
            rec = dict(label=label, M=M, Phi_gl=phi_gl_of(h, V, M, g1, g2, l2), pp=pp, dress=dress, ph=ph)
            (test if "t=2" in label else train).append(rec)
    A, b, names = system(train)
    At, bt, _ = system(test)
    scale = np.linalg.norm(A, axis=0)
    scale[scale == 0] = 1.0
    As = A / scale
    print("train states: %s ; held out: %s" % ([r["label"] for r in train], [r["label"] for r in test]))
    # IRLS for min ||c||_1 s.t. As c ~ b  (weighted ridge iterations)
    w = np.ones(As.shape[1])
    lam = 1e-6
    for it in range(40):
        Wd = 1.0 / (w + 1e-9)
        # solve (As^T As + lam diag(1/w)) c = As^T b  via the normal equations in the reduced space
        Mx = As.T @ As + lam * np.diag(Wd)
        c = np.linalg.solve(Mx, As.T @ b)
        w = np.abs(c)
        res = np.linalg.norm(As @ c - b) / np.linalg.norm(b)
        nnz = int((np.abs(c) > 1e-6 * np.abs(c).max()).sum())
        if it % 10 == 9:
            print("  IRLS it %2d: residual %.1e, nnz(>1e-6 rel) %d" % (it + 1, res, nnz))
    # select, refit, verify
    for thr in (1e-3, 1e-4, 1e-5):
        sel = np.where(np.abs(c) > thr * np.abs(c).max())[0]
        cs_, *_ = np.linalg.lstsq(As[:, sel], b, rcond=None)
        res_tr = np.linalg.norm(As[:, sel] @ cs_ - b) / np.linalg.norm(b)
        res_te = np.linalg.norm((At / scale)[:, sel] @ cs_ - bt) / np.linalg.norm(bt)
        print("  threshold %.0e: %d columns, train residual %.1e, HELD-OUT residual %.1e" % (thr, len(sel), res_tr, res_te))
        if res_te < 1e-8:
            print("  identity (coefficients on column-normalised monomials; P+P^dag or i(P-P^dag) as marked):")
            for j, v in sorted(zip(sel, cs_ / scale[sel]), key=lambda t: -abs(t[1])):
                print("    %+.5f  %-11s P = %s" % (v, "(P+P^dag)" if j % 2 == 0 else "i(P-P^dag)", names[j // 2]))
            break


if __name__ == "__main__":
    main()
