#!/usr/bin/env python3
r"""Route F, episode 6b: is the lambda_3 contribution to the flow the minimal D/Q/G tangency correction?

Symbols as in tools/channel_cones.py. At a pure state the cones D, Q, G have kernels (H3: 9, 9, 17; H4: G only, 3).
The exact flow is tangent to all three (XR theorem); the lambda_3 = 0 closed-form flow F0 leaks; the exact
remainder R3 = Phi[lambda_3] = Xdot_exact - F0 restores tangency.

QUESTION. Let delta be the MINIMUM-NORM Hermitian pair-space correction such that F0 + delta has vanishing
kernel blocks on D, Q and G simultaneously (a linear constraint: for each cone, ker^dag Mdot(F0 + delta) ker = 0).
How close is delta to the true R3?  Reported: |R3|, |delta|, |delta - R3|/|R3|, the overlap <delta, R3>/(|delta||R3|),
the leak of F0 + delta on each cone (must be ~0), and the rank-2 error of F0 + delta against the exact flow versus
F0's. Two variants: delta free, and delta constrained to have zero contraction c(delta) = 0 (R3's own contraction
is ~1e-3 of |gamma_1-dot|). Also the same with only the D-kernel constraint (the XR prescription) and only D+Q.

If |delta - R3|/|R3| is small the 3-cumulant's effect on the 2-RDM flow is fixed by N-representability; if it is
~1 the tangency correction is a small and different object and R3 carries physics beyond positivity.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
import commutator_structure as cs
from commutator_k3 import from_mat2
from closed_form_sep import closed_form
from channel_cones import q_matrix, g_matrix, KER
from qdiss import ansatz as az
from qdiss.eom import fullci, wick


def herm_basis(n):
    """Real basis of n x n Hermitian matrices as a (n*n, n*n) real->complex map."""
    B = []
    for i in range(n):
        E = np.zeros((n, n), dtype=complex); E[i, i] = 1; B.append(E)
    for i in range(n):
        for j in range(i + 1, n):
            E = np.zeros((n, n), dtype=complex); E[i, j] = E[j, i] = 1 / np.sqrt(2); B.append(E)
            E = np.zeros((n, n), dtype=complex); E[i, j] = 1j / np.sqrt(2); E[j, i] = -1j / np.sqrt(2); B.append(E)
    return B


def main():
    for atom, N, spin in (cs.H3, cs.H4):
        tag = "H3" if N == 3 else "H4"
        for label, dm, dets, M, h, V in cs.states(atom, N, spin, tag, quench_times=(1.0, 2.0)):
            if "t=0" in label or "thermal" in label:
                continue
            H = fullci.build_H(h, V, M)
            (g1, g2, g3, d1, d2), rk = cs.ensemble_rdms(dm, dets, M, H)
            Dm_, Qm, Gm = cs.to_mat2(g2, M) / 2.0, q_matrix(g1, g2, M), g_matrix(g1, g2, M)
            kers = {}
            for nm, Mm in (("D", Dm_), ("Q", Qm), ("G", Gm)):
                x, U = np.linalg.eigh(0.5 * (Mm + Mm.conj().T))
                kers[nm] = U[:, x < KER]
            drift, sandwich, qterm, Phi_ph, _ = closed_form(h, V, g1, g2, M)
            F0 = drift + sandwich + qterm + Phi_ph
            Xdot = cs.to_mat2(d2, M)
            R3 = Xdot - F0
            nP = len(cs.pairs(M))
            basis = herm_basis(nP)

            def cone_rates(F):
                g2dot = from_mat2(F, M)
                g1dot = np.einsum("pvqv->pq", g2dot) / (N - 1)
                return {"D": F / 2.0, "Q": q_matrix(g1dot, g2dot, M) - q_matrix(0 * g1, 0 * g2, M),
                        "G": g_matrix(g1dot, g2dot, M)}, g1dot

            def kernel_blocks(F, cones):
                rates, g1dot = cone_rates(F)
                out = []
                for c in cones:
                    K = kers[c]
                    if K.shape[1]:
                        B = K.conj().T @ (0.5 * (rates[c] + rates[c].conj().T)) @ K
                        out.append(B.reshape(-1))
                return np.concatenate(out) if out else np.zeros(0), g1dot

            print("=== %s   kernels D %d, Q %d, G %d   |Xdot|=%.4f  |F0|=%.4f  |R3|=%.4f  |R3|/|Xdot|=%.3f"
                  % (label, kers["D"].shape[1], kers["Q"].shape[1], kers["G"].shape[1], np.linalg.norm(Xdot), np.linalg.norm(F0), np.linalg.norm(R3), np.linalg.norm(R3) / np.linalg.norm(Xdot)))
            for cones in (("D",), ("D", "Q"), ("D", "Q", "G")):
                if sum(kers[c].shape[1] for c in cones) == 0:
                    continue
                for constrain_c in (False, True):
                    # linear map delta -> (kernel blocks, [contraction])
                    b0, g10 = kernel_blocks(F0, cones)
                    cols = []
                    for E in basis:
                        kb, g1d = kernel_blocks(E, cones)
                        col = kb if not constrain_c else np.concatenate([kb, g1d.reshape(-1)])
                        cols.append(np.concatenate([col.real, col.imag]))
                    A = np.array(cols).T
                    rhs = -b0 if not constrain_c else np.concatenate([-b0, np.zeros(M * M)])
                    rhs = np.concatenate([rhs.real, rhs.imag])
                    coef, *_ = np.linalg.lstsq(A, rhs, rcond=1e-10)
                    delta = sum(c_ * E for c_, E in zip(coef, basis))
                    F1 = F0 + delta
                    leaks = []
                    for c in ("D", "Q", "G"):
                        K = kers[c]
                        r, _ = cone_rates(F1)
                        leaks.append(np.abs(K.conj().T @ (0.5 * (r[c] + r[c].conj().T)) @ K).max() / np.linalg.norm(r[c]) if K.shape[1] else float("nan"))
                    ov = np.real(np.vdot(delta.reshape(-1), R3.reshape(-1))) / (np.linalg.norm(delta) * np.linalg.norm(R3))
                    print("  tangent to %-7s %s : |delta|=%.4f  |delta-R3|/|R3|=%.3f  overlap=%.3f   leaks D %.1e Q %.1e G %.1e   rank-2 err F0 %.4f -> F0+delta %.4f"
                          % ("+".join(cones), "c(delta)=0" if constrain_c else "free      ", np.linalg.norm(delta), np.linalg.norm(delta - R3) / np.linalg.norm(R3), ov,
                             *leaks, np.linalg.norm(F0 - Xdot) / np.linalg.norm(Xdot), np.linalg.norm(F1 - Xdot) / np.linalg.norm(Xdot)))
            print()


if __name__ == "__main__":
    main()
