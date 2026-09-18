#!/usr/bin/env python3
r"""Route F, episode 6c: is the lambda_3 remainder R3 = Phi[lambda_3] "inward pointing"?

For each cone (D, Q, G; symbols in tools/channel_cones.py) take the near-boundary eigenvalues
0 < x_i < 0.05 x_max and compare the eigenvalue rates from F0 (the lambda_3 = 0 closed form) and from R3:
mean[sign(F0 rate) * R3 rate], the fraction of small eigenvalues where R3 opposes F0, the fraction it raises,
and the spectral entropy production dS/dt = -sum_i xdot_i (ln x_i + 1) of F0, R3 and the exact flow.
An inward (GKSL-like) model of R3 would raise small eigenvalues and produce entropy systematically.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
import commutator_structure as cs
from commutator_k3 import from_mat2
from closed_form_sep import closed_form
from channel_cones import q_matrix, g_matrix
from qdiss.eom import fullci


def main():
    print("%-38s cone n_small  mean[sign(F0)*R3]  frac(R3 opposes F0)  frac(R3 raises)   dS/dt: F0 / R3 / exact" % "state")
    for atom, N, spin in (cs.H3, cs.H4):
        for label, dm, dets, M, h, V in cs.states(atom, N, spin, "H3" if N == 3 else "H4", quench_times=(1.0, 2.0)):
            if "t=0" in label:
                continue
            H = fullci.build_H(h, V, M)
            (g1, g2, g3, d1, d2), rk = cs.ensemble_rdms(dm, dets, M, H)
            drift, sandwich, qterm, Phi_ph, _ = closed_form(h, V, g1, g2, M)
            F0 = drift + sandwich + qterm + Phi_ph
            Xd = cs.to_mat2(d2, M)
            R3 = Xd - F0
            for cone in ("D", "Q", "G"):
                def mat(gg1, gg2, affine=True):
                    if cone == "D":
                        return cs.to_mat2(gg2, M) / 2
                    if cone == "Q":
                        return q_matrix(gg1, gg2, M) - (0 if affine else q_matrix(0 * g1, 0 * g2, M))
                    return g_matrix(gg1, gg2, M)
                Mm = mat(g1, g2)
                x, U = np.linalg.eigh(0.5 * (Mm + Mm.conj().T))
                rates = {}
                for nm, F in (("F0", F0), ("R3", R3), ("ex", Xd)):
                    gd = from_mat2(F, M)
                    g1d = np.einsum("pvqv->pq", gd) / (N - 1)
                    Md = mat(g1d, gd, affine=False) if cone != "D" else F / 2
                    rates[nm] = np.real(np.diag(U.conj().T @ (0.5 * (Md + Md.conj().T)) @ U))
                small = (x > 1e-9) & (x < 0.05 * x.max())
                xs = np.clip(x, 1e-300, None)
                dS = {nm: -np.sum(r * (np.log(xs) + 1)) for nm, r in rates.items()}
                lab = label if cone == "D" else ""
                if small.sum() == 0:
                    print("%-38s %s   none %60s %+.1e / %+.1e / %+.1e" % (lab, cone, "", dS["F0"], dS["R3"], dS["ex"]))
                    continue
                f0, r3 = rates["F0"][small], rates["R3"][small]
                print("%-38s %s   %3d     %+.2e          %.2f                 %.2f          %+.1e / %+.1e / %+.1e"
                      % (lab, cone, small.sum(), np.mean(np.sign(f0) * r3), np.mean(np.sign(f0) * np.sign(r3) < 0),
                         np.mean(r3 > 0), dS["F0"], dS["R3"], dS["ex"]))


if __name__ == "__main__":
    main()
