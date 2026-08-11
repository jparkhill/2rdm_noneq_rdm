#!/usr/bin/env python3
"""What rank does the true gamma_j actually need? The singular spectrum of R^(j).

WHY THIS IS THE RIGHT QUESTION FOR A CLOSURE. In the one-rung BBGKY recipe (docs/G-TRANSFER.md)
the closure parameter is chi', the target dimension of the single modelled ladder
frak_a_x^(k) : C^chi -> C^chi'. Since Mat_{k+1} = W^dag W with W of chi' rows, chi' IS THE RANK WE
ALLOW THE RECONSTRUCTED gamma_{k+1} TO HAVE. So the question that decides whether a small chi' is
admissible is: what effective rank does the EXACT gamma_{k+1} carry?

That is measurable a priori. With u_I = a_{i_j}..a_{i_1}|Psi> and Mat_j = j! Gram(u), the singular
values of R^(j) = sqrt(j!) * (coefficients of u) give the spectrum of Mat_j, and

    chi(f) = smallest chi whose top-chi squared singular values reach fraction f of the total.

Bounds for context: r_j = rank Mat_j <= min(C(M,j), C(M,N-j)), the first because there are only
C(M,j) spanning vectors, the second because they all lie in the (N-j)-particle sector. The second
factor only bites when N-j is small, i.e. k close to N.

Ground truth qdiss.eom.fullci at the FCI ground state.
"""
import sys
import numpy as np
from math import comb, factorial
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from molecular_2rdm import build_H_sector, molecule, sector
from qdiss import ansatz as az
from qdiss.eom import fullci


def ann(p, t):
    for x in t:
        p = fullci.aa(x, p)
    return p


CASES = (("H 0.00 0.00 0.0; H 0.80 0.00 0.0; H 1.05 0.72 0.0; H -0.15 0.72 0.0",
          "sto-3g", 4, 0, "H4 trap sto-3g"),
         ("Li 0 0 0; H 0 0 1.6", "sto-3g", 4, 0, "LiH sto-3g"),
         ("H 0 0 0; H 0 0 1.4; H 0 0 2.8; H 0 0 4.2; H 0 0 5.6; H 0 0 7.0",
          "sto-3g", 6, 0, "H6 chain sto-3g"))

print(__doc__.split("Ground truth")[0].strip()[:0] or "", end="")
print("Singular spectrum of R^(j): how many directions carry the weight?")
print("chi(f) = smallest chi whose top-chi squared singular values reach fraction f of the total.")
print()
print("%-17s %3s %3s | %2s %5s %6s | %6s %6s %8s %9s"
      % ("system", "M", "N", "j", "r_j", "bound", "chi90", "chi99", "chi99.9", "s_last/s_1"))
for atom, bas, N, spin, lab in (CASES if __name__ == "__main__" else ()):
    mol, mf, M, h, V, enuc = molecule(atom, bas, 0, spin)
    dets = sector(M, N)
    es, vec = np.linalg.eigh(build_H_sector(h, V, M, dets))
    psi = np.zeros(1 << M, dtype=complex)
    psi[dets] = vec[:, 0]
    for j in range(1, min(N, 3) + 1):
        HI, TUP = sector(M, N - j), az.tuples(M, j)
        U = np.zeros((len(HI), len(TUP)), dtype=complex)
        for a, I in enumerate(TUP):
            U[:, a] = ann(psi, I)[HI]
        S = np.linalg.svd(np.sqrt(float(factorial(j))) * U, compute_uv=False)
        w = S ** 2
        cw = np.cumsum(w / w.sum())
        rj = int((S > 1e-10 * S[0]).sum())
        f = lambda q: int(np.searchsorted(cw, q) + 1)
        print("%-17s %3d %3d | %2d %5d %6d | %6d %6d %8d %9.1e"
              % (lab, M, N, j, rj, min(comb(M, j), comb(M, N - j)),
                 f(.90), f(.99), f(.999), S[rj - 1] / S[0]))
