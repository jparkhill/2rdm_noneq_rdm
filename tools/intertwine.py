#!/usr/bin/env python3
r"""Route 2: can the ladder be fixed DYNAMICALLY, by the intertwining relation?

THE RELATION. With V^{(k)} the span of the k-fold annihilations of |Psi>, P_k the isometry onto it,
frak_a_x^{(k)} = P_{k+1}^dag a_x P_k the ladder, and K^{(j)} = P_j^dag H P_j the compressed
Hamiltonian, the proposal is to determine the ladder from

    frak_a_x K^(k)  -  K^(k+1) frak_a_x   =   P_{k+1}^dag [a_x, H] P_k

which is linear in the ladder, so it fixes it given the compressed Hamiltonians.

IT IS NOT EXACT, AND THE ERROR HAS A CLOSED FORM. Writing Pi_j = P_j P_j^dag,

    frak_a_x K^(k)     = P_{k+1}^dag a_x Pi_k H P_k
    K^(k+1) frak_a_x   = P_{k+1}^dag H Pi_{k+1} a_x P_k ,  and  P^dag (P P^dag - I) = 0

so the second term is EXACT and the whole residual comes from the first:

    R = frak_a_x K^(k) - K^(k+1) frak_a_x - P_{k+1}^dag [a_x,H] P_k
      = - P_{k+1}^dag a_x (I - Pi_k) H P_k

the leakage of H out of V^{(k)}, seen through a_x. If that is small the relation is a usable
dynamical closure for the ladder; if it is the size of the relation itself, route 2 fails for the
same kind of reason route 1 did.

WHAT IS REPORTED, per level k, over all M orbitals: the norm of the relation's right-hand side, the
norm of the residual, and their ratio. A NON-STATIONARY state is used -- at an eigenstate H P_0 is
parallel to P_0 and the k = 0 leakage would be spuriously zero.

Ground truth qdiss.eom.fullci.
"""
import sys
from math import factorial
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from molecular_2rdm import build_H_sector, molecule, sector
from qdiss import ansatz as az
from qdiss.eom import fullci

CASES = (("H 0.00 0.00 0.0; H 1.00 0.00 0.0; H 0.35 0.90 0.0", "sto-3g", 3, 1, "H3 bent"),
         ("H 0.00 0.00 0.0; H 0.80 0.00 0.0; H 1.05 0.72 0.0; H -0.15 0.72 0.0",
          "sto-3g", 4, 0, "H4 trapezoid"))


def run(atom, basis, N, spin, label):
    mol, mf, M, h, V, enuc = molecule(atom, basis, 0, spin)
    dets = sector(M, N)
    es, vec = np.linalg.eigh(build_H_sector(h, V, M, dets))
    rng = np.random.default_rng(7)
    c = vec[:, 0] + 0.4 * (rng.normal(size=len(dets)) + 1j * rng.normal(size=len(dets)))
    c /= np.linalg.norm(c)
    psi = np.zeros(1 << M, dtype=complex)
    psi[dets] = c
    H = fullci.build_H(h, V, M)

    def annihilate(p, tup):
        out = p
        for t in tup:
            out = fullci.aa(t, out)
        return out

    # the tower
    P, SEC = {}, {}
    for j in range(N + 1):
        SEC[j] = sector(M, N - j)
        TUP = az.tuples(M, j)
        U = np.zeros((len(SEC[j]), max(len(TUP), 1)), dtype=complex)
        for a, I in enumerate(TUP):
            U[:, a] = annihilate(psi, I)[SEC[j]]
        Q, sv, _ = np.linalg.svd(np.sqrt(float(factorial(j))) * U, full_matrices=False)
        P[j] = Q[:, sv > 1e-10 * max(sv[0], 1e-30)]

    def embed(col, j):
        b = np.zeros(1 << M, dtype=complex)
        b[SEC[j]] = col
        return b

    print("=== %s  M=%d N=%d   tower dims %s ==="
          % (label, M, N, [P[j].shape[1] for j in range(N + 1)]))
    print("    %-6s %12s %12s %12s   %s"
          % ("level k", "|RHS|", "|residual|", "ratio", "verdict"))
    for k in range(N):
        rk, rk1 = P[k].shape[1], P[k + 1].shape[1]
        Kk = np.zeros((rk, rk), dtype=complex)
        for al in range(rk):
            Hv = H @ embed(P[k][:, al], k)
            Kk[:, al] = P[k].conj().T @ Hv[SEC[k]]
        Kk1 = np.zeros((rk1, rk1), dtype=complex)
        for al in range(rk1):
            Hv = H @ embed(P[k + 1][:, al], k + 1)
            Kk1[:, al] = P[k + 1].conj().T @ Hv[SEC[k + 1]]
        lad = np.zeros((M, rk1, rk), dtype=complex)
        for x in range(M):
            for al in range(rk):
                lad[x, :, al] = P[k + 1].conj().T @ fullci.aa(x, embed(P[k][:, al], k))[SEC[k + 1]]
        nrhs = nres = 0.0
        for x in range(M):
            rhs = np.zeros((rk1, rk), dtype=complex)
            for al in range(rk):
                v = embed(P[k][:, al], k)
                comm = fullci.aa(x, H @ v) - H @ fullci.aa(x, v)
                rhs[:, al] = P[k + 1].conj().T @ comm[SEC[k + 1]]
            res = lad[x] @ Kk - Kk1 @ lad[x] - rhs
            nrhs += float(np.linalg.norm(rhs)) ** 2
            nres += float(np.linalg.norm(res)) ** 2
        nrhs, nres = np.sqrt(nrhs), np.sqrt(nres)
        ratio = nres / max(nrhs, 1e-30)
        print("    %-6d %12.5f %12.5f %12.4f   %s"
              % (k, nrhs, nres, ratio,
                 "usable" if ratio < 0.1 else ("marginal" if ratio < 0.5 else "FAILS")))
    print()


for a, bs, n, s, l in (CASES if __name__ == "__main__" else ()):
    run(a, bs, n, s, l)
