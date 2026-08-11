#!/usr/bin/env python3
"""How does the semidefinite CP constraint scale with M? Measure the Choi rank.

THE QUESTION. The positive-generator constraint is Kossakowski = Q Choi(L) Q >= 0 with
    Choi(L) = sum_ij |i><j| (x) L[|i><j|],   Q = 1 - |Omega><Omega|,  |Omega> = n^-1/2 sum_i |i>|i>
on Mat_k of dimension n_k = C(M,k). Dense, that is an n_k^2 x n_k^2 eigenproblem: O(M^{4k}) memory
and O(M^{6k}) time, which dies at about M = 16 for rank 2 and M = 10 for rank 3.

THE HOPED-FOR COLLAPSE. For a single sandwich term L[X] = A X B,
    Choi[(i,k),(j,l)] = A_ki B_jl
which is RANK ONE. Every equation-of-motion term is a sandwich, so Choi should have rank at most
T = (number of sandwich terms) = 2k + 2k(k-1) + 4k plus one per closure structure -- O(k^2), and
INDEPENDENT OF M. If so, lambda_min follows from a T x T Gram matrix and the cost falls to
O(M^{2k}), cheaper than one rank-2 EOM evaluation.

THE CAVEAT BEING TESTED. Not every term is a clean sandwich. gamma_1 = Tr_2(gamma_2)/(N-1) is
recovered by contraction at every stage, and the trace-identity dressing calls Tr_k. A partial
trace is not a sandwich, and it may raise the Choi rank by O(1) or by O(n_{k-1}). Which of those
happens decides between the cheap and the expensive regime. So: measure the numerical rank of
Choi(L) for each generator, at two values of M, and see whether it grows.

Symbols and conventions as in tools/exactness_limit.py, whose machinery this reuses and which
validated it: the exact rank-N generator there gave Kossakowski +-1.1e-15 and the rank-2 drift
+-1.5e-15.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from ansatz_hierarchy import rank3_set
from molecular_2rdm import build_H_sector, molecule, sector
from qdiss import ansatz as az
from qdiss.eom import fullci, wick

SYS = (("H 0.00 0.00 0.0; H 1.00 0.00 0.0; H 0.35 0.90 0.0", "sto-3g", 3, 0, 1, "H3 bent"),
       ("H 0.00 0.00 0.0; H 0.80 0.00 0.0; H 1.05 0.72 0.0; H -0.15 0.72 0.0",
        "sto-3g", 4, 0, 0, "H4 trapezoid"))


def choi_of(L, n):
    C = np.zeros((n * n, n * n), dtype=complex)
    for i in range(n):
        for j in range(n):
            LE = L[:, i * n + j].reshape(n, n)
            C[i * n:i * n + n, j * n:j * n + n] = LE
    return C


def kossakowski_from_choi(C, n):
    om = np.zeros(n * n, dtype=complex)
    for i in range(n):
        om[i * n + i] = 1.0 / np.sqrt(n)
    Q = np.eye(n * n, dtype=complex) - np.outer(om, om.conj())
    K = Q @ C @ Q
    return np.linalg.eigvalsh(0.5 * (K + K.conj().T))


print("Choi rank of the rank-2 generator. If it is independent of M, the CP constraint reduces")
print("to a T x T Gram problem and the cost is O(M^4) rather than O(M^12).")
print()
print("%-14s %3s %5s %7s | %-14s %8s %8s %10s"
      % ("system", "M", "n_2", "dim", "generator", "choiRank", "kossRank", "min eig"))
for atom, basis, N, chg, spin, label in SYS:
    mol, mf, M, h, V, enuc = molecule(atom, basis, chg, spin)
    dets = sector(M, N)
    es, vec = np.linalg.eigh(build_H_sector(h, V, M, dets))
    rng = np.random.default_rng(7)
    c = vec[:, 0] + 0.4 * (rng.normal(size=len(dets)) + 1j * rng.normal(size=len(dets)))
    c /= np.linalg.norm(c)
    psi = np.zeros(1 << M, dtype=complex)
    psi[dets] = c
    g2 = fullci.rdm(psi, M, 2)
    g1 = fullci.rdm(psi, M, 1)
    n2 = len(az.tuples(M, 2))
    EYE = np.eye(M, dtype=complex)
    Z3 = np.zeros((M,) * 6, dtype=complex)

    def mat2(T):
        return 2.0 * az.tensor_to_coord(T, M, 2)

    def unmat2(A):
        return az.coord_to_tensor(A / 2.0, M, 2)

    def drift(T):
        return wick.gamma_dot(2, h, V, {("D", 1): EYE, ("G", 1): g1, ("G", 2): T,
                                        ("G", 3): Z3}, M)

    def closed(T, nm):
        gg1 = np.einsum("pxqx->pq", T) / (N - 1)
        lam2 = az.cumulant(2, T, {("G", 1): gg1})
        g3 = az.separable(3, {("G", 1): gg1, ("L", 2): lam2}) + rank3_set(gg1, lam2, N, M)[nm]
        return wick.gamma_dot(2, h, V, {("D", 1): EYE, ("G", 1): gg1, ("G", 2): T,
                                        ("G", 3): g3}, M)

    gens = [("drift only", lambda T: drift(T), True)]
    for nm in ("null", "trace", "scaled+trace"):
        gens.append((nm, (lambda T, k=nm: closed(T, k)), False))

    X0 = mat2(g2)
    eps = 1e-5
    for nm, flow, linear in gens:
        L = np.zeros((n2 * n2, n2 * n2), dtype=complex)
        f0 = None if linear else mat2(flow(unmat2(X0)))
        for i in range(n2):
            for j in range(n2):
                E = np.zeros((n2, n2), dtype=complex)
                E[i, j] = 1.0
                if linear:
                    col = mat2(flow(unmat2(E)))
                else:
                    col = (mat2(flow(unmat2(X0 + eps * E))) - f0) / eps
                L[:, i * n2 + j] = col.reshape(-1)
        C = choi_of(L, n2)
        sv = np.linalg.svd(C, compute_uv=False)
        cr = int((sv > 1e-8 * sv[0]).sum())
        ev = kossakowski_from_choi(C, n2)
        kr = int((np.abs(ev) > 1e-8 * max(np.abs(ev).max(), 1e-30)).sum())
        print("%-14s %3d %5d %7d | %-14s %8d %8d %10.3e"
              % (label, M, n2, n2 * n2, nm, cr, kr, ev.min()))
    print()
