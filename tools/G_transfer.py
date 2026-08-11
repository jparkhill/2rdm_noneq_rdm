#!/usr/bin/env python3
"""The G-transfer form: gamma_{k+1} from R_c and a small PSD block matrix, all polynomial.

WHERE THIS COMES FROM. tools/W_construction.py showed Mat_{k+1} = W^dag W with W built from R by
the single-orbital annihilators, but that Mat_{k+1} is NOT determined by X = R^dag R: a gauge
rotation R -> U R left X unchanged at 4.4e-16 while moving Mat_3 by 1.444 against its own norm of
1.362. The gauge lives in an isometry P onto range(R) inside the (N-k)-particle sector, which is
exponentially tall. This file shows P enters ONLY through small blocks.

DERIVATION. With u_I = a_{i_k}..a_{i_1}|Psi> the k-annihilated states and Mat_k = k! Gram(u), write
R = sqrt(k!) * (coefficients of u) and factor R = P R_c with P an isometry onto range(R), r = rank R,
so R_c is r x C(M,k) and X = R^dag R = R_c^dag R_c. The (k+1)-annihilated states are
w_L = a_{l_{k+1}} u_{(l_1..l_k)} with L sorted, so l_{k+1} is the largest element, and

    Mat_{k+1}[(I,x),(J,y)]  =  (k+1) * R_c(:,I)^dag  G_xy  R_c(:,J),
    G_xy  =  P^dag a_x^dag a_y P                                  (an r x r block)

So the ONLY thing needed beyond R_c is the M x M array of r x r blocks G_xy. Assembled as an
(M r) x (M r) matrix, G is the Gram matrix of the vectors {a_x e_alpha} for {e_alpha} a basis of
range(R), hence

    G >= 0                             positivity of gamma_{k+1} is STRUCTURAL
    sum_x G_xx = (N-k) I_r             because the e_alpha lie in the (N-k)-particle sector, and
                                       sum_x a_x^dag a_x = Nhat. This is the contraction identity
                                       Tr_{k+1} gamma_{k+1} = (N-k) gamma_k, i.e. the moment form
                                       of the single-trace cumulant identity.
    at N = k the prefactor is zero, so G >= 0 with zero trace forces G = 0, hence gamma_{k+1} = 0:
                                       EXACTNESS AT THE TOP, automatic.

SIZES. r = rank(Mat_k) <= C(M,k), so G is at most M*C(M,k) square -- polynomial in M at fixed k.
The exponentially large (N-k)-particle sector never appears.

WHAT IS VERIFIED HERE, on exact FCI states at a NON-STATIONARY point (an eigenstate would make the
dynamical comparisons vacuous):
  1. R = P R_c reproduces X exactly.
  2. Mat_{k+1} equals the G-transfer expression exactly.
  3. G >= 0 as an (M r) x (M r) matrix.
  4. sum_x G_xx = (N-k) I_r.
  5. sizes, to confirm nothing exponential is required.

Ground truth qdiss.eom.fullci. Symbols: gamma_k[p_1..p_k, q_1..q_k] = <a+_p1..a+_pk a_qk..a_q1>,
Mat_k(T)[I,J] = k! T[I+J] over sorted k-tuples.
"""
import sys
from math import comb, factorial
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from molecular_2rdm import build_H_sector, molecule, sector
from qdiss import ansatz as az
from qdiss.eom import fullci

CASES = (("H 0.00 0.00 0.0; H 1.00 0.00 0.0; H 0.35 0.90 0.0", 3, 1, 2, "H3 bent   k=2"),
         ("H 0.00 0.00 0.0; H 0.80 0.00 0.0; H 1.05 0.72 0.0; H -0.15 0.72 0.0",
          4, 0, 2, "H4 trap   k=2"),
         ("H 0.00 0.00 0.0; H 1.00 0.00 0.0; H 0.35 0.90 0.0", 3, 1, 1, "H3 bent   k=1"))


def annihilate(p, tup):
    out = p
    for t in tup:
        out = fullci.aa(t, out)
    return out


for atom, N, spin, K in [(a, n, s, k) for a, n, s, k, _ in CASES]:
    label = [lab for a, n, s, k, lab in CASES if (a, n, s, k) == (atom, N, spin, K)][0]
    mol, mf, M, h, V, enuc = molecule(atom, "sto-3g", 0, spin)
    dets = sector(M, N)
    es, vec = np.linalg.eigh(build_H_sector(h, V, M, dets))
    rng = np.random.default_rng(7)
    c = vec[:, 0] + 0.4 * (rng.normal(size=len(dets)) + 1j * rng.normal(size=len(dets)))
    c /= np.linalg.norm(c)
    psi = np.zeros(1 << M, dtype=complex)
    psi[dets] = c

    TUP_K = az.tuples(M, K)
    TUP_K1 = az.tuples(M, K + 1)
    SEC_HI = sector(M, N - K)
    SEC_LO = sector(M, N - K - 1)

    # u_I coefficients in the (N-k)-particle sector, and R
    U = np.zeros((len(SEC_HI), len(TUP_K)), dtype=complex)
    for a, I in enumerate(TUP_K):
        U[:, a] = annihilate(psi, I)[SEC_HI]
    R = np.sqrt(float(factorial(K))) * U
    X = az.tensor_to_coord(fullci.rdm(psi, M, K), M, K) * factorial(K)

    # R = P R_c with P an isometry onto range(R)
    Q, S, Vh = np.linalg.svd(R, full_matrices=False)
    keep = S > 1e-10 * max(S[0], 1e-30)
    r = int(keep.sum())
    P = Q[:, keep]
    Rc = (np.diag(S[keep]) @ Vh[keep, :])

    # G_xy = (a_x P)^dag (a_y P), assembled as an (M r) x (M r) Gram matrix
    aP = np.zeros((M, len(SEC_LO), r), dtype=complex)
    for x in range(M):
        for al in range(r):
            b = np.zeros(1 << M, dtype=complex)
            b[SEC_HI] = P[:, al]
            aP[x, :, al] = fullci.aa(x, b)[SEC_LO]
    Gbig = np.zeros((M * r, M * r), dtype=complex)
    for x in range(M):
        for y in range(M):
            Gbig[x * r:(x + 1) * r, y * r:(y + 1) * r] = aP[x].conj().T @ aP[y]

    # the G-transfer expression for Mat_{k+1}
    Yex = az.tensor_to_coord(fullci.rdm(psi, M, K + 1), M, K + 1) * factorial(K + 1)
    posK = {t: i for i, t in enumerate(TUP_K)}
    Yg = np.zeros_like(Yex)
    for a, L in enumerate(TUP_K1):
        I, x = L[:K], L[K]
        for b, Lp in enumerate(TUP_K1):
            J, y = Lp[:K], Lp[K]
            Gxy = Gbig[x * r:(x + 1) * r, y * r:(y + 1) * r]
            Yg[a, b] = (K + 1) * (Rc[:, posK[I]].conj() @ (Gxy @ Rc[:, posK[J]]))

    ev = np.linalg.eigvalsh(0.5 * (Gbig + Gbig.conj().T))
    Ssum = sum(Gbig[x * r:(x + 1) * r, x * r:(x + 1) * r] for x in range(M))
    print("=== %s   M=%d N=%d   (N-k)-sector dim %d, (N-k-1) dim %d ==="
          % (label, M, N, len(SEC_HI), len(SEC_LO)))
    print("    1. |X - R_c^dag R_c|                  = %.3e   (|X| = %.5f)"
          % (np.abs(X - Rc.conj().T @ Rc).max(), np.abs(X).max()))
    print("    2. |Mat_%d - G-transfer expression|    = %.3e   (|Mat_%d| = %.5f)"
          % (K + 1, np.abs(Yex - Yg).max(), K + 1, np.abs(Yex).max()))
    print("    3. min eig G (as %d x %d Gram)         = %.3e" % (M * r, M * r, ev.min()))
    print("    4. |sum_x G_xx - (N-k) I_r|           = %.3e   (N-k = %d)"
          % (np.abs(Ssum - (N - K) * np.eye(r)).max(), N - K))
    print("    5. sizes: r = rank(Mat_%d) = %d of %d,  G is %d x %d,  C(M,k) = %d"
          % (K, r, len(TUP_K), M * r, M * r, comb(M, K)))
    print()
