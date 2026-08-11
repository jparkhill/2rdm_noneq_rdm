#!/usr/bin/env python3
"""Cast the exact case into the XR form: does W exist, and what is it made of?

THE CLAIM TO TEST. In the XR construction X = Mat_k = R^dag R. The k-annihilated states
    u_I = a_{i_k} ... a_{i_1} |Psi>,   I = (i_1 < ... < i_k)
live in the (N-k)-particle sector, and Mat_k = k! * Gram(u). So R is (proportional to) the
coefficient matrix of the u_I in a basis of that sector.

The (k+1)-annihilated states are ONE MORE annihilator applied to those:
    w_L = a_{l_{k+1}} u_{(l_1..l_k)}       for a sorted (k+1)-tuple L
so Mat_{k+1} = (k+1)! * Gram(w) = W^dag W with W built from R by applying the single-orbital
annihilators. If that holds, then:

  (a) POSITIVITY of gamma_{k+1} is structural, since Mat_{k+1} = W^dag W.
  (b) The CONTRACTION IDENTITY is automatic: the partial trace of W^dag W over the last index
      contracts sum_x a_x^dag a_x = Nhat, which acts as (N-k) on the (N-k)-particle sector, giving
      Tr_{k+1} gamma_{k+1} = (N-k) gamma_k exactly -- and that is the moment form of the
      single-trace cumulant identity.
  (c) EXACTNESS AT THE TOP follows from (b): at N = k the factor (N-k) is zero, so
      sum_x A_x^dag A_x = 0, so every A_x kills the space, so W = 0 and gamma_{k+1} = 0 exactly.

WHAT IS MEASURED, on H3 bent / sto-3g, M = 6, N = 3, at a non-stationary state, with k = 2 so that
the reconstructed object is gamma_3:
  1. Mat_2 == 2 * Gram(u), fixing the normalisation of R.
  2. Mat_3 == 6 * Gram(w), fixing the normalisation of W.
  3. W built ONLY from R and the known matrices of a_x on the (N-k)-particle sector reproduces
     Mat_3 exactly. This is the constructive claim.
  4. sum_x A_x^dag A_x == (N-k) * identity on that sector.
  5. THE GAUGE OBSTRUCTION. X = R^dag R fixes R only up to R -> U R with U unitary on the
     (N-k)-particle sector. Applying a random such U and rebuilding W tests whether Mat_{k+1} is
     determined by X alone. If it is not, the exact W needs information the reduced description
     does not carry, and that gap IS the closure problem.

Symbols: gamma_k[p_1..p_k, q_1..q_k] = <a+_p1 .. a+_pk a_qk .. a_q1>;
Mat_k(T)[I,J] = k! T[I+J] over sorted k-tuples. Ground truth qdiss.eom.fullci.
"""
import itertools
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from molecular_2rdm import build_H_sector, molecule, sector
from qdiss import ansatz as az
from qdiss.eom import fullci

H3_BENT = "H 0.00 0.00 0.0; H 1.00 0.00 0.0; H 0.35 0.90 0.0"
N, K = 3, 2

mol, mf, M, h, V, enuc = molecule(H3_BENT, "sto-3g", 0, 1)
dets = sector(M, N)
es, vec = np.linalg.eigh(build_H_sector(h, V, M, dets))
rng = np.random.default_rng(7)
c = vec[:, 0] + 0.4 * (rng.normal(size=len(dets)) + 1j * rng.normal(size=len(dets)))
c /= np.linalg.norm(c)
psi = np.zeros(1 << M, dtype=complex)
psi[dets] = c

PAIR = az.tuples(M, K)
TRIP = az.tuples(M, K + 1)
SEC_LO = sector(M, N - K)          # (N-k)-particle sector: where the u_I live
SEC_LO2 = sector(M, N - K - 1)     # one particle fewer: where the w_L live
print("M=%d N=%d k=%d   %d pairs, %d triples, (N-k)-sector dim %d, (N-k-1)-sector dim %d"
      % (M, N, K, len(PAIR), len(TRIP), len(SEC_LO), len(SEC_LO2)))


def annihilate_tuple(p, tup):
    """a_{t_last} ... a_{t_first} applied to a Fock vector, in ascending-index order."""
    out = p
    for t in tup:
        out = fullci.aa(t, out)
    return out


# ---- 1. Mat_2 = 2 * Gram(u), giving R
U = np.zeros((len(SEC_LO), len(PAIR)), dtype=complex)
for a, I in enumerate(PAIR):
    v = annihilate_tuple(psi, I)
    U[:, a] = v[SEC_LO]
X = az.tensor_to_coord(fullci.rdm(psi, M, K), M, K) * 2.0
print("1. |Mat_2 - 2 Gram(u)|            = %.3e   (|Mat_2| = %.5f)"
      % (np.abs(X - 2.0 * (U.conj().T @ U)).max(), np.abs(X).max()))
R = np.sqrt(2.0) * U                       # so that X = R^dag R

# ---- 2. Mat_3 = 6 * Gram(w)
Wv = np.zeros((len(SEC_LO2), len(TRIP)), dtype=complex)
for a, L in enumerate(TRIP):
    v = annihilate_tuple(psi, L)
    Wv[:, a] = v[SEC_LO2]
Y = az.tensor_to_coord(fullci.rdm(psi, M, K + 1), M, K + 1) * 6.0
print("2. |Mat_3 - 6 Gram(w)|            = %.3e   (|Mat_3| = %.5f)"
      % (np.abs(Y - 6.0 * (Wv.conj().T @ Wv)).max(), np.abs(Y).max()))

# ---- the matrices of a_x between the two sectors
A = {}
for x in range(M):
    Ax = np.zeros((len(SEC_LO2), len(SEC_LO)), dtype=complex)
    for j, d in enumerate(SEC_LO):
        b = np.zeros(1 << M, dtype=complex)
        b[d] = 1.0
        Ax[:, j] = fullci.aa(x, b)[SEC_LO2]
    A[x] = Ax

# ---- 3. W from R alone: w_L = a_{l_{k+1}} u_{(l_1..l_k)}, last index appended
pos = {t: i for i, t in enumerate(PAIR)}
W_from_R = np.zeros((len(SEC_LO2), len(TRIP)), dtype=complex)
for a, L in enumerate(TRIP):
    I, x = L[:K], L[K]
    W_from_R[:, a] = A[x] @ (R[:, pos[I]] / np.sqrt(2.0))
print("3. |Mat_3 - 6 W^dag W| with W from R and the known a_x  = %.3e"
      % np.abs(Y - 6.0 * (W_from_R.conj().T @ W_from_R)).max())

# ---- 4. sum_x A_x^dag A_x = (N-k) * identity
S = sum(A[x].conj().T @ A[x] for x in range(M))
print("4. |sum_x A_x^dag A_x - (N-k) I| = %.3e   (N-k = %d, sector dim %d)"
      % (np.abs(S - (N - K) * np.eye(len(SEC_LO))).max(), N - K, len(SEC_LO)))

# ---- 5. the gauge obstruction
G = rng.normal(size=(len(SEC_LO),) * 2) + 1j * rng.normal(size=(len(SEC_LO),) * 2)
Q, _ = np.linalg.qr(G)
Rg = Q @ R
print()
print("5. GAUGE. X = R^dag R fixes R only up to R -> U R with U unitary on the (N-k)-sector.")
print("   |X - (U R)^dag (U R)|         = %.3e   (X is unchanged, as it must be)"
      % np.abs(X - Rg.conj().T @ Rg).max())
Wg = np.zeros((len(SEC_LO2), len(TRIP)), dtype=complex)
for a, L in enumerate(TRIP):
    I, x = L[:K], L[K]
    Wg[:, a] = A[x] @ (Rg[:, pos[I]] / np.sqrt(2.0))
Yg = 6.0 * (Wg.conj().T @ Wg)
print("   |Mat_3(gauged) - Mat_3|       = %.3e   (|Mat_3| = %.5f)"
      % (np.abs(Yg - Y).max(), np.abs(Y).max()))
print("   -> if this is NOT zero, Mat_{k+1} is not determined by X alone, and the missing")
print("      information is exactly the closure problem.")
