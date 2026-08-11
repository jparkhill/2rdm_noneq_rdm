#!/usr/bin/env python3
"""The computable hierarchy: never store R, generate it by ladder products.

THE CLAIM THAT MAKES IT COMPUTABLE. In the G-transfer form the state at level k is
R_c^{(k)}, of size r_k x C(M,k) -- and C(M,k) is exponential at k ~ N/2, so storing it for every
level would be hopeless. But it need not be stored. Since u_I = a_{i_k}..a_{i_1}|Psi> and each
annihilator is the compressed ladder,

    R_c^{(k)}(:, I)  =  frak_a_{i_k}^{(k-1)} ... frak_a_{i_1}^{(0)} v_0

where v_0 spans the one-dimensional V^{(0)} = span{|Psi>}. So the ONLY stored data is the tower of
ladders {frak_a_x^{(j)}} for j = 0..N-1, and every k-tuple column is a matrix product along the
tuple. That is a matrix-product factorisation in the PARTICLE-NUMBER direction, with the orbital
label as the physical index on each level.

CONSEQUENCE. The convergence parameter is not the rank k but the compressed dimension chi_j. Capping
chi_j <= chi gives storage O(N M chi^2) and makes the scheme exact as chi -> the true ranks
r_j = rank(Mat_j) <= min(C(M,j), C(M,N-j)). Rank k only says which reduced matrix you look at; chi
says how accurate the whole tower is.

WHAT IS VERIFIED HERE, on H3 bent / sto-3g, M = 6, N = 3, at a non-stationary state:
  1. the tower dimensions r_j = rank(Mat_j) against the bound min(C(M,j), C(M,N-j));
  2. that ladder products reproduce R_c^{(k)} exactly, for every k -- the load-bearing claim;
  3. that Mat_k assembled from those products equals the exact Mat_k;
  4. Sum_x frak_a_x^{(j)dag} frak_a_x^{(j)} = (N-j) I at every level;
  5. the anticommutation {frak_a_x^{(j+1)}, frak_a_y^{(j)}} = 0 across levels, which is the CAR
     tower-consistency condition;
  6. storage of the tower against the wavefunction, so the crossover is stated honestly rather than
     assumed.

Symbols: Mat_k(T)[I,J] = k! T[I+J] over sorted k-tuples, T = gamma_k. Ground truth qdiss.eom.fullci.
"""
import sys
from math import comb, factorial
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from molecular_2rdm import build_H_sector, molecule, sector
from qdiss import ansatz as az
from qdiss.eom import fullci

H3_BENT = "H 0.00 0.00 0.0; H 1.00 0.00 0.0; H 0.35 0.90 0.0"
N = 3

mol, mf, M, h, V, enuc = molecule(H3_BENT, "sto-3g", 0, 1)
dets = sector(M, N)
es, vec = np.linalg.eigh(build_H_sector(h, V, M, dets))
rng = np.random.default_rng(7)
c = vec[:, 0] + 0.4 * (rng.normal(size=len(dets)) + 1j * rng.normal(size=len(dets)))
c /= np.linalg.norm(c)
psi = np.zeros(1 << M, dtype=complex)
psi[dets] = c


def annihilate(p, tup):
    out = p
    for t in tup:
        out = fullci.aa(t, out)
    return out


# ---- build the tower: for each level j, an isometry P_j onto V^{(j)} inside the (N-j)-sector
P, SEC, Rc = {}, {}, {}
for j in range(N + 1):
    SEC[j] = sector(M, N - j)
    TUP = az.tuples(M, j)
    U = np.zeros((len(SEC[j]), max(len(TUP), 1)), dtype=complex)
    for a, I in enumerate(TUP):
        U[:, a] = annihilate(psi, I)[SEC[j]]
    Rfull = np.sqrt(float(factorial(j))) * U
    Q, S, Vh = np.linalg.svd(Rfull, full_matrices=False)
    keep = S > 1e-10 * max(S[0], 1e-30)
    P[j] = Q[:, keep]
    Rc[j] = np.diag(S[keep]) @ Vh[keep, :]

print("=== ladder tower, H3 bent / sto-3g, M=%d N=%d ===" % (M, N))
print("1. tower dimensions")
print("   %-6s %8s %10s %14s" % ("level j", "r_j", "C(M,j)", "min(C(M,j),C(M,N-j))"))
for j in range(N + 1):
    print("   %-6d %8d %10d %14d"
          % (j, P[j].shape[1], comb(M, j), min(comb(M, j), comb(M, N - j))))

# ---- the ladders
lad = {}
for j in range(N):
    A = np.zeros((M, P[j + 1].shape[1], P[j].shape[1]), dtype=complex)
    for x in range(M):
        for al in range(P[j].shape[1]):
            b = np.zeros(1 << M, dtype=complex)
            b[SEC[j]] = P[j][:, al]
            A[x, :, al] = P[j + 1].conj().T @ fullci.aa(x, b)[SEC[j + 1]]
    lad[j] = A

# ---- 2 & 3: ladder products reproduce R_c and Mat_k
print()
print("2/3. ladder products against the exact R_c and Mat_k")
print("     (the product gives raw u_I coefficients; R_c carries sqrt(k!) from Mat_k = k! Gram(u))")
print("   %-4s %14s %16s" % ("k", "|Rc - product|", "|Mat_k - assembled|"))
v0 = Rc[0][:, 0]                       # the one-dimensional top state
for k in range(1, N + 1):
    TUP = az.tuples(M, k)
    Rp = np.zeros((P[k].shape[1], len(TUP)), dtype=complex)
    for a, I in enumerate(TUP):
        col = v0.copy()
        for lev, x in enumerate(I):
            col = lad[lev][x] @ col
        Rp[:, a] = col
    # The ladder product gives the raw coefficients of u_I; R_c carries the sqrt(k!) of
    # Mat_k = k! Gram(u), so the two differ by exactly that factor.
    Rp = np.sqrt(float(factorial(k))) * Rp
    Xex = az.tensor_to_coord(fullci.rdm(psi, M, k), M, k) * factorial(k)
    print("   %-4d %14.3e %16.3e"
          % (k, np.abs(Rp - Rc[k]).max(), np.abs(Rp.conj().T @ Rp - Xex).max()))

# ---- 4: the number condition at every level
print()
print("4. number condition  sum_x a_x^dag a_x = (N-j) I  at each level")
for j in range(N):
    S = sum(lad[j][x].conj().T @ lad[j][x] for x in range(M))
    print("   level %d: |sum - (N-j)I| = %.3e   (N-j = %d, dim %d)"
          % (j, np.abs(S - (N - j) * np.eye(P[j].shape[1])).max(), N - j, P[j].shape[1]))

# ---- 5: CAR tower consistency across levels
print()
print("5. CAR tower consistency  {a_x^{(j+1)}, a_y^{(j)}} = 0")
for j in range(N - 1):
    worst = 0.0
    for x in range(M):
        for y in range(M):
            Z = lad[j + 1][x] @ lad[j][y] + lad[j + 1][y] @ lad[j][x]
            worst = max(worst, float(np.abs(Z).max()))
    print("   levels %d->%d: max |anticommutator| = %.3e" % (j, j + 2, worst))

# ---- 6: storage, honestly
print()
tower = sum(M * lad[j].shape[1] * lad[j].shape[2] for j in range(N))
print("6. storage: tower %d complex numbers, wavefunction C(M,N) = %d."
      % (tower, comb(M, N)))
print("   At this size the tower is LARGER; it pays only when the true r_j are far below")
print("   min(C(M,j),C(M,N-j)), or when chi_j is capped below r_j. The exact tower is")
print("   exponential, since r_j peaks near C(M,N/2) -- capping chi is the approximation.")
