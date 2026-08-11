#!/usr/bin/env python3
"""At the top of its own hierarchy, how wrong is a cumulant truncation?

THE QUESTION. For an N-electron state gamma_{N+1} == 0 identically, so at k = N the rank-k
equation of motion closes with NO approximation. A hierarchy can only be CONVERGENT if it becomes
exact there. In cumulant language the exactness condition reads

    lambda_{N+1} = -separable(gamma_1, lambda_2, ..., lambda_N)      which is NOT zero

so a truncation lambda_{k+1} = 0 gives gamma_{N+1} = +separable instead of 0, and is wrong at the
top of its own hierarchy. This file measures by how much. If the number is large, the whole
lambda_k = 0 family is non-convergent and the exactness condition has to be built into any
candidate closure from the start.

WHAT IS COMPARED, at rank k = N = 3, so that gamma_4 is the reconstructed object:

    truth        qdiss.eom.fullci.rdm_dot(psi, H, M, 3), exact d/dt gamma_3 by linear algebra
    gamma_4 = 0  the EXACT closure at k = N. Must reproduce truth to round-off.
    lambda_4 = 0 the cumulant truncation, i.e. gamma_4 = separable(gamma_1, lambda_2, lambda_3),
                 the MKNO cumulant expansion of gamma_4 with the connected 4-cumulant deleted.

SYMBOLS. gamma_k[p_1..p_k, q_1..q_k] = <a+_p1 .. a+_pk a_qk .. a_q1>. lambda_k is the MKNO
k-cumulant. `separable` at rank 4 is every term of mkno.gamma_to_cumulants over 4 upper and 4 lower
indices that does NOT contain ("L",4); it is built from gamma_1, lambda_2 and lambda_3, all taken
EXACT here, so the only approximation in the third row is the deletion of lambda_4 itself.

A NON-STATIONARY state is used. An eigenstate gives d/dt gamma_3 = 0 and every comparison would be
zero against zero.

Ground truth: qdiss.eom.fullci. Term lists: qdiss.eom.wick, validated against fullci at 8.7e-16.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from molecular_2rdm import build_H_sector, molecule, sector
from qdiss import ansatz as az
from qdiss.eom import fullci, hamiltonian, wick
from qdiss.symbolic import mkno

UP4, LO4 = ("p", "v", "x", "z"), ("q", "w", "y", "t")
_ALL4 = mkno.gamma_to_cumulants(UP4, LO4)
SEP4 = [(s, ",".join("".join(u) + "".join(l) for n, r, u, l in f) + "->pvxzqwyt",
         [(n, r) for n, r, u, l in f])
        for s, f in _ALL4 if not any(n == "L" and r == 4 for n, r, u, l in f)]
print("MKNO expansion of gamma_4: %d terms total, %d with lambda_4 deleted"
      % (len(_ALL4), len(SEP4)))


def separable4(g1, l2, l3, M):
    tab = {("G", 1): g1, ("L", 2): l2, ("L", 3): l3}
    out = np.zeros((M,) * 8, dtype=complex)
    for s, sub, keys in SEP4:
        out += s * np.einsum(sub, *[tab[k] for k in keys], optimize=True)
    return out


H3_BENT = "H 0.00 0.00 0.0; H 1.00 0.00 0.0; H 0.35 0.90 0.0"
H3_LIN = "H 0 0 0; H 0 0 1.0; H 0 0 2.0"
N = 3

print()
print("%-22s | %11s %11s | %11s %11s %9s"
      % ("system", "|d g3/dt|", "|g4 exact|", "g4=0 resid", "lam4=0 resid", "rel err"))

CASES = [(H3_BENT, "H3 bent / sto-3g"), (H3_LIN, "H3 linear / sto-3g")]
for atom, label in CASES:
    mol, mf, M, h, V, enuc = molecule(atom, "sto-3g", 0, 1)
    dets = sector(M, N)
    es, vec = np.linalg.eigh(build_H_sector(h, V, M, dets))
    rng = np.random.default_rng(7)
    c = vec[:, 0] + 0.4 * (rng.normal(size=len(dets)) + 1j * rng.normal(size=len(dets)))
    c /= np.linalg.norm(c)
    psi = np.zeros(1 << M, dtype=complex)
    psi[dets] = c
    H = fullci.build_H(h, V, M)
    g = {k: fullci.rdm(psi, M, k) for k in (1, 2, 3, 4)}
    l2 = az.cumulant(2, g[2], {("G", 1): g[1]})
    l3 = az.cumulant(3, g[3], {("G", 1): g[1], ("L", 2): l2})
    truth = fullci.rdm_dot(psi, H, M, 3)
    tab = {("D", 1): np.eye(M, dtype=complex), ("G", 1): g[1], ("G", 2): g[2], ("G", 3): g[3]}
    tab[("G", 4)] = np.zeros((M,) * 8, dtype=complex)
    r_exact = np.abs(wick.gamma_dot(3, h, V, tab, M) - truth).max()
    sep4 = separable4(g[1], l2, l3, M)
    tab[("G", 4)] = sep4
    r_null = np.abs(wick.gamma_dot(3, h, V, tab, M) - truth).max()
    print("%-22s | %11.5f %11.3e | %11.3e %11.3e %8.1f%%"
          % (label, np.abs(truth).max(), np.abs(g[4]).max(), r_exact, r_null,
             100.0 * r_null / np.abs(truth).max()))
    print("%-22s | |separable_4| = %.5f, which the truncation puts where the exact gamma_4 = 0"
          % ("", np.abs(sep4).max()))

# An independent data point with a random Hamiltonian, so the result is not geometry-specific.
print()
rng = np.random.default_rng(23)
M = 6
h = hamiltonian.random_h(M, rng)
V = hamiltonian.random_V(M, rng)
psi = fullci.random_state(M, N, rng)
H = fullci.build_H(h, V, M)
g = {k: fullci.rdm(psi, M, k) for k in (1, 2, 3, 4)}
l2 = az.cumulant(2, g[2], {("G", 1): g[1]})
l3 = az.cumulant(3, g[3], {("G", 1): g[1], ("L", 2): l2})
truth = fullci.rdm_dot(psi, H, M, 3)
tab = {("D", 1): np.eye(M, dtype=complex), ("G", 1): g[1], ("G", 2): g[2], ("G", 3): g[3]}
tab[("G", 4)] = np.zeros((M,) * 8, dtype=complex)
r_exact = np.abs(wick.gamma_dot(3, h, V, tab, M) - truth).max()
sep4 = separable4(g[1], l2, l3, M)
tab[("G", 4)] = sep4
r_null = np.abs(wick.gamma_dot(3, h, V, tab, M) - truth).max()
print("%-22s | %11.5f %11.3e | %11.3e %11.3e %8.1f%%"
      % ("random h,V  M=6 N=3", np.abs(truth).max(), np.abs(g[4]).max(),
         r_exact, r_null, 100.0 * r_null / np.abs(truth).max()))
print("%-22s | |separable_4| = %.5f" % ("", np.abs(sep4).max()))
