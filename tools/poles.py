#!/usr/bin/env python3
"""Poles of the closed 1-RDM EOMs, against exact excitation energies.

Derivation: qdiss.eom.gwt cumulant rank-1 EOM, validated against qdiss.eom.fullci at 1e-15.
Propagated state gamma_1 alone with lambda_2 := R[gamma_1].

Procedure. Find a Hartree-Fock reference, an idempotent gamma_1 commuting with the Fock
generator, which is stationary under TDHF. Linearize gamma_1-dot about it by finite
differences over a Hermitian basis, and read the poles off the Jacobian eigenvalues, which
come in +-i omega pairs. Compare against the exact excitation energies of H in the N-electron
sector.

A caution checked first. The blocked reconstruction is lambda_2 = A^f + f^A with
A = gamma_1 - gamma_1^2, and at a Hartree-Fock reference A vanishes identically. Worse, its
coefficient f_i = -a_i/[(1+Q)(T-2a_i)] is 0/0 there, and approaching along a_i = eps b_i gives
f_i -> -b_i/[(1+Q)(sum b - 2 b_i)], which depends on the DIRECTION b. So lambda_2 is O(eps) but
not differentiable at the HF point: linear response about HF is ill-defined for it. Measured
below by approaching HF along several directions.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
from qdiss.eom import fullci, gwt
from qdiss.eom.evaluate import assemble
from qdiss.eom.hamiltonian import random_V, random_h

CUM1 = gwt.terms_cumulant(1)
M, N = 6, 3
EYE = np.eye(M, dtype=complex)


def wedge(a, b):
    return np.einsum("pq,vw->pvqw", a, b) - np.einsum("pw,vq->pvqw", a, b)


def sym_wedge(a, b):
    return wedge(a, b) + wedge(b, a)


def blocked(g1):
    A = 0.5 * ((g1 - g1 @ g1) + (g1 - g1 @ g1).conj().T)
    a, U = np.linalg.eigh(A)
    T = a.sum()
    den = T - 2 * a
    if np.abs(den).min() < 1e-12:
        return np.zeros((M,) * 4, dtype=complex)
    Q = (a / den).sum()
    f = -a / ((1.0 + Q) * den)
    return sym_wedge(A, U @ np.diag(f.astype(complex)) @ U.conj().T)


rng = np.random.default_rng(3)
h, V = random_h(M, rng), random_V(M, rng)
H = fullci.build_H(h, V, M)


def fock(g1):
    return h + np.einsum("prqs,rs->pq", V, g1)


g1 = np.zeros((M, M), dtype=complex)
ev, vec = np.linalg.eigh(fock(g1).conj())
g1 = vec[:, :N] @ vec[:, :N].conj().T
for _ in range(300):
    ev, vec = np.linalg.eigh(0.5 * (fock(g1).conj() + fock(g1).conj().conj().T))
    new = vec[:, :N] @ vec[:, :N].conj().T
    if np.abs(new - g1).max() < 1e-13:
        g1 = new
        break
    g1 = 0.5 * new + 0.5 * g1
F = fock(g1).conj()
print("HF reference: |g1^2 - g1| %.2e   |[F, g1]| %.2e   tr g1 %.4f"
      % (np.abs(g1 @ g1 - g1).max(), np.abs(F @ g1 - g1 @ F).max(), np.trace(g1).real))


def d1(x, kind):
    l2 = np.zeros((M,) * 4, dtype=complex) if kind == "tdhf" else blocked(x)
    t = {("D", 1): EYE, ("G", 1): x, ("L", 2): l2}
    return assemble(CUM1[0], CUM1[1], h, V, t, "pq", M)


print("   |gamma_1-dot| at HF, TDHF    %.3e" % np.abs(d1(g1, "tdhf")).max())
print("   |lambda_2 blocked| at HF     %.3e" % np.abs(blocked(g1)).max())
print()
print("approaching HF along random directions, is blocked's lambda_2/eps direction-dependent?")
r2 = np.random.default_rng(11)
for trial in range(3):
    Wd = r2.normal(size=(M, M)) + 1j * r2.normal(size=(M, M))
    Wd = 0.5 * (Wd + Wd.conj().T)
    Wd = Wd - np.trace(Wd) / M * np.eye(M)
    vals = []
    for eps in (1e-3, 1e-4):
        gx = g1 + eps * Wd
        vals.append(np.abs(blocked(gx)).max() / eps)
    print("   direction %d : |lambda_2|/eps = %.5f, %.5f" % (trial, vals[0], vals[1]))
print()

basis = []
for i in range(M):
    E0 = np.zeros((M, M), dtype=complex); E0[i, i] = 1.0
    basis.append(E0)
for i in range(M):
    for j in range(i + 1, M):
        E0 = np.zeros((M, M), dtype=complex); E0[i, j] = 1.0; E0[j, i] = 1.0
        basis.append(E0)
        E1 = np.zeros((M, M), dtype=complex); E1[i, j] = 1j; E1[j, i] = -1j
        basis.append(E1)


def jacobian(kind, eps=1e-6):
    n = len(basis)
    J = np.zeros((n, n))
    f0 = d1(g1, kind)
    for k, Bk in enumerate(basis):
        dp = d1(g1 + eps * Bk, kind) - f0
        dp = dp / eps
        for m, Bm in enumerate(basis):
            J[m, k] = float(np.real(np.trace(Bm.conj().T @ dp))) / float(
                np.real(np.trace(Bm.conj().T @ Bm)))
    return J


exact = np.linalg.eigvalsh(H)
pop = np.array([bin(s).count("1") for s in range(1 << M)])
sub = np.where(pop == N)[0]
Hs = H[np.ix_(sub, sub)]
es = np.linalg.eigvalsh(Hs)
exc = np.sort(es - es[0])[1:9]
print("exact excitation energies (N-electron sector): %s"
      % np.array2string(exc, precision=4))
for kind in ("tdhf", "blocked"):
    J = jacobian(kind)
    ev = np.linalg.eigvals(J)
    om = np.sort(np.abs(ev.imag[ev.imag > 1e-8]))
    print("%-8s poles: %s" % (kind, np.array2string(om[:8], precision=4)))
print()
print("blocked's lambda_2 vanishes at the HF reference, so its poles must coincide with TDHF")
print("there. A reconstruction that shifts poles has to be nonzero at the reference.")
