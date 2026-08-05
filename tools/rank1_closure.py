#!/usr/bin/env python3
"""Apply the MK_gamma30 construction one rank down. Is it TDHF?

Derivation: qdiss.eom.wick / qdiss.eom.gwt (identical term lists), validated against
qdiss.eom.fullci at 1e-15. All statements are about the rank-1 equation.

The MK_gamma30 closure cannot act on the 1-RDM EOM: that equation contains no lambda_3,
which is exactly why propagating (gamma_1, lambda_2) with lambda_3 solved from gamma_3 = 0
was exact at N = 2.

The analogue one rank down replaces gamma_3 = 0, exact for N = 2, with gamma_2 = 0, exact
for N = 1, and solves it for the top cumulant:

    lambda_2 := -(gamma_1 wedge gamma_1)      so the implied gamma_2 is exactly zero

Three rank-1 equations are compared:

    (a) exact         gamma_2 from the wavefunction
    (b) lambda_2 = 0  gamma_2 = gamma_1 wedge gamma_1        the TDHF candidate
    (c) gamma_2 = 0   the MK_gamma30 analogue

The generator is EXTRACTED rather than guessed. gamma_1[p,q] = <a+_p a_q> puts the bra
index on the row, so a unitary flow appears as gamma_1-dot = i[A, gamma_1] with A Hermitian.
In the eigenbasis of gamma_1 that reads d_ij = i A_ij (lam_j - lam_i), so the flow is a
commutator if and only if d vanishes wherever lam_i = lam_j, and A is then determined off
the degenerate blocks. Extract A, verify it reproduces d, and only then compare it to
h-transpose and to a mean-field form.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
from qdiss.eom import fullci, gwt, wick
from qdiss.eom.hamiltonian import random_V, random_h

M, N = 6, 3
EYE = np.eye(M, dtype=complex)


def wedge(a, b):
    return np.einsum("pq,vw->pvqw", a, b) - np.einsum("pw,vq->pvqw", a, b)


c1, c2 = gwt.terms_cumulant(1)
inv = {}
for s, f in list(c1) + list(c2):
    for n, r, u, l in f:
        inv[(n, r)] = inv.get((n, r), 0) + 1
print("rank-1 EOM, MKNO cumulant form: %d + %d terms, factors %s"
      % (len(c1), len(c2), ", ".join("%s%d x%d" % (k[0], k[1], v)
                                     for k, v in sorted(inv.items()))))
print("lambda_3 present: %s" % ("YES" if ("L", 3) in inv else
                                "NO -- the MK_gamma30 closure cannot act on the 1-RDM EOM"))
print()

rng = np.random.default_rng(21)
h, V = random_h(M, rng), random_V(M, rng)
psi = fullci.random_state(M, N, rng)
g1 = fullci.rdm(psi, M, 1)
g2 = fullci.rdm(psi, M, 2)
lam, U = np.linalg.eigh(0.5 * (g1 + g1.conj().T))


def extract(d):
    """If d = i[A, gamma_1] with A Hermitian, return (degeneracy residual, A, rebuild error)."""
    de = U.conj().T @ d @ U
    n = len(lam)
    mask = np.abs(lam[None, :] - lam[:, None]) < 1e-9
    resid = float(np.abs(de[mask]).max())
    A = np.zeros_like(de)
    for i in range(n):
        for j in range(n):
            if not mask[i, j]:
                A[i, j] = de[i, j] / (1j * (lam[j] - lam[i]))
    Af = U @ A @ U.conj().T
    rebuild = float(np.abs(1j * (Af @ g1 - g1 @ Af) - d).max())
    return resid, Af, rebuild


def d1_of(g2u):
    return wick.gamma_dot(1, h, V, {("D", 1): EYE, ("G", 1): g1, ("G", 2): g2u}, M)


W = np.einsum("prqs,sr->pq", V, g1)
cases = (("(a) exact", g2), ("(b) lambda_2 = 0", wedge(g1, g1)),
         ("(c) gamma_2 = 0", np.zeros_like(g2)))
print("=== is each rank-1 flow a commutator, and what is its generator? ===")
print("%-18s | %9s %10s %10s | %11s %11s %8s"
      % ("case", "|g1dot|", "degen res", "rebuild", "A vs h^T", "A vs f^T", "kappa"))
for name, g2u in cases:
    d = d1_of(g2u)
    resid, A, rebuild = extract(d)
    ehT = float(np.abs(A - h.conj()).max())
    best, bk = None, None
    for kappa in (0.25, 0.5, 1.0, 2.0):
        e = float(np.abs(A - (h + kappa * W).conj()).max())
        if best is None or e < best:
            best, bk = e, kappa
    print("%-18s | %9.5f %10.2e %10.2e | %11.3e %11.3e %8.2f"
          % (name, np.abs(d).max(), resid, rebuild, ehT, best, bk))
print()
print("degen res near zero means the flow moves no occupation, i.e. it IS a commutator;")
print("rebuild is the error of i[A, gamma_1] against the term list.")
print()

print("=== rank-2 comparison, generator extracted the same way ===")
g3 = fullci.rdm(psi, M, 3)
PAIRS = [(p, v) for p in range(M) for v in range(p + 1, M)]


def to_pair(T):
    return np.array([[2.0 * T[p, v, q, w] for (q, w) in PAIRS] for (p, v) in PAIRS])


rho = 0.5 * (to_pair(g2) + to_pair(g2).conj().T)
lam2, U2 = np.linalg.eigh(rho)
H2 = np.zeros((len(PAIRS), len(PAIRS)), dtype=complex)
for i, (p, v) in enumerate(PAIRS):
    for j, (q, w) in enumerate(PAIRS):
        H2[i, j] = (h[p, q] * (v == w) + h[v, w] * (p == q)
                    - h[p, w] * (v == q) - h[v, q] * (p == w) + V[p, v, q, w])
for nm, g3u in (("gamma_3 = 0", np.zeros_like(g3)), ("exact", g3)):
    d2 = wick.gamma_dot(2, h, V, {("D", 1): EYE, ("G", 1): g1,
                                  ("G", 2): g2, ("G", 3): g3u}, M)
    X = 0.5 * (to_pair(d2) + to_pair(d2).conj().T)
    de = U2.conj().T @ X @ U2
    mask = np.abs(lam2[None, :] - lam2[:, None]) < 1e-9
    print("   %-12s degen residual %.3e   |X - i[H_2^T, rho]| %.3e   (|X| %.4f)"
          % (nm, np.abs(de[mask]).max(),
             np.abs(X - 1j * (H2.conj() @ rho - rho @ H2.conj())).max(), np.abs(X).max()))
