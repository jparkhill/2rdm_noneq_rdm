#!/usr/bin/env python3
"""Direct test: does each rank-1 closure equal i[A, gamma_1] for a named generator A?

Comparing an extracted generator to a candidate is wrong, since a commutator fixes its
generator only off the degenerate blocks. Test the flow itself instead. Convention, fixed by
the rank-2 result |X - i[H_2^T, rho]| = 1.2e-15: gamma_1[p,q] = <a+_p a_q> puts the bra index
on the row, so a unitary flow reads gamma_1-dot = i[A^T, gamma_1].

Ground truth is the validated term list from qdiss.eom.wick; the candidates are identified
against it, never the reverse.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
from qdiss.eom import fullci, wick
from qdiss.eom.hamiltonian import random_V, random_h

M, N = 6, 3
EYE = np.eye(M, dtype=complex)


def wedge(a, b):
    return np.einsum("pq,vw->pvqw", a, b) - np.einsum("pw,vq->pvqw", a, b)


rng = np.random.default_rng(21)
h, V = random_h(M, rng), random_V(M, rng)
psi = fullci.random_state(M, N, rng)
g1 = fullci.rdm(psi, M, 1)
g2 = fullci.rdm(psi, M, 2)


def d1_of(g2u):
    return wick.gamma_dot(1, h, V, {("D", 1): EYE, ("G", 1): g1, ("G", 2): g2u}, M)


def resid(d, A):
    return float(np.abs(d - 1j * (A.conj() @ g1 - g1 @ A.conj())).max())


MF = {
    "V[p,r,q,s] g1[s,r]": np.einsum("prqs,sr->pq", V, g1),
    "V[p,r,q,s] g1[r,s]": np.einsum("prqs,rs->pq", V, g1),
    "V[r,p,s,q] g1[s,r]": np.einsum("rpsq,sr->pq", V, g1),
}
for name, g2u in (("(c) gamma_2 = 0", np.zeros_like(g2)),
                  ("(b) lambda_2 = 0", wedge(g1, g1)),
                  ("(a) exact", g2)):
    d = d1_of(g2u)
    print("%-18s |gamma_1-dot| = %.5f" % (name, np.abs(d).max()))
    print("      A = h                        residual %.3e" % resid(d, h))
    best = None
    for mname, W in MF.items():
        for kappa in (0.25, 0.5, 1.0):
            r = resid(d, h + kappa * W)
            if best is None or r < best[0]:
                best = (r, mname, kappa)
    print("      A = h + kappa*W  best        residual %.3e   (%s, kappa=%.2f)"
          % (best[0], best[1], best[2]))
print()
print("TDHF is i[f^T, gamma_1] with f the Fock operator, so case (b) matching a mean-field")
print("generator at machine precision IS TDHF; case (c) matching bare h is free evolution.")
