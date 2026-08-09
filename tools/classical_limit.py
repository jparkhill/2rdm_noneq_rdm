#!/usr/bin/env python3
"""The strong-bath classical limit: does dephasing rescue a closure, and by which mechanism?

THE ARGUMENT BEING TESTED. Pure dephasing in the occupation-number basis damps the off-diagonal
elements of rho between determinants. In the limit rho is diagonal -- a classical probability
distribution over occupation configurations -- and then gamma_2[p,v,q,w] = <a+_p a+_v a_w a_q> is
nonzero only when the index sets {p,v} and {q,w} coincide, so Mat_2(gamma_2)[I,J] = 2 gamma_2[I+J]
becomes DIAGONAL with entries 2<n_p n_v>. The positive semidefinite cone degenerates to the
non-negative orthant, where inwardness is automatic for a classical master equation: at a face
with <n_p n_v> = 0 every configuration containing that pair has zero probability, so the loss term
vanishes identically and only non-negative gain survives. That is the Metzler condition, and it
requires nothing of the closure. Hence the weakest-bath strength gamma* is finite: full dephasing
certainly suffices.

TWO CANDIDATE MECHANISMS, which this file separates.

  (A) the closure's defect becomes harmless -- the compression C(X) = P_K F[X] P_K stops having
      negative eigenvalues as coherences are damped, so tangency is satisfied rather than evaded;
  (B) the kernel is LIFTED -- Mat_2(gamma_2) becomes strictly positive definite, the boundary is
      simply never reached, and the tangency condition becomes vacuous.

(B) would parallel the earlier finding that mixing raises the positivity clearance by five orders
of magnitude, making the constraint inactive rather than satisfied.

THE FAMILY. rho(s) has the Fock-basis matrix elements of the pure ground state with every
off-diagonal multiplied by (1 - s), so s = 0 is the pure state and s = 1 is the classical
distribution over determinants with the same populations. gamma_1, gamma_2 and gamma_3 are linear
in rho and are accumulated directly from that matrix. System (M,N) = (6,3) with a random
Hamiltonian, matching tools/project_tangency.py, chosen because rank Mat_2 <= C(M,N-2) = 6 of 15
forces an exact 9-dimensional kernel at s = 0.

REPORTED per s: the kernel dimension and clearance of Mat_2(gamma_2), and the extreme eigenvalues
of the compression of the closure's gamma_3 terms, for `separable` (lambda_3 = 0) and
`muller2405`. The exact gamma_3 is shown as the control and must stay at zero.
"""
import itertools
import sys
import numpy as np
sys.path.insert(0, "."); sys.path.insert(0, "tools")
from qdiss import muller2405
from qdiss.eom import fullci, wick
from qdiss.eom.hamiltonian import random_h, random_V
from qdiss.symbolic import mkno

M, N = 6, 3
C2 = list(itertools.combinations(range(M), 2))
EYE = np.eye(M, dtype=complex)
rng = np.random.default_rng(3)
h, V = random_h(M, rng), random_V(M, rng)
dets = [d for d in range(1 << M) if bin(d).count("1") == N]
es, vec = np.linalg.eigh(fullci.build_H(h, V, M)[np.ix_(dets, dets)])
c0 = vec[:, 0]

SEP = [(s, ",".join("".join(u) + "".join(l) for n, r, u, l in f) + "->pvxqwy",
        [(n, r) for n, r, u, l in f])
       for s, f in mkno.gamma_to_cumulants(("p", "v", "x"), ("q", "w", "y"))
       if not any(n == "L" and r == 3 for n, r, u, l in f)]
COLL = [(s, ["".join(u) + "".join(l) for n, r, u, l in f], [(n, r) for n, r, u, l in f])
        for s, f in wick.terms(2)[1]
        if ("G", 3) in [(n, r) for n, r, u, l in f]]


def mat2(T):
    return np.array([[2.0 * T[I + J] for J in C2] for I in C2])


def coll2(g3):
    t = {("D", 1): EYE, ("G", 3): g3}
    tot = np.zeros((M,) * 4, dtype=complex)
    for s, subs, keys in COLL:
        tot += s * 0.25 * np.einsum(",".join(["rtsu"] + subs) + "->pvqw", V,
                                    *[t[k] for k in keys], optimize=True)
    return 1j * tot


def rdms(s):
    """gamma_1..3 of rho(s): the pure-state density matrix with off-diagonals damped by (1-s)."""
    R = np.outer(c0, c0.conj()) * (1.0 - s)
    np.fill_diagonal(R, np.abs(c0) ** 2)
    g = {}
    w, U = np.linalg.eigh(R)
    keep = np.where(w > 1e-14)[0]
    for r in (1, 2, 3):
        acc = np.zeros((M,) * (2 * r), dtype=complex)
        for a in keep:
            psi = np.zeros(1 << M, dtype=complex); psi[dets] = U[:, a]
            acc += w[a] * fullci.rdm(psi, M, r)
        g[r] = acc
    return g


print("strong-bath classical limit, (M,N) = (%d,%d);  s = 0 pure, s = 1 fully dephased" % (M, N))
print("   %6s %8s %12s | %22s | %22s | %12s"
      % ("s", "ker dim", "clearance", "separable  min / max", "muller2405 min / max", "exact max"))
for s in (0.0, 1e-4, 1e-2, 0.1, 0.5, 0.9, 0.99, 1.0):
    g = rdms(s)
    X = mat2(g[2]); ev, U = np.linalg.eigh(0.5 * (X + X.conj().T))
    Km = U[:, ev < 1e-9 * max(ev.max(), 1.0)]
    d = Km.shape[1]
    tab = {("D", 1): EYE, ("G", 1): g[1], ("G", 2): g[2]}
    tab[("L", 2)] = g[2] - (np.einsum("pq,vw->pvqw", g[1], g[1])
                            - np.einsum("pw,vq->pvqw", g[1], g[1]))
    g3s = sum(sg * np.einsum(sub, *[tab[k] for k in keys], optimize=True)
              for sg, sub, keys in SEP)
    g3m = muller2405.gamma3(g[1], g[2], N)

    def comp(g3):
        if d == 0:
            return None
        B = Km.conj().T @ mat2(coll2(g3)) @ Km
        return np.linalg.eigvalsh(0.5 * (B + B.conj().T))

    cs, cm, ce = comp(g3s), comp(g3m), comp(g[3])
    fmt = (lambda e: "  vacuous (no kernel) " if e is None
           else "%+10.3e %+10.3e" % (e.min(), e.max()))
    print("   %6.4f %8d %12.3e | %22s | %22s | %12s"
          % (s, d, ev.min(), fmt(cs), fmt(cm),
             "n/a" if ce is None else "%.1e" % np.abs(ce).max()))
