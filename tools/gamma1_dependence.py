#!/usr/bin/env python3
"""Does the bare-Wick gamma_2 EOM contain gamma_1 at all, and where does gamma_1 enter?

John asks how gamma_2 could be propagated alone when the EOM invokes gamma_1. Count the
factors by (name, rank) in each form of the equation, and check the normalization of the
contraction actually used in the propagating tests.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from qdiss.symbolic import mkno
from qdiss.representability.dqg import contract_to_1rdm
from tools.eom2_coherent import ONE, TWO


def inventory(terms):
    d = {}
    for s, f in terms:
        for n, r, u, l in f:
            d[(n, r)] = d.get((n, r), 0) + 1
    return d


def fmt(d):
    return ", ".join("%s%d x%d" % (k[0], k[1], v) for k, v in sorted(d.items())) or "none"


print("=== factor inventory of the gamma_2 EOM, bare-moment form ===")
print("   one-body piece : %s" % fmt(inventory(ONE)))
print("   two-body piece : %s" % fmt(inventory(TWO)))
tot = inventory(ONE); [tot.update({k: tot.get(k, 0) + v}) for k, v in inventory(TWO).items()]
print("   combined       : %s" % fmt(tot))
print("   gamma_1 factors present: %s" % ("YES" if ("G", 1) in tot else "NO"))
print()
print("So the bare-moment EOM couples gamma_2 to gamma_2 and gamma_3 only. That is the BBGKY")
print("structure: rank k couples to rank k and rank k+1, never downward. No gamma_1 appears.")
print()

cum_two = mkno.to_cumulant_form(TWO)
cum_one = mkno.to_cumulant_form(ONE)
print("=== same equation after cumulant substitution ===")
print("   one-body : %d terms, %s" % (len(cum_one), fmt(inventory(cum_one))))
print("   two-body : %d terms, %s" % (len(cum_two), fmt(inventory(cum_two))))
print("   gamma_1 is INTRODUCED by the substitution, not present beforehand.")
print()

M, N = 6, 3
dim = 1 << M


def ad(p, v):
    out = np.zeros_like(v); bit = 1 << p
    for s in np.nonzero(v)[0]:
        if s & bit: continue
        out[s | bit] += ((-1.0) ** bin(s & (bit - 1)).count("1")) * v[s]
    return out


def aa(p, v):
    out = np.zeros_like(v); bit = 1 << p
    for s in np.nonzero(v)[0]:
        if not (s & bit): continue
        out[s ^ bit] += ((-1.0) ** bin(s & (bit - 1)).count("1")) * v[s]
    return out


r = np.random.default_rng(11)
popc = np.array([bin(s).count("1") for s in range(dim)])
psi = np.where(popc == N, r.normal(size=dim) + 1j * r.normal(size=dim), 0.0)
psi = psi / np.linalg.norm(psi)
Q1 = np.array([aa(p, psi) for p in range(M)])
g1 = Q1.conj() @ Q1.T
Q2 = np.array([aa(v, aa(p, psi)) for p in range(M) for v in range(M)])
g2 = (Q2.conj() @ Q2.T).reshape(M, M, M, M)

print("=== the contraction used in the propagating tests ===")
c = contract_to_1rdm(g2, N)
print("   |contract_to_1rdm(gamma_2, N) - exact gamma_1| = %.3e" % np.abs(c - g1).max())
print("   manual (1/(N-1)) sum_x gamma_2[px,qx]          = %.3e"
      % np.abs(np.einsum("pxqx->pq", g2) / (N - 1) - g1).max())
print("   so gamma_1 in the propagation IS the trace of gamma_2, exact at any instant.")
