#!/usr/bin/env python3
"""Which physical properties of the closed RDM EOM survive which reconstructions.

The question is inverted from the usual one: for each known physical condition on the closing
equation -- hermiticity, antisymmetry, particle-number/trace conservation, contraction
consistency with the lower rank, energy conservation, positivity, and existence of the exact
ground state as a fixed point -- what does it actually demand of the cumulant reconstruction?
Some are automatic from the structure of the commutator and cost the reconstruction nothing;
others are genuine constraints. Measuring separates the two.

System: H3 linear / sto-3g, M=6 N=3, PySCF integrals checked against FCI. Derivation
qdiss.eom.gwt cumulant form, validated against qdiss.eom.fullci at 1e-15. State (gamma_1,
lambda_2); the reconstructed object is lambda_3. Evaluated at the exact FCI ground state, where
every listed quantity vanishes for the exact lambda_3, so any nonzero entry is closure error.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from math import comb
from poles_real import build as build_so
from qdiss.eom import fullci, gwt
from qdiss.eom.evaluate import assemble
from qdiss.eom.hamiltonian import TWO_BODY_PREFACTOR
from qdiss.symbolic import mkno

M, N = 6, 3
EYE = np.eye(M, dtype=complex)
PAIRS = [(p, v) for p in range(M) for v in range(p + 1, M)]
NP = len(PAIRS)
UP3, LO3 = ("p", "v", "x"), ("q", "w", "y")
ORB = {"A": [], "B": []}
for s, f in mkno.gamma_to_cumulants(UP3, LO3):
    if any(n == "L" and r == 3 for n, r, u, l in f):
        continue
    k = "B" if any(n == "L" and r == 2 for n, r, u, l in f) else "A"
    ORB[k].append((s, ",".join("".join(u) + "".join(l) for n, r, u, l in f) + "->pvxqwy",
                   [(n, r) for n, r, u, l in f]))
CUM1, CUM2 = gwt.terms_cumulant(1), gwt.terms_cumulant(2)

_, h, V, enuc = build_so("H 0 0 0; H 0 0 1.6; H 0 0 3.2", "sto-3g", 0, 1)
Hm = fullci.build_H(h, V, M)
popc = np.array([bin(s).count("1") for s in range(1 << M)])
sub = np.where(popc == N)[0]
es, vec = np.linalg.eigh(Hm[np.ix_(sub, sub)])
psi0 = np.zeros(1 << M, dtype=complex)
psi0[sub] = vec[:, 0]


def wedge(a, b):
    return np.einsum("pq,vw->pvqw", a, b) - np.einsum("pw,vq->pvqw", a, b)


def orbit(key, a, X):
    t = {("G", 1): a, ("L", 2): X}
    out = np.zeros((M,) * 6, dtype=complex)
    for s, sb, keys in ORB[key]:
        out = out + s * np.einsum(sb, *[t[k] for k in keys], optimize=True)
    return out


def to_pair(T):
    A = np.array([[2.0 * T[p, v, q, w] for (q, w) in PAIRS] for (p, v) in PAIRS])
    return 0.5 * (A + A.conj().T)


def invert(R):
    S = np.einsum("pxqx->pq", R)
    tr = np.trace(S) / (3 * M - 6)
    SX = (S - tr * EYE) / (2 * M - 6)
    return (R - wedge(EYE, SX) - wedge(SX, EYE)) / (M - 4)


g1 = fullci.rdm(psi0, M, 1)
g2 = fullci.rdm(psi0, M, 2)
g3 = fullci.rdm(psi0, M, 3)
l2 = g2 - wedge(g1, g1)
l3x = g3 - orbit("A", g1, l2) - orbit("B", g1, l2)
E = g1 - 0.5 * EYE
R3 = (np.einsum("py,yvqw->pvqw", E, l2) + np.einsum("vy,pyqw->pvqw", E, l2)
      + np.einsum("pvyw,yq->pvqw", l2, E) + np.einsum("pvqy,yw->pvqw", l2, E))
l3_tc = orbit("B", EYE, invert(R3))
KDIM = NP - comb(M, N - 2)
lam_, U_ = np.linalg.eigh(to_pair(g2))
P = U_[:, :KDIM] @ U_[:, :KDIM].conj().T


def flows(l3):
    t = {("D", 1): EYE, ("G", 1): g1, ("L", 2): l2, ("L", 3): l3}
    d1 = assemble(CUM1[0], CUM1[1], h, V, t, "pq", M)
    d2 = assemble(CUM2[0], CUM2[1], h, V, t, "pvqw", M)
    return d1, d2


print("H3 linear / sto-3g, M=%d N=%d, at the exact FCI ground state" % (M, N))
print("every column must be zero for the exact lambda_3; nonzero is closure error")
print()
hdr = ("%-22s | %9s %9s | %9s %9s | %9s %9s %9s"
       % ("lambda_3", "herm g2d", "antisym", "tr g1dot", "tr g2dot", "dE/dt",
          "contract", "min X_KK"))
print(hdr)
print("-" * len(hdr))
for nm, l3 in (("exact", l3x),
               ("zero", np.zeros((M,) * 6, dtype=complex)),
               ("-(A+B)  MK_gamma30", -(orbit("A", g1, l2) + orbit("B", g1, l2))),
               ("closed-form trace-exact", l3_tc)):
    d1, d2 = flows(l3)
    herm = np.abs(d2 - np.transpose(d2.conj(), (2, 3, 0, 1))).max()
    anti = np.abs(d2 + np.transpose(d2, (1, 0, 2, 3))).max()
    tr1 = abs(np.trace(d1))
    tr2 = abs(np.einsum("pvpv->", d2))
    dE = abs(float(np.real(np.einsum("pq,pq->", h, d1)
                           + TWO_BODY_PREFACTOR * np.einsum("pvqw,pvqw->", V, d2))))
    contr = np.abs(np.einsum("pxqx->pq", d2) / (N - 1) - d1).max()
    X = to_pair(d2)
    xk = np.linalg.eigvalsh(0.5 * (P @ X @ P + (P @ X @ P).conj().T)).min()
    print("%-22s | %9.2e %9.2e | %9.2e %9.2e | %9.2e %9.2e %9.2e"
          % (nm, herm, anti, tr1, tr2, dE, contr, xk))
print()
print("herm/antisym: structural properties of gamma_2-dot")
print("tr g1dot, tr g2dot: particle-number conservation at each rank")
print("dE/dt: energy conservation of the exact functional on the propagated RDMs")
print("contract: |contract(gamma_2-dot)/(N-1) - gamma_1-dot|, BBGKY consistency")
print("min X_KK: Nagumo / CP condition on the kernel of gamma_2")
