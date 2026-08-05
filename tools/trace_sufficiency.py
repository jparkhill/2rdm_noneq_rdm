#!/usr/bin/env python3
"""Is preserving the partial-trace identity SUFFICIENT for positivity?

John's objection is correct: the exact lambda_3 satisfies both the trace identity and
positivity, so asking whether the intersection is non-empty is vacuous. The useful question
is whether Rdot = 0 is sufficient, i.e. whether every lambda_3 preserving the identity also
keeps the flow inside the representable set.

That is decidable directly. Start from a representable state and the EXACT lambda_3, so both
conditions hold. Then perturb lambda_3 along directions in the kernel of the map

    L(lambda_3) = sum_x [ lambda_3-carrying terms of gamma_2-dot ][px,qx]

which is the only way lambda_3 enters Rdot, since gamma_1-dot is lambda_3-free. A perturbation
in ker L leaves Rdot exactly zero. If it can still drive the flow out of the PSD cone, then
Rdot = 0 is necessary but not sufficient and the search is real; if positivity is robust
against all of ker L, the identity is doing the whole job.

Positivity is measured by the Nagumo condition, which for the PSD cone is exact: with K the
kernel of gamma_2, the flow stays PSD iff the kernel block X_KK of gamma_2-dot is PSD. At
M=6 N=3 that kernel is 9-dimensional since rank gamma_2 <= C(M, N-2) = 6.

Derivation: qdiss.eom.gwt cumulant form, validated against qdiss.eom.fullci at 1e-15.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
from qdiss.eom import fullci, gwt
from qdiss.eom.evaluate import assemble
from qdiss.eom.hamiltonian import random_V, random_h
from qdiss.symbolic import mkno

M, N = 6, 3
EYE = np.eye(M, dtype=complex)
PAIRS = [(p, v) for p in range(M) for v in range(p + 1, M)]
UP3, LO3 = ("p", "v", "x"), ("q", "w", "y")
ORB = {"A": [], "B": []}
for s, f in mkno.gamma_to_cumulants(UP3, LO3):
    if any(n == "L" and r == 3 for n, r, u, l in f):
        continue
    k = "B" if any(n == "L" and r == 2 for n, r, u, l in f) else "A"
    ORB[k].append((s, ",".join("".join(u) + "".join(l) for n, r, u, l in f) + "->pvxqwy",
                   [(n, r) for n, r, u, l in f]))
CUM1, CUM2 = gwt.terms_cumulant(1), gwt.terms_cumulant(2)
L3T = [t for t in CUM2[1] if any(n == "L" and r == 3 for n, r, u, l in t[1])]


def wedge(a, b):
    return np.einsum("pq,vw->pvqw", a, b) - np.einsum("pw,vq->pvqw", a, b)


def orbit(key, g1, l2):
    t = {("G", 1): g1, ("L", 2): l2}
    out = np.zeros((M,) * 6, dtype=complex)
    for s, sub, keys in ORB[key]:
        out = out + s * np.einsum(sub, *[t[k] for k in keys], optimize=True)
    return out


def antisymmetrize(T):
    import itertools
    out = np.zeros_like(T)
    for pu in itertools.permutations(range(3)):
        su = np.linalg.det(np.eye(3)[list(pu)])
        for pl in itertools.permutations(range(3)):
            sl = np.linalg.det(np.eye(3)[list(pl)])
            out += su * sl * np.transpose(T, tuple(pu) + tuple(3 + i for i in pl))
    return out / 36.0


rng = np.random.default_rng(7)
h, V = random_h(M, rng), random_V(M, rng)
psi = fullci.random_state(M, N, rng)
g1 = fullci.rdm(psi, M, 1)
g2 = fullci.rdm(psi, M, 2)
g3 = fullci.rdm(psi, M, 3)
l2 = g2 - wedge(g1, g1)
l3x = g3 - orbit("A", g1, l2) - orbit("B", g1, l2)


def g2dot(l3):
    t = {("D", 1): EYE, ("G", 1): g1, ("L", 2): l2, ("L", 3): l3}
    return assemble(CUM2[0], CUM2[1], h, V, t, "pvqw", M)


def l3_only(l3):
    t = {("D", 1): EYE, ("G", 1): g1, ("L", 2): l2, ("L", 3): l3}
    return assemble([], L3T, h, V, t, "pvqw", M)


def Lmap(l3):
    return np.einsum("pxqx->pq", l3_only(l3))


def to_pair(T):
    A = np.array([[2.0 * T[p, v, q, w] for (q, w) in PAIRS] for (p, v) in PAIRS])
    return 0.5 * (A + A.conj().T)


rho = to_pair(g2)
lam, U = np.linalg.eigh(rho)
Kb = U[:, lam < 1e-9]
print("pair dim %d   rank gamma_2 %d   kernel dim %d"
      % (len(PAIRS), int((lam > 1e-9).sum()), Kb.shape[1]))


def kernel_block(l3):
    X = to_pair(g2dot(l3))
    XKK = Kb.conj().T @ X @ Kb
    return np.linalg.eigvalsh(0.5 * (XKK + XKK.conj().T))


print()
print("exact lambda_3 :  |Rdot contribution| %.3e   min eig X_KK %.3e"
      % (np.abs(Lmap(l3x) - Lmap(l3x)).max(), kernel_block(l3x).min()))
print()

NB = 80
basis, cols = [], []
r2 = np.random.default_rng(99)
for k in range(NB):
    T = r2.normal(size=(M,) * 6) + 1j * r2.normal(size=(M,) * 6)
    T = antisymmetrize(T)
    T = T / np.abs(T).max()
    basis.append(T)
    cols.append(Lmap(T).ravel())
Amat = np.stack(cols, axis=1)
Amat = np.concatenate([Amat.real, Amat.imag], axis=0)
u, sv, vt = np.linalg.svd(Amat, full_matrices=True)
nullvecs = vt[np.sum(sv > 1e-8 * sv.max()):]
print("basis %d, rank of L on it %d, kernel dim %d"
      % (NB, int((sv > 1e-8 * sv.max()).sum()), nullvecs.shape[0]))
print()
print("%10s | %14s %14s | %14s %14s"
      % ("eps", "kerL |Rdot|", "kerL min X_KK", "generic |Rdot|", "gen min X_KK"))
Kdir = sum(nullvecs[0][k] * basis[k] for k in range(NB))
Kdir = Kdir / np.abs(Kdir).max()
Gdir = basis[0]
for eps in (0.0, 0.05, 0.2, 1.0):
    lk = l3x + eps * Kdir
    lg = l3x + eps * Gdir
    print("%10.2f | %14.3e %14.6f | %14.3e %14.6f"
          % (eps, np.abs(Lmap(lk) - Lmap(l3x)).max(), kernel_block(lk).min(),
             np.abs(Lmap(lg) - Lmap(l3x)).max(), kernel_block(lg).min()))
print()
print("kerL perturbations keep Rdot unchanged by construction. If their min X_KK goes")
print("negative, the trace identity is necessary but NOT sufficient for positivity.")
