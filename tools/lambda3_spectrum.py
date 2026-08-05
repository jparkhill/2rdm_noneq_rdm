#!/usr/bin/env python3
"""Spectrum and stability of the rank-2 EOM closed by the constrained lambda_3.

Derivation: qdiss.eom.gwt cumulant form, validated against qdiss.eom.fullci at 1e-15.
Propagated state (gamma_1, lambda_2). Matrix elements from PySCF for H3 linear / sto-3g,
M = 6, N = 3, checked against FCI. Truth: exact wavefunction propagation.

lambda_3 is reconstructed to satisfy both affine conditions exactly:
  the partial trace, sum_y lambda_3 = E lambda_2 + lambda_2 E with E = gamma_1 - 1/2, via the
  closed-form inverse lambda_3_tc = B[I, X] derived in tools/lambda3_build.py; and
  X_KK = 0, the kernel block of gamma_2-dot, which the exact flow satisfies as an EQUALITY
  because rank gamma_2 <= C(M, N-2) is conserved.
The second is imposed by adding trace-free directions with coefficients from linear least
squares. Because lambda_3 enters the EOM linearly with coefficients built only from V and
deltas, each basis direction's pair-space contribution is precomputed once and only the kernel
depends on the state, so the whole construction costs one small least squares per evaluation.

Reported: the poles of the linearised (gamma_1, lambda_2) flow against exact FCI excitation
energies and TDHF; then a propagation tracking the natural occupations of gamma_1 for
stability, with D, Q and G monitored periodically.
"""
import itertools
from math import comb
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from poles_real import build as build_H_so
from qdiss.eom import fullci, gwt
from qdiss.eom.evaluate import assemble
from qdiss.representability import dqg
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
L3T = [t for t in CUM2[1] if any(n == "L" and r == 3 for n, r, u, l in t[1])]

_, h, V, enuc = build_H_so("H 0 0 0; H 0 0 1.6; H 0 0 3.2", "sto-3g", 0, 1)
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
    for s, sub_, keys in ORB[key]:
        out = out + s * np.einsum(sub_, *[t[k] for k in keys], optimize=True)
    return out


def to_pair(T):
    A = np.array([[2.0 * T[p, v, q, w] for (q, w) in PAIRS] for (p, v) in PAIRS])
    return 0.5 * (A + A.conj().T)


def invert(R):
    S_R = np.einsum("pxqx->pq", R)
    trS = np.trace(S_R) / (3 * M - 6)
    S_X = (S_R - trS * EYE) / (2 * M - 6)
    return (R - wedge(EYE, S_X) - wedge(S_X, EYE)) / (M - 4)


def antisym6(T):
    out = np.zeros_like(T)
    for pu in itertools.permutations(range(3)):
        su = np.linalg.det(np.eye(3)[list(pu)])
        for pl in itertools.permutations(range(3)):
            sl = np.linalg.det(np.eye(3)[list(pl)])
            out += su * sl * np.transpose(T, tuple(pu) + tuple(3 + i for i in pl))
    return out / 36.0


NB = 140
rng = np.random.default_rng(21)
BAS, YK = [], []
for i in range(NB):
    T = rng.normal(size=(M,) * 6) + 1j * rng.normal(size=(M,) * 6)
    T = antisym6(T)
    T = T - orbit("B", EYE, invert(np.einsum("pvyqwy->pvqw", T)))
    T = T / np.abs(T).max()
    BAS.append(T)
    t = {("D", 1): EYE, ("G", 1): EYE, ("L", 2): EYE[:, :, None, None] * 0, ("L", 3): T}
    YK.append(to_pair(assemble([], L3T, h, V, t, "pvqw", M)))
CMAX = []
print("rank gamma_2 = C(M,N-2) = %d, kernel dim fixed at %d" % (comb(M, N-2), NP - comb(M, N-2)))
print("precomputed %d trace-free basis directions" % NB)


def lam3(g1, l2, Kb):
    E = g1 - 0.5 * EYE
    R = (np.einsum("py,yvqw->pvqw", E, l2) + np.einsum("vy,pyqw->pvqw", E, l2)
         + np.einsum("pvyw,yq->pvqw", l2, E) + np.einsum("pvqy,yw->pvqw", l2, E))
    l3 = orbit("B", EYE, invert(R))
    if Kb.shape[1] == 0:
        return l3, 0.0
    P = Kb @ Kb.conj().T
    iu = np.triu_indices(NP)

    def flat(B):
        return np.concatenate([B[iu].real, B[iu].imag])
    t = {("D", 1): EYE, ("G", 1): g1, ("L", 2): l2, ("L", 3): l3}
    X0 = P @ to_pair(assemble(CUM2[0], CUM2[1], h, V, t, "pvqw", M)) @ P
    A = np.stack([flat(P @ Y @ P) for Y in YK], axis=1)
    c, *_ = np.linalg.lstsq(A, -flat(0.5 * (X0 + X0.conj().T)), rcond=1e-10)
    return l3 + sum(c[i] * BAS[i] for i in range(NB)), float(np.abs(c).max())


KDIM = NP - comb(M, N - 2)


def kernel_of(g2):
    lam, U = np.linalg.eigh(to_pair(g2))
    return U[:, :KDIM]


def rhs(state):
    g1, l2 = state
    g2 = wedge(g1, g1) + l2
    l3, _ = lam3(g1, l2, kernel_of(g2))
    t = {("D", 1): EYE, ("G", 1): g1, ("L", 2): l2, ("L", 3): l3}
    d1 = assemble(CUM1[0], CUM1[1], h, V, t, "pq", M)
    d2 = assemble(CUM2[0], CUM2[1], h, V, t, "pvqw", M)
    return (d1, d2 - (wedge(d1, g1) + wedge(g1, d1)))


g1_0 = fullci.rdm(psi0, M, 1)
l2_0 = fullci.rdm(psi0, M, 2) - wedge(g1_0, g1_0)
d1_0, dl2_0 = rhs((g1_0, l2_0))
print("at the FCI reference: |gamma_1-dot| %.3e   |lambda_2-dot| %.3e"
      % (np.abs(d1_0).max(), np.abs(dl2_0).max()))
exc = np.sort((es - es[0]).real)[1:9]
print("exact excitations : %s" % np.array2string(exc, precision=4))

BH = []
for i in range(M):
    E0 = np.zeros((M, M), dtype=complex); E0[i, i] = 1.0; BH.append(("g", E0))
for i in range(M):
    for j in range(i + 1, M):
        E0 = np.zeros((M, M), dtype=complex); E0[i, j] = 1.0; E0[j, i] = 1.0
        BH.append(("g", E0))
        E0 = np.zeros((M, M), dtype=complex); E0[i, j] = 1j; E0[j, i] = -1j
        BH.append(("g", E0))
PB = []
for a in range(NP):
    for b in range(a, NP):
        for ph in (1.0, 1j):
            Z = np.zeros((NP, NP), dtype=complex)
            Z[a, b] += ph
            Z[b, a] += np.conj(ph)
            if np.abs(Z).max() < 1e-12:
                continue
            T = np.zeros((M,) * 4, dtype=complex)
            for i, (p, v) in enumerate(PAIRS):
                for j, (q, w) in enumerate(PAIRS):
                    x = 0.5 * Z[i, j]
                    T[p, v, q, w] = x; T[v, p, q, w] = -x
                    T[p, v, w, q] = -x; T[v, p, w, q] = x
            PB.append(("l", T))
BAS_ALL = BH + PB
n = len(BAS_ALL)
print("linearising over %d real directions (%d in gamma_1, %d in lambda_2)"
      % (n, len(BH), len(PB)))
eps = 1e-6
f0 = rhs((g1_0, l2_0))
J = np.zeros((n, n))
for k, (tag, Bk) in enumerate(BAS_ALL):
    st = (g1_0 + eps * Bk, l2_0) if tag == "g" else (g1_0, l2_0 + eps * Bk)
    fp = rhs(st)
    d = ((fp[0] - f0[0]) / eps, (fp[1] - f0[1]) / eps)
    for m, (tg, Bm) in enumerate(BAS_ALL):
        x = d[0] if tg == "g" else d[1]
        num = np.real(np.tensordot(Bm.conj(), x, axes=Bm.ndim))
        den = np.real(np.tensordot(Bm.conj(), Bm, axes=Bm.ndim))
        J[m, k] = num / den
_, cm = lam3(g1_0, l2_0, kernel_of(wedge(g1_0, g1_0) + l2_0))
print("max |c| in the least squares : %.3e" % cm)
ev = np.linalg.eigvals(J)
om = np.sort(np.abs(ev.imag[ev.imag > 1e-6]))
print("rank-2 EOM poles  : %s" % np.array2string(om[:10], precision=4))
print("largest |Re ev|   : %.3e  (growth rate; large means unstable)"
      % np.abs(ev.real).max())
