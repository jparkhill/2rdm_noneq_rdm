#!/usr/bin/env python3
"""Is the tangency-projected closure Lipschitz in the state, or does it jump when the rank does?

THE WORRY. The correction Delta that enforces exact tangency is solved against Kmat, the matrix
whose columns are the eigenvectors of Mat_2(gamma_2)[I,J] = 2 gamma_2[i1,i2,j1,j2] (sorted pairs
I, J) having numerically zero eigenvalue. Kmat is defined by a THRESHOLD on that spectrum. Where
an eigenvalue crosses the threshold the number of columns changes discontinuously and Delta can
jump with it -- exactly the failure that disqualified the unregularised adiabatic closure `pt`,
whose pseudo-inverse changed rank as gamma_1 varied.

TWO SWEEPS, chosen so one cannot change the rank and the other must.

  pure path   psi(theta) = cos(theta) Psi_0 + sin(theta) Psi_1, a PURE state throughout. Then
              Mat_2(gamma_2) is the Gram matrix of the doubly annihilated state, which lives in
              the (N-2)-electron sector, so rank <= C(M,N-2) = 6 for every theta and the kernel
              dimension is pinned at 15 - 6 = 9. Expect smooth behaviour.

  mixed path  rho(p) = (1-p)|Psi_0><Psi_0| + p|Psi_1><Psi_1|. Now Mat_2 is the Gram matrix of the
              operators (a a) sqrt(rho), so rank <= C(M,N-2) * rank(rho) = 12, and the kernel
              collapses from 9 to 3 as p leaves zero. This is where a jump must appear if the
              construction has one.

MEASURED per point: the kernel dimension, the norm of the correction, and a finite-difference
Lipschitz estimate |Delta(x+eps) - Delta(x)| / |gamma_2(x+eps) - gamma_2(x)| between consecutive
points. A bounded estimate on the pure path and a divergent one on the mixed path would localise
the problem precisely at the rank change.

System (M,N) = (6,3) with a random Hamiltonian, matching tools/project_tangency.py. Starting
closure `separable`, i.e. gamma_3 from the MKNO cumulant expansion with lambda_3 dropped.
Helper routines are duplicated from project_tangency.py rather than imported, because that file
executes its analysis at module level.
"""
import itertools
import sys
import numpy as np
sys.path.insert(0, "."); sys.path.insert(0, "tools")
from qdiss.eom import fullci, wick
from qdiss.eom.hamiltonian import random_h, random_V
from qdiss.symbolic import mkno

M, N = 6, 3
C2 = list(itertools.combinations(range(M), 2))
C3 = list(itertools.combinations(range(M), 3))
EYE = np.eye(M, dtype=complex)


def perm_sign(p):
    s, p = 1, list(p)
    for i in range(len(p)):
        for j in range(len(p) - 1 - i):
            if p[j] > p[j + 1]:
                p[j], p[j + 1] = p[j + 1], p[j]; s = -s
    return s


def ten3(Z):
    T = np.zeros((M,) * 6, dtype=complex)
    for i, I in enumerate(C3):
        for j, J in enumerate(C3):
            z = Z[i, j] / 6.0
            for su in itertools.permutations(range(3)):
                pu = perm_sign(su); Iu = tuple(I[a] for a in su)
                for sl in itertools.permutations(range(3)):
                    T[Iu + tuple(J[b] for b in sl)] = pu * perm_sign(sl) * z
    return T


def mat2(T):
    return np.array([[2.0 * T[I + J] for J in C2] for I in C2])


rng = np.random.default_rng(3)
h, V = random_h(M, rng), random_V(M, rng)
dets = [d for d in range(1 << M) if bin(d).count("1") == N]
es, vec = np.linalg.eigh(fullci.build_H(h, V, M)[np.ix_(dets, dets)])
PSI = []
for k in (0, 1):
    p = np.zeros(1 << M, dtype=complex); p[dets] = vec[:, k]; PSI.append(p)

SEP = [(s, ",".join("".join(u) + "".join(l) for n, r, u, l in f) + "->pvxqwy",
        [(n, r) for n, r, u, l in f])
       for s, f in mkno.gamma_to_cumulants(("p", "v", "x"), ("q", "w", "y"))
       if not any(n == "L" and r == 3 for n, r, u, l in f)]
COLL = [(s, ["".join(u) + "".join(l) for n, r, u, l in f], [(n, r) for n, r, u, l in f])
        for s, f in wick.terms(2)[1]
        if ("G", 3) in [(n, r) for n, r, u, l in f]]
BAS = []
n3 = len(C3)
for a in range(n3):
    E = np.zeros((n3, n3), dtype=complex); E[a, a] = 1.0; BAS.append(E)
for a in range(n3):
    for b in range(a + 1, n3):
        E = np.zeros((n3, n3), dtype=complex); E[a, b] = 1.0; E[b, a] = 1.0; BAS.append(E)
        E = np.zeros((n3, n3), dtype=complex); E[a, b] = 1j; E[b, a] = -1j; BAS.append(E)
TEN_BAS = [ten3(Z) for Z in BAS]


def coll2(g3):
    t = {("D", 1): EYE, ("G", 3): g3}
    tot = np.zeros((M,) * 4, dtype=complex)
    for s, subs, keys in COLL:
        tot += s * 0.25 * np.einsum(",".join(["rtsu"] + subs) + "->pvqw", V,
                                    *[t[k] for k in keys], optimize=True)
    return 1j * tot


COLL_BAS = [mat2(coll2(T)) for T in TEN_BAS]
TR_BAS = [np.einsum("pvyqwy->pvqw", T) for T in TEN_BAS]


def state(weights, kets):
    g = {}
    for r in (1, 2, 3):
        acc = np.zeros((M,) * (2 * r), dtype=complex)
        for wgt, psi in zip(weights, kets):
            acc += wgt * fullci.rdm(psi, M, r)
        g[r] = acc
    return g


def solve(g, tol=1e-9):
    tab = {("D", 1): EYE, ("G", 1): g[1], ("G", 2): g[2]}
    tab[("L", 2)] = g[2] - (np.einsum("pq,vw->pvqw", g[1], g[1])
                            - np.einsum("pw,vq->pvqw", g[1], g[1]))
    g3s = sum(s * np.einsum(sub, *[tab[k] for k in keys], optimize=True) for s, sub, keys in SEP)
    X = mat2(g[2]); ev, U = np.linalg.eigh(0.5 * (X + X.conj().T))
    Km = U[:, ev < tol * max(ev.max(), 1.0)]
    d = Km.shape[1]
    if d == 0:
        return np.zeros((M,) * 6, dtype=complex), 0, 0.0
    iu = np.triu_indices(d)
    iu1 = np.triu_indices(d, 1)

    def cz(A):
        B = Km.conj().T @ A @ Km; B = 0.5 * (B + B.conj().T)
        return np.concatenate([np.real(B[iu]), np.imag(B[iu1])])

    A = np.vstack([np.array([cz(Cb) for Cb in COLL_BAS]).T,
                   np.array([np.concatenate([np.real(T).ravel(), np.imag(T).ravel()])
                             for T in TR_BAS]).T])
    c0 = cz(mat2(coll2(g3s)))
    b = np.concatenate([-c0, np.zeros(A.shape[0] - len(c0))])
    x = np.linalg.lstsq(A, b, rcond=None)[0]
    return sum(xi * T for xi, T in zip(x, TEN_BAS)), d, np.abs(cz(mat2(coll2(g3s)))).max()


def sweep(tag, params, make):
    print("%s" % tag)
    print("   %10s %8s %12s %12s %14s"
          % ("param", "ker dim", "|Delta|", "resid before", "Lipschitz est"))
    prev = None
    for t in params:
        g = make(t)
        D, d, r0 = solve(g)
        lip = ""
        if prev is not None:
            dg = np.abs(g[2] - prev[1][2]).max()
            lip = "%.3e" % (np.abs(D - prev[0]).max() / dg) if dg > 0 else "-"
        print("   %10.3e %8d %12.4e %12.3e %14s" % (t, d, np.abs(D).max(), r0, lip))
        prev = (D, g)
    print()


sweep("PURE path  psi = cos(th) Psi_0 + sin(th) Psi_1   (rank cannot change)",
      [0.0, 1e-4, 2e-4, 1e-3, 2e-3, 1e-2],
      lambda th: state([1.0], [(np.cos(th) * PSI[0] + np.sin(th) * PSI[1])
                               / np.linalg.norm(np.cos(th) * PSI[0] + np.sin(th) * PSI[1])]))
sweep("MIXED path  rho = (1-p) Psi_0 + p Psi_1          (rank must change)",
      [0.0, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1e-1],
      lambda p: state([1.0 - p, p], PSI))
