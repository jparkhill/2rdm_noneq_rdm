#!/usr/bin/env python3
"""Project a gamma_3 closure onto the exact tangency condition, and measure what it costs.

THE CONDITION. Write Mat_2(T)[I,J] = 2 T[i1,i2,j1,j2] over sorted index pairs I=(i1<i2),
J=(j1<j2), and let Kmat be the C(M,2) x d matrix whose columns are the eigenvectors of
Mat_2(gamma_2) with numerically zero eigenvalue (d = dim of that kernel). The rank-2 collision
term, i.e. the 8 two-body terms of the exact EOM that carry gamma_3, is

  Coll_2[p,v,q,w] = i (1/4) sum_{r,t,s,u} V[r,t,s,u] (
        delta[p,s] gamma_3[r,t,v,q,u,w] - delta[p,u] gamma_3[r,t,v,q,s,w]
      - delta[r,q] gamma_3[p,t,v,s,u,w] + delta[r,w] gamma_3[p,t,v,q,s,u]
      + delta[t,q] gamma_3[p,r,v,s,u,w] - delta[t,w] gamma_3[p,r,v,q,s,u]
      - delta[v,s] gamma_3[p,r,t,q,u,w] + delta[v,u] gamma_3[p,r,t,q,s,w] )

and the exact flow satisfies, to machine precision,

  Kmat^dag Mat_2( Coll_2[gamma_3_exact] ) Kmat = 0                (d^2 real equations)

because the rank of Mat_2(gamma_2) is conserved under unitary evolution, so its kernel is
preserved. The drift terms compress to zero identically (Lyapunov form), so this is the whole
positivity requirement at this rank.

THE PROJECTION. Coll_2 is LINEAR in gamma_3, so writing gamma_3 = gamma_3_closure + Delta the
condition is a linear system for Delta:

  Kmat^dag Mat_2( Coll_2[Delta] ) Kmat = - Kmat^dag Mat_2( Coll_2[gamma_3_closure] ) Kmat

Delta is parameterised by its own matrix Mat_3(Delta), a C(M,3) x C(M,3) Hermitian matrix, which
automatically enforces antisymmetry and hermiticity of the tensor. A second block of linear
constraints requires Delta to have vanishing partial trace, sum_y Delta[p,v,y,q,w,y] = 0, so that
the rank-3 contraction identity satisfied by the starting closure is not broken by the fix. The
two blocks are solved together by minimum-norm least squares: either the intersection is nonempty
or it is not, and the residual says which.

Starting closure: `separable`, the MKNO cumulant expansion of gamma_3 with the 3-cumulant dropped
(lambda_3 = 0). System (M,N) = (6,3), chosen because rank Mat_2(gamma_2) <= C(M,N-2) = 6 of 15
forces an exact 9-dimensional kernel. Term list from qdiss.eom.wick, validated against
qdiss.eom.fullci at 1e-15.
"""
import itertools
import sys
import numpy as np
sys.path.insert(0, "."); sys.path.insert(0, "tools")
from qdiss.eom import fullci, wick
from qdiss.eom.hamiltonian import random_h, random_V
from qdiss.symbolic import mkno

M, N, K = 6, 3, 2
C2 = list(itertools.combinations(range(M), 2))
C3 = list(itertools.combinations(range(M), 3))


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
Hm = fullci.build_H(h, V, M)
es, vec = np.linalg.eigh(Hm[np.ix_(dets, dets)])
psi = np.zeros(1 << M, dtype=complex); psi[dets] = vec[:, 0]
g = {r: fullci.rdm(psi, M, r) for r in (1, 2, 3)}
EYE = np.eye(M, dtype=complex)

tab = {("D", 1): EYE, ("G", 1): g[1], ("G", 2): g[2]}
tab[("L", 2)] = g[2] - (np.einsum("pq,vw->pvqw", g[1], g[1])
                        - np.einsum("pw,vq->pvqw", g[1], g[1]))
SEP = [(s, ",".join("".join(u) + "".join(l) for n, r, u, l in f) + "->pvxqwy",
        [(n, r) for n, r, u, l in f])
       for s, f in mkno.gamma_to_cumulants(("p", "v", "x"), ("q", "w", "y"))
       if not any(n == "L" and r == 3 for n, r, u, l in f)]
g3_sep = sum(s * np.einsum(sub, *[tab[k] for k in keys], optimize=True) for s, sub, keys in SEP)

COLL = []
for s, f in wick.terms(2)[1]:
    keys = [(n, r) for n, r, u, l in f]
    if ("G", 3) in keys:
        COLL.append((s, ["".join(u) + "".join(l) for n, r, u, l in f], keys))


def coll2(g3):
    t = {("D", 1): EYE, ("G", 3): g3}
    tot = np.zeros((M,) * 4, dtype=complex)
    for s, subs, keys in COLL:
        tot += s * 0.25 * np.einsum(",".join(["rtsu"] + subs) + "->pvqw", V,
                                    *[t[k] for k in keys], optimize=True)
    return 1j * tot


X = mat2(g[2]); ev, U = np.linalg.eigh(0.5 * (X + X.conj().T))
Kmat = U[:, ev < 1e-9]
d = Kmat.shape[1]
print("(M,N) = (%d,%d);  Mat_2 is %dx%d, numerical rank %d, kernel dimension d = %d"
      % (M, N, len(C2), len(C2), int((ev > 1e-9).sum()), d))


def compress(T):
    B = Kmat.conj().T @ mat2(T) @ Kmat
    return 0.5 * (B + B.conj().T)


print("tangency residual |Kmat^dag Mat_2(Coll_2) Kmat|, exact vs starting closure:")
print("   exact gamma_3        %.3e" % np.abs(compress(coll2(g[3]))).max())
print("   separable (l3 = 0)   %.3e" % np.abs(compress(coll2(g3_sep))).max())

# basis for Delta: Hermitian C(M,3) x C(M,3) matrices
n3 = len(C3)
BAS = []
for a in range(n3):
    E = np.zeros((n3, n3), dtype=complex); E[a, a] = 1.0; BAS.append(E)
for a in range(n3):
    for b in range(a + 1, n3):
        E = np.zeros((n3, n3), dtype=complex); E[a, b] = 1.0; E[b, a] = 1.0; BAS.append(E)
        E = np.zeros((n3, n3), dtype=complex); E[a, b] = 1j; E[b, a] = -1j; BAS.append(E)
NB = len(BAS)

rowsT, rowsTr = [], []
for Z in BAS:
    T = ten3(Z)
    c = compress(coll2(T))
    rowsT.append(np.concatenate([np.real(c[np.triu_indices(d)]),
                                 np.imag(c[np.triu_indices(d, 1)])]))
    tr = np.einsum("pvyqwy->pvqw", T)
    rowsTr.append(np.concatenate([np.real(tr).ravel(), np.imag(tr).ravel()]))
A_tan = np.array(rowsT).T
A_tr = np.array(rowsTr).T
c0 = compress(coll2(g3_sep))
b_tan = -np.concatenate([np.real(c0[np.triu_indices(d)]),
                         np.imag(c0[np.triu_indices(d, 1)])])
A = np.vstack([A_tan, A_tr])
b = np.concatenate([b_tan, np.zeros(A_tr.shape[0])])
x, _, rk, sv = np.linalg.lstsq(A, b, rcond=None)
Delta = ten3(sum(xi * Z for xi, Z in zip(x, BAS)))
g3_fix = g3_sep + Delta
print()
print("least-squares projection: %d unknowns, %d tangency + %d trace equations, rank %d"
      % (NB, A_tan.shape[0], A_tr.shape[0], rk))
print("   |residual of the combined system|      %.3e" % np.abs(A @ x - b).max())
print("   tangency residual after correction     %.3e" % np.abs(compress(coll2(g3_fix))).max())
print("   partial trace of Delta (must be 0)     %.3e"
      % np.abs(np.einsum("pvyqwy->pvqw", Delta)).max())
print("   |Delta| / |gamma_3_exact|              %.4f"
      % (np.abs(Delta).max() / np.abs(g[3]).max()))
print()
print("accuracy against the exact gamma_3:")
for tag, g3 in (("separable", g3_sep), ("separable + Delta", g3_fix)):
    print("   %-20s |g3 - g3_exact| = %.5f   |Tr_3(g3)/(N-2) - g2| = %.3e"
          % (tag, np.abs(g3 - g[3]).max(),
             np.abs(np.einsum("pvxqwx->pvqw", g3) / (N - 2) - g[2]).max()))
