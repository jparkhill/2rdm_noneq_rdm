#!/usr/bin/env python3
"""Which part of the EOM supplies positivity, and do candidate closures reproduce it?

THE SETUP. At rank k the propagated object is Mat_k(gamma_k)[I,J] = k! gamma_k[I+J] over sorted
k-tuples I, J of spin-orbitals; it is a Gram matrix of the k-fold annihilated state, hence
positive semidefinite. At a boundary point with kernel projector P_K, Nagumo's no-outward-flow
condition is P_K (d/dt Mat_k) P_K >= 0.

The rank-k EOM splits into 2k^2 terms carrying gamma_k (the drift) and 4k two-body terms carrying
gamma_{k+1} (what a closure supplies). The drift was shown to be of Lyapunov form A X + X A^dag,
whose kernel compression is identically zero. So ALL of the required non-negativity must be
supplied by the gamma_{k+1} terms, and this file measures that directly.

EXACT KERNELS. For a pure N-electron state, Mat_k(gamma_k) is the Gram matrix of the vectors
a_{j1}..a_{jk}|Psi>, which live in the (N-k)-electron sector, so rank <= C(M, N-k). Choosing
N - k small forces a large exact kernel: (M,N) = (6,3) at rank 2 gives rank <= C(6,1) = 6 out of
C(6,2) = 15, a 9-dimensional kernel; (M,N) = (6,4) at rank 3 gives rank <= 6 out of 20.

WHAT IS REPORTED, per rank:
  drift compression      P_K (drift) P_K, which must be zero to machine precision
  exact gamma_{k+1}      P_K (gamma_{k+1} terms, exact) P_K, whose eigenvalues must be >= 0
  closure gamma_{k+1}    the same with gamma_{k+1} replaced by a reconstruction. A negative
                         eigenvalue here is a closure that drives gamma_k out of the PSD cone at
                         that state -- a structural failure, independent of accuracy.

Closures tested. At rank 2: `separable`, i.e. lambda_3 = 0, gamma_3 built from the MKNO cumulant
expansion with the 3-cumulant dropped; and `muller2405`, the Mueller et al. arXiv:2405.05093
3-RDM reconstruction with the Lackner contraction correction. At rank 3: `separable`, i.e.
lambda_4 = 0. Term lists from qdiss.eom.wick, validated against qdiss.eom.fullci at 8.7e-16.
"""
import itertools
import sys
import numpy as np
sys.path.insert(0, "."); sys.path.insert(0, "tools")
from qdiss import muller2405
from qdiss.eom import fullci, wick
from qdiss.eom.hamiltonian import random_h, random_V
from qdiss.symbolic import mkno


def perm_sign(p):
    s, p = 1, list(p)
    for i in range(len(p)):
        for j in range(len(p) - 1 - i):
            if p[j] > p[j + 1]:
                p[j], p[j + 1] = p[j + 1], p[j]; s = -s
    return s


def mat_from_ten(T, k, M):
    C = list(itertools.combinations(range(M), k))
    f = float(np.math.factorial(k))
    return np.array([[f * T[I + J] for J in C] for I in C])


def partition(k):
    one, two = wick.terms(k)
    lo, hi = ([], []), ([], [])
    for which, lst in ((0, one), (1, two)):
        for s, f in lst:
            keys = [(n, r) for n, r, u, l in f]
            subs = ["".join(u) + "".join(l) for n, r, u, l in f]
            (hi if ("G", k + 1) in keys else lo)[which].append((s, subs, keys))
    return lo, hi


def apply_terms(terms, k, tab, h, V, M):
    out = "".join(wick.OPERATOR[k][0]) + "".join(wick.OPERATOR[k][1])
    tot = np.zeros((M,) * (2 * k), dtype=complex)
    for s, subs, keys in terms[0]:
        tot += s * np.einsum(",".join(["rs"] + subs) + "->" + out, h,
                             *[tab[key] for key in keys], optimize=True)
    for s, subs, keys in terms[1]:
        tot += s * 0.25 * np.einsum(",".join(["rtsu"] + subs) + "->" + out, V,
                                    *[tab[key] for key in keys], optimize=True)
    return 1j * tot


def separable(k, tab, M):
    """gamma_k from the MKNO cumulant expansion with the top cumulant lambda_k dropped."""
    up = tuple("pvxz"[:k]); lo = tuple("qwya"[:k])
    out = np.zeros((M,) * (2 * k), dtype=complex)
    for s, facs in mkno.gamma_to_cumulants(up, lo):
        if any(n == "L" and r == k for n, r, u, l in facs):
            continue
        sub = ",".join("".join(u) + "".join(l) for n, r, u, l in facs) + "->" + "".join(up + lo)
        out += s * np.einsum(sub, *[tab[(n, r)] for n, r, u, l in facs], optimize=True)
    return out


def cumulant_table(g, M):
    t = {("D", 1): np.eye(M, dtype=complex)}
    for r, v in g.items():
        t[("G", r)] = v
    for r in sorted(g)[1:]:
        t[("L", r)] = g[r] - separable(r, t, M)
    return t


for k, (M, N) in ((2, (6, 3)), (3, (6, 4))):
    rng = np.random.default_rng(3)
    h, V = random_h(M, rng), random_V(M, rng)
    dets = [d for d in range(1 << M) if bin(d).count("1") == N]
    Hm = fullci.build_H(h, V, M)
    es, vec = np.linalg.eigh(Hm[np.ix_(dets, dets)])
    psi = np.zeros(1 << M, dtype=complex); psi[dets] = vec[:, 0]
    g = {r: fullci.rdm(psi, M, r) for r in range(1, k + 2)}
    tab = cumulant_table(g, M)

    X = mat_from_ten(g[k], k, M)
    ev, U = np.linalg.eigh(0.5 * (X + X.conj().T))
    ker = U[:, ev < 1e-9 * max(ev.max(), 1.0)]
    print("=== rank %d,  (M,N) = (%d,%d) ===" % (k, M, N))
    print("    Mat_%d is %dx%d;  rank bound C(M,N-k) = %d;  numerical rank %d;  kernel dim %d"
          % (k, X.shape[0], X.shape[0],
             len(list(itertools.combinations(range(M), N - k))),
             int((ev > 1e-9).sum()), ker.shape[1]))

    lo, hi = partition(k)

    def compress(T):
        Y = mat_from_ten(T, k, M)
        B = ker.conj().T @ (0.5 * (Y + Y.conj().T)) @ ker
        return np.linalg.eigvalsh(B)

    d = compress(apply_terms(lo, k, tab, h, V, M))
    print("    drift  P_K . P_K :  max |eig| = %.2e   (must be zero: Lyapunov form)"
          % np.abs(d).max())

    cands = [("exact", g[k + 1]), ("separable", separable(k + 1, tab, M))]
    if k == 2:
        cands.append(("muller2405", muller2405.gamma3(g[1], g[2], N)))
    print("    %-12s | %12s %12s | %s" % ("gamma_%d" % (k + 1), "min eig", "max eig", "verdict"))
    for tag, gk1 in cands:
        t2 = dict(tab); t2[("G", k + 1)] = gk1
        e = compress(apply_terms(hi, k, t2, h, V, M))
        print("    %-12s | %+12.4e %+12.4e | %s"
              % (tag, e.min(), e.max(), "PSD" if e.min() > -1e-9 else "VIOLATES"))
    print()
