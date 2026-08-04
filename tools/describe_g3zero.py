#!/usr/bin/env python3
"""Exactly what gamma_3 = 0 is, and exactly where it enters the 2-RDM EOM.

Prints the full term list of the exact equation, marks which terms the approximation
deletes, checks the contraction sum rule that ties gamma_3 to gamma_2, and tests whether
what survives is a commutator flow.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from tools.eom2_coherent import ONE, TWO

M, N = 6, 3
dim = 1 << M
PAIRS = [(p, v) for p in range(M) for v in range(p + 1, M)]
SYM = {"D": "delta", "G": "gamma", "L": "lambda"}


def show(label, terms, pref):
    print("%s   (prefactor %s)" % (label, pref))
    for i, (s, facs) in enumerate(terms):
        parts = ["%s%d[%s;%s]" % (SYM[n], r, "".join(u), "".join(l)) for n, r, u, l in facs]
        rk = max([r for n, r, u, l in facs if n == "G"] + [0])
        print("   %2d  %+d  %-34s %s" % (i, s, " ".join(parts),
                                         "<-- DELETED by gamma_3 = 0" if rk == 3 else ""))
    print()


print("gamma_2[pv,qw] = <a+_p a+_v a_w a_q>,   gamma_2-dot = +i <[H, a+_p a+_v a_w a_q]>")
print("H = sum h^r_s a+_r a_s  +  (1/4) sum V^rt_su a+_r a+_t a_u a_s")
print()
show("one-body commutator <[a^r_s, a^pv_qw]>, contracted with h", ONE, "1")
show("two-body commutator <[a^rt_su, a^pv_qw]>, contracted with V", TWO, "1/4")
kept = [t for t in TWO if max([r for n, r, u, l in t[1] if n == "G"] + [0]) != 3]
print("exact EOM: %d + %d = %d terms.  gamma_3 = 0 keeps %d + %d = %d, all linear in gamma_2."
      % (len(ONE), len(TWO), len(ONE) + len(TWO), len(ONE), len(kept), len(ONE) + len(kept)))
print()


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
Q2 = np.array([aa(v, aa(p, psi)) for p in range(M) for v in range(M)])
g2 = (Q2.conj() @ Q2.T).reshape(M, M, M, M)
Q3 = np.array([aa(x, aa(v, aa(p, psi)))
               for p in range(M) for v in range(M) for x in range(M)])
g3 = (Q3.conj() @ Q3.T).reshape((M,) * 6)

print("=== the sum rule gamma_3 must satisfy ===")
c = np.einsum("pvxqwx->pvqw", g3)
print("   |sum_x gamma_3[pvx,qwx] - (N-2) gamma_2[pv,qw]| = %.3e   (exact gamma_3)"
      % np.abs(c - (N - 2) * g2).max())
print("   with gamma_3 = 0 the same expression is |(N-2) gamma_2| = %.4f"
      % np.abs((N - 2) * g2).max())
print("   so gamma_3 = 0 forces (N-2) gamma_2 = 0: consistent ONLY for N <= 2.")
print("   it is not merely a truncation, it contradicts the contraction that defines")
print("   gamma_2 as a partial trace of gamma_3.")
print()

A_ = r.normal(size=(M, M)) + 1j * r.normal(size=(M, M))
h = 0.25 * (A_ + A_.conj().T)
Bt = r.normal(size=(M,) * 4) + 1j * r.normal(size=(M,) * 4)
V = Bt - np.transpose(Bt, (1, 0, 2, 3))
V = V - np.transpose(V, (0, 1, 3, 2))
V = 0.25 * (V + np.transpose(V.conj(), (2, 3, 0, 1)))
g1 = np.array([[np.vdot(psi, ad(p, aa(q, psi))) for q in range(M)] for p in range(M)])


def to_pair(T):
    return np.array([[2.0 * T[p, v, q, w] for (q, w) in PAIRS] for (p, v) in PAIRS])


def g2dot(g3u):
    eye = np.eye(M, dtype=complex)
    gam = {1: g1, 2: g2, 3: g3u}
    tot = np.zeros((M,) * 4, dtype=complex)
    for tl, Hop, pf, hs in ((ONE, h, 1.0, "rs"), (TWO, V, 0.25, "rtsu")):
        for s, facs in tl:
            subs = [hs]; ops = [Hop]
            for n, rk, u, l in facs:
                subs.append("".join(u) + "".join(l))
                ops.append(eye if n == "D" else gam[rk])
            tot = tot + s * pf * np.einsum(",".join(subs) + "->pvqw", *ops, optimize=True)
    return 1j * tot


rho = 0.5 * (to_pair(g2) + to_pair(g2).conj().T)
lam, U = np.linalg.eigh(rho)
print("=== is the surviving gamma_3 = 0 flow a commutator? ===")
for name, g3u in (("gamma_3 = 0", np.zeros((M,) * 6, dtype=complex)), ("exact", g3)):
    X = to_pair(g2dot(g3u)); X = 0.5 * (X + X.conj().T)
    Xe = U.conj().T @ X @ U
    dg = np.abs(np.diag(Xe)).max()
    deg = np.abs(np.array([[Xe[i, j] for j in range(len(lam)) if abs(lam[i] - lam[j]) < 1e-9]
                           for i in range(len(lam))], dtype=object).sum())
    G = np.zeros_like(Xe)
    for i in range(len(lam)):
        for j in range(len(lam)):
            if abs(lam[j] - lam[i]) > 1e-9:
                G[i, j] = 1j * Xe[i, j] / (lam[j] - lam[i])
    Xr = U @ (-1j * (G @ np.diag(lam) - np.diag(lam) @ G)) @ U.conj().T
    print("   %-12s |tr Xdot| %.2e   max|diag in eigenbasis| %.3e   commutator residual %.3e"
          % (name, abs(np.trace(X)), dg, np.abs(Xr - X).max()))
print()
print("A zero diagonal in rho's own eigenbasis means no eigenvalue moves: the flow is")
print("isospectral, hence expressible as -i[X, rho], hence a rotation of the geminals")
print("with their occupations fixed.")
