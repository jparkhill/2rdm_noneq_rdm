#!/usr/bin/env python3
"""Is the CLOSURE's own flow CP-generated, and what is the nearest CP closure?

Every CP tool in this repo fits a free GKSL generator to the EXACT gamma_2-dot
(sdp_cpfit.py line 104 sets target from the exact gamma_3). That answers whether the
true dynamics is CP-representable. It never asks the constructive question: is the flow
produced by the MKNO closure itself CP-generated, and if not, what is the closest flow
that is?

At a fixed rho the question is convex. Ask for the nearest point of the GKSL cone image,

    L(rho) = -i[K, rho] + sum_ab c_ab ( F_a rho F_b^dag - 1/2 {F_b^dag F_a, rho} )

over K Hermitian and c >= 0, minimizing || L(rho) - target ||_F. Zero residual means the
closure's instantaneous flow already lies in the cone and is D-preserving. Nonzero means
it does not, and the minimizer defines a CP-projected closure -- derived by projection
rather than by inserting a factor.

Then the question that decides whether the projection is worth anything: does the
projected flow still RELAX? d(tr rho^2)/dt is zero iff the flow is isospectral. A CP
projection that lands on an isospectral generator has bought positivity at the price of
frozen occupations, which is gamma_3 = 0 again by another route.
"""
import itertools
import sys
import numpy as np
import cvxpy as cp
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from qdiss.symbolic import mkno
from tools.eom2_coherent import ONE, TWO

M, N = 6, 3
dim = 1 << M
PAIRS = [(p, v) for p in range(M) for v in range(p + 1, M)]
NP = len(PAIRS)
UP3, LO3 = ("p", "v", "x"), ("q", "w", "y")
SEP3 = [(s, ",".join("".join(u) + "".join(l) for n, r, u, l in f) + "->pvxqwy",
         [(n, r) for n, r, u, l in f])
        for s, f in mkno.gamma_to_cumulants(UP3, LO3)
        if not any(n == "L" and r == 3 for n, r, u, l in f)]


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


def rdms(psi):
    Q1 = np.array([aa(p, psi) for p in range(M)])
    Q2 = np.array([aa(v, aa(p, psi)) for p in range(M) for v in range(M)])
    Q3 = np.array([aa(x, aa(v, aa(p, psi)))
                   for p in range(M) for v in range(M) for x in range(M)])
    return (Q1.conj() @ Q1.T, (Q2.conj() @ Q2.T).reshape(M, M, M, M),
            (Q3.conj() @ Q3.T).reshape((M,) * 6))


def to_pair(T):
    return np.array([[2.0 * T[p, v, q, w] for (q, w) in PAIRS] for (p, v) in PAIRS])


def sep(g1, g2):
    l2 = g2 - (np.einsum("pq,vw->pvqw", g1, g1) - np.einsum("pw,vq->pvqw", g1, g1))
    t = {("G", 1): g1, ("L", 2): l2}
    out = np.zeros((M,) * 6, dtype=complex)
    for s, sub, keys in SEP3:
        out = out + s * np.einsum(sub, *[t[k] for k in keys], optimize=True)
    return out


def g2dot(h, V, g1, g2, g3):
    eye = np.eye(M, dtype=complex)
    gam = {1: g1, 2: g2, 3: g3}
    tot = np.zeros((M,) * 4, dtype=complex)
    for tl, Hop, pref, hsub in ((ONE, h, 1.0, "rs"), (TWO, V, 0.25, "rtsu")):
        for sign, facs in tl:
            subs = [hsub]; ops = [Hop]
            for n, r, u, l in facs:
                subs.append("".join(u) + "".join(l))
                ops.append(eye if n == "D" else gam[r])
            tot = tot + sign * pref * np.einsum(",".join(subs) + "->pvqw", *ops,
                                                optimize=True)
    return 1j * tot


r = np.random.default_rng(11)
A_ = r.normal(size=(M, M)) + 1j * r.normal(size=(M, M))
h = 0.25 * (A_ + A_.conj().T)
Bt = r.normal(size=(M,) * 4) + 1j * r.normal(size=(M,) * 4)
V = Bt - np.transpose(Bt, (1, 0, 2, 3))
V = V - np.transpose(V, (0, 1, 3, 2))
V = 0.25 * (V + np.transpose(V.conj(), (2, 3, 0, 1)))
popc = np.array([bin(s).count("1") for s in range(dim)])
psi = np.where(popc == N, r.normal(size=dim) + 1j * r.normal(size=dim), 0.0)
psi = psi / np.linalg.norm(psi)
g1, g2, g3 = rdms(psi)
rho = 0.5 * (to_pair(g2) + to_pair(g2).conj().T)
ev = np.linalg.eigvalsh(rho)
print("pair dim %d   rho eigenvalues %.3e .. %.3e   tr %.4f" % (NP, ev[0], ev[-1], ev.sum()))
print()

FB = []
for i in range(NP):
    for j in range(NP):
        E = np.zeros((NP, NP), dtype=complex); E[i, j] = 1.0
        FB.append(E)
FB = [F - np.trace(F) / NP * np.eye(NP) for F in FB]
FB = [F for F in FB if np.abs(F).max() > 1e-12]
NB = len(FB)


def project(target, nb):
    basis = FB[:nb]
    K = cp.Variable((NP, NP), hermitian=True)
    c = cp.Variable((nb, nb), hermitian=True)
    expr = -1j * (K @ rho - rho @ K)
    for a in range(nb):
        for b in range(nb):
            Fa, Fb = basis[a], basis[b]
            expr = expr + c[a, b] * (Fa @ rho @ Fb.conj().T
                                     - 0.5 * (Fb.conj().T @ Fa @ rho
                                              + rho @ Fb.conj().T @ Fa))
    prob = cp.Problem(cp.Minimize(cp.norm(expr - target, "fro")), [c >> 0])
    prob.solve(solver=cp.SCS, max_iters=8000, eps=1e-7)
    return prob.value, expr.value


NBUSE = int(sys.argv[1]) if len(sys.argv) > 1 else 40
print("jump basis: %d of %d" % (NBUSE, NB))
print("%-12s %12s %12s %12s %12s"
      % ("closure", "|target|", "CP resid", "rel resid", "d(tr rho^2)"))
for name, g3u in (("exact", g3), ("lambda_3=0", sep(g1, g2)),
                  ("gamma_3=0", np.zeros((M,) * 6, dtype=complex))):
    tgt = to_pair(g2dot(h, V, g1, g2, g3u))
    tgt = 0.5 * (tgt + tgt.conj().T)
    nt = np.linalg.norm(tgt)
    res, Lval = project(tgt, NBUSE)
    dpur = 2 * float(np.real(np.trace(rho @ tgt)))
    print("%-12s %12.4f %12.4f %12.4f %12.4e"
          % (name, nt, res, res / nt, dpur))

print()
print("d(tr rho^2) is zero iff the flow is isospectral and cannot relax occupations.")
