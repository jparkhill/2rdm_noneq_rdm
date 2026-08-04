#!/usr/bin/env python3
"""How does the exact EOM protect a saturated geminal?

The right place to see blocking is a state that already sits at the ceiling, not a
gamma_2 deformed out to it. Deforming rho -> rho + (s-n_g)|g><g| leaves the
N-representable set almost immediately, so the vector field there is not something
forward invariance has any opinion about.

The AGP condensate saturates Coleman's bound exactly, and that bound is a maximum over
all N-representable states. So along exact unitary dynamics n_g(t) <= c with equality
at t = 0, which forces

    ndot = 0   at the saturated geminal, for EVERY Hamiltonian.

That is a prediction, and if it holds the term-by-term decomposition shows how: which
sources vanish outright (a blocking factor) and which survive only to be cancelled by a
drain. Then the same decomposition under gamma_3 = 0 and lambda_3 = 0 shows which
approximation destroys the protection.
"""
import itertools
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from qdiss.symbolic import mkno
from qdiss.representability.dqg import contract_to_1rdm
from tools.eom2_coherent import ONE, TWO

UP3, LO3 = ("p", "v", "x"), ("q", "w", "y")
SEP3 = [(s, ",".join("".join(u) + "".join(l) for n, r, u, l in f) + "->pvxqwy",
         [(n, r) for n, r, u, l in f])
        for s, f in mkno.gamma_to_cumulants(UP3, LO3)
        if not any(n == "L" and r == 3 for n, r, u, l in f)]


def run(M, N, ntrial=6):
    CEIL = N * (M - N + 2) / M
    dim = 1 << M
    PAIRS = [(p, v) for p in range(M) for v in range(p + 1, M)]

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

    vac = np.zeros(dim, dtype=complex); vac[0] = 1.0
    psi = vac
    for _ in range(N // 2):
        acc = np.zeros_like(psi)
        for k in range(M // 2):
            acc = acc + ad(2 * k, ad(2 * k + 1, psi))
        psi = acc
    psi = psi / np.linalg.norm(psi)

    P1 = np.array([aa(p, psi) for p in range(M)])
    g1 = (P1.conj() @ P1.T).reshape(M, M)
    P2 = np.array([aa(v, aa(p, psi)) for p in range(M) for v in range(M)])
    g2 = (P2.conj() @ P2.T).reshape(M, M, M, M).transpose(0, 1, 2, 3)
    P3 = np.array([aa(x, aa(v, aa(p, psi)))
                   for p in range(M) for v in range(M) for x in range(M)])
    g3 = (P3.conj() @ P3.T).reshape(M, M, M, M, M, M)

    def to_pair(T):
        return np.array([[2.0 * T[p, v, q, w] for (q, w) in PAIRS] for (p, v) in PAIRS])

    rho = 0.5 * (to_pair(g2) + to_pair(g2).conj().T)
    w, U = np.linalg.eigh(rho)
    occ, gv = w[-1], U[:, -1]
    l2 = g2 - (np.einsum("pq,vw->pvqw", g1, g1) - np.einsum("pw,vq->pvqw", g1, g1))
    t = {("G", 1): g1, ("L", 2): l2}
    g3sep = np.zeros((M,) * 6, dtype=complex)
    for sign, sub, keys in SEP3:
        g3sep = g3sep + sign * np.einsum(sub, *[t[k] for k in keys], optimize=True)

    print("M=%d N=%d   ceiling %.4f   top geminal %.6f   (gap %.2e)   |g3 - g3sep| %.4f"
          % (M, N, CEIL, occ, abs(occ - CEIL), np.abs(g3 - g3sep).max()))

    eye = np.eye(M, dtype=complex)
    labels, ranks = [], []
    for tag, terms, _, _, _ in (("h", ONE, 0, 0, 0), ("V", TWO, 0, 0, 0)):
        for ti, (sign, facs) in enumerate(terms):
            labels.append("%s%d" % (tag, ti))
            ranks.append(max([r for n, r, u, l in facs if n == "G"] + [0]))

    def decompose(g3use):
        gam = {1: g1, 2: g2, 3: g3use}
        vals = []
        for terms, Hop, pref, hsub in ((ONE, h, 1.0, "rs"), (TWO, V, 0.25, "rtsu")):
            for sign, facs in terms:
                subs = [hsub]; ops = [Hop]
                for n, r, u, l in facs:
                    subs.append("".join(u) + "".join(l))
                    ops.append(eye if n == "D" else gam[r])
                T = sign * pref * np.einsum(",".join(subs) + "->pvqw", *ops, optimize=True)
                Tp = to_pair(1j * T); Tp = 0.5 * (Tp + Tp.conj().T)
                vals.append(float(np.real(gv.conj() @ Tp @ gv)))
        return np.array(vals)

    tot = {"exact": [], "g3=0": [], "l3=0": []}
    perterm = None
    for tr in range(ntrial):
        r = np.random.default_rng(900 + tr)
        A_ = r.normal(size=(M, M)) + 1j * r.normal(size=(M, M))
        h = 0.25 * (A_ + A_.conj().T)
        Bt = r.normal(size=(M,) * 4) + 1j * r.normal(size=(M,) * 4)
        V = Bt - np.transpose(Bt, (1, 0, 2, 3))
        V = V - np.transpose(V, (0, 1, 3, 2))
        V = 0.25 * (V + np.transpose(V.conj(), (2, 3, 0, 1)))
        ve = decompose(g3)
        tot["exact"].append(ve.sum())
        tot["g3=0"].append(decompose(np.zeros((M,) * 6, dtype=complex)).sum())
        tot["l3=0"].append(decompose(g3sep).sum())
        if tr == 0:
            perterm = ve
    for k in ("exact", "g3=0", "l3=0"):
        a = np.abs(np.array(tot[k]))
        print("   ndot at the saturated geminal, %-6s : max %.3e  mean %.3e"
              % (k, a.max(), a.mean()))
    print("   term-by-term for one H (exact gamma_3):")
    src = [(labels[i], ranks[i], perterm[i]) for i in range(len(perterm))
           if abs(perterm[i]) > 1e-10]
    zero = [labels[i] for i in range(len(perterm)) if abs(perterm[i]) <= 1e-10]
    for nm, rk, v in src:
        print("      %-4s carries gamma_%d : %+0.6f" % (nm, rk, v))
    print("      identically zero (%d terms): %s" % (len(zero), " ".join(zero)))
    print("      sum of nonzero terms: %+0.3e" % sum(v for _, _, v in src))
    print()


for M, N in ((8, 4), (10, 4)):
    run(M, N)
