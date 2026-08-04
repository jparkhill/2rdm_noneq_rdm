#!/usr/bin/env python3
"""Does the flow PRESERVE pinning? The consequence test, done correctly.

tools/pinning_gamma3.py contracted gamma_2-dot down to gamma_1-dot and found all
closures agreeing at 1e-16. That measures nothing: gamma_1-dot depends on h, V,
gamma_1 and gamma_2 only, never on gamma_3, so the contraction is gamma_3-blind by
construction. That test is void.

The real consequence is whether the pinned structure SURVIVES the flow. At n_p = 1
representability forces

    A[v,w]   = gamma_2[p v, p w] - gamma_1[v, w]  = 0
    B[v,q,w] = gamma_2[p v, q w]                  = 0

so a flow that keeps the state representable must hold these at zero, i.e.

    Adot[v,w]   = gamma_2-dot[p v, p w] - gamma_1-dot[v, w] = 0
    Bdot[v,q,w] = gamma_2-dot[p v, q w]                     = 0

where gamma_1-dot comes from the rank-1 EOM (gamma_2 only) and gamma_2-dot from the
rank-2 EOM (gamma_3 included). Now gamma_3 enters, and a closure that mis-pins gamma_3
should break the structure immediately while n_p sits at exactly 1 -- which is the
state leaving the N-representable set without the occupation ever moving.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from qdiss.symbolic import mkno
from tools.eom1 import one as ONE1, two as TWO1
from tools.eom2_coherent import ONE, TWO

M, N = 8, 4
dim = 1 << M
P = 0
OTH = [q for q in range(M) if q != P]
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


def sep(g1, g2):
    l2 = g2 - (np.einsum("pq,vw->pvqw", g1, g1) - np.einsum("pw,vq->pvqw", g1, g1))
    t = {("G", 1): g1, ("L", 2): l2}
    out = np.zeros((M,) * 6, dtype=complex)
    for s, sub, keys in SEP3:
        out = out + s * np.einsum(sub, *[t[k] for k in keys], optimize=True)
    return out


def evolve(tl_h, tl_v, out_idx, h, V, gam):
    eye = np.eye(M, dtype=complex)
    shape = (M,) * len(out_idx)
    tot = np.zeros(shape, dtype=complex)
    for tl, Hop, pref, hsub in ((tl_h, h, 1.0, "rs"), (tl_v, V, 0.25, "rtsu")):
        for sign, facs in tl:
            subs = [hsub]; ops = [Hop]
            for n, r, u, l in facs:
                subs.append("".join(u) + "".join(l))
                ops.append(eye if n == "D" else gam[r])
            tot = tot + sign * pref * np.einsum(",".join(subs) + "->" + out_idx, *ops,
                                                optimize=True)
    return 1j * tot


r = np.random.default_rng(808)
popc = np.array([bin(s).count("1") for s in range(dim)])
chi = np.where((popc == N - 1) & (np.arange(dim) & 1 == 0),
               r.normal(size=dim) + 1j * r.normal(size=dim), 0.0)
chi = chi / np.linalg.norm(chi)
sat = ad(P, chi) / np.linalg.norm(ad(P, chi))
g1, g2, g3 = rdms(sat)
g3s = sep(g1, g2)
g3z = np.zeros((M,) * 6, dtype=complex)
print("saturated orbital occupation %.15f" % np.linalg.eigvalsh(g1)[-1])
print()
print("static pinning residual of gamma_3 : exact %.2e   lambda_3=0 %.2e   gamma_3=0 %.2e"
      % (max(abs(g3[P, v, x, P, w, y] - g2[v, x, w, y])
             for v in OTH for x in OTH for w in OTH for y in OTH),
         max(abs(g3s[P, v, x, P, w, y] - g2[v, x, w, y])
             for v in OTH for x in OTH for w in OTH for y in OTH),
         max(abs(g3z[P, v, x, P, w, y] - g2[v, x, w, y])
             for v in OTH for x in OTH for w in OTH for y in OTH)))
print()
print("RATE at which the flow breaks pinning, with n_p pinned at 1 throughout")
print("%4s | %23s | %23s | %23s"
      % ("H", "exact  |Adot|  |Bdot|", "lambda_3=0", "gamma_3=0"))
acc = {"exact": [], "l3": [], "g3": []}
for tr in range(6):
    rr = np.random.default_rng(2200 + tr)
    A_ = rr.normal(size=(M, M)) + 1j * rr.normal(size=(M, M))
    h = 0.25 * (A_ + A_.conj().T)
    Bt = rr.normal(size=(M,) * 4) + 1j * rr.normal(size=(M,) * 4)
    V = Bt - np.transpose(Bt, (1, 0, 2, 3))
    V = V - np.transpose(V, (0, 1, 3, 2))
    V = 0.25 * (V + np.transpose(V.conj(), (2, 3, 0, 1)))
    d1 = evolve(ONE1, TWO1, "pq", h, V, {1: g1, 2: g2})
    cells = []
    for key, g3u in (("exact", g3), ("l3", g3s), ("g3", g3z)):
        d2 = evolve(ONE, TWO, "pvqw", h, V, {1: g1, 2: g2, 3: g3u})
        Ad = max(abs(d2[P, v, P, w] - d1[v, w]) for v in OTH for w in OTH)
        Bd = max(abs(d2[P, v, q, w]) for v in OTH for q in OTH for w in OTH)
        acc[key].append((Ad, Bd)); cells.append((Ad, Bd))
    print("%4d | %10.3e %10.3e   | %10.3e %10.3e   | %10.3e %10.3e"
          % (tr, cells[0][0], cells[0][1], cells[1][0], cells[1][1],
             cells[2][0], cells[2][1]))
print()
for k, nm in (("exact", "exact"), ("l3", "lambda_3=0"), ("g3", "gamma_3=0")):
    a = np.array(acc[k])
    print("   %-11s worst |Adot| %.3e   worst |Bdot| %.3e" % (nm, a[:, 0].max(), a[:, 1].max()))
