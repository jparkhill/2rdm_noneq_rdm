#!/usr/bin/env python3
"""Does the lambda_3 = 0 reconstruction respect pinning, and does violating it break Pauli?

Rank-2 pinning at Coleman's ceiling needs the extremality theorem. This sidesteps that:
use rank-1 saturation, which factorizes elementarily as psi = a+_p |chi>, but ask what
it pins in the rank-2 ancilla gamma_3:

    gamma_3[p v x, p w y] = gamma_2[v x, w y]        v,x,w,y != p
    gamma_3[p v x, q w y] = 0                        p not in {q, w, y}

Order of business:
  0. verify the exact gamma_3 satisfies both, which validates the index conventions
  1. the rank-1 control: does the HF closure gamma_2 = gamma_1 wedge gamma_1 pin?
  2. does the lambda_3 = 0 reconstruction of gamma_3 pin?
  3. the consequence: propagate gamma_2 under each closure, contract gamma_2-dot back to
     gamma_1-dot, and read ndot at the saturated orbital. Exact must give zero. A closure
     that violates pinning should give nonzero, i.e. it pushes n_p above 1 and off the
     N-representable set.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from qdiss.symbolic import mkno
from qdiss.representability.dqg import contract_to_1rdm
from tools.eom2_coherent import ONE, TWO

M, N = 8, 4
dim = 1 << M
P = 0
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
    P1 = np.array([aa(p, psi) for p in range(M)])
    P2 = np.array([aa(v, aa(p, psi)) for p in range(M) for v in range(M)])
    P3 = np.array([aa(x, aa(v, aa(p, psi)))
                   for p in range(M) for v in range(M) for x in range(M)])
    return (P1.conj() @ P1.T, (P2.conj() @ P2.T).reshape(M, M, M, M),
            (P3.conj() @ P3.T).reshape((M,) * 6))


OTH = [q for q in range(M) if q != P]


def pin2(g1, g2):
    a = max(abs(g2[P, v, P, w] - g1[v, w]) for v in OTH for w in OTH)
    b = max(abs(g2[P, v, q, w]) for v in OTH for q in OTH for w in OTH)
    return a, b


def pin3(g2, g3):
    a = max(abs(g3[P, v, x, P, w, y] - g2[v, x, w, y])
            for v in OTH for x in OTH for w in OTH for y in OTH)
    b = max(abs(g3[P, v, x, q, w, y])
            for v in OTH for x in OTH for q in OTH for w in OTH for y in OTH)
    return a, b


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


r = np.random.default_rng(808)
popc = np.array([bin(s).count("1") for s in range(dim)])
chi = np.where((popc == N - 1) & (np.arange(dim) & 1 == 0),
               r.normal(size=dim) + 1j * r.normal(size=dim), 0.0)
chi = chi / np.linalg.norm(chi)
sat = ad(P, chi) / np.linalg.norm(ad(P, chi))
g1, g2, g3 = rdms(sat)
w, U = np.linalg.eigh(g1)
u = U[:, -1]
print("saturated orbital occupation: %.15f" % w[-1])
print()

print("=== 0. exact gamma_3, to validate the index conventions ===")
a, b = pin3(g2, g3)
print("   |g3[pvx,pwy] - g2[vx,wy]| = %.3e     |g3[pvx,qwy]| = %.3e" % (a, b))
print()

print("=== 1. rank-1 control: does the HF closure pin gamma_2? ===")
g2hf = np.einsum("pq,vw->pvqw", g1, g1) - np.einsum("pw,vq->pvqw", g1, g1)
a, b = pin2(g1, g2)
print("   exact gamma_2      : %.3e   %.3e" % (a, b))
a, b = pin2(g1, g2hf)
print("   HF gamma_1^gamma_1 : %.3e   %.3e   |g2hf - g2| = %.4f"
      % (a, b, np.abs(g2hf - g2).max()))
print()

print("=== 2. does lambda_3 = 0 pin gamma_3? ===")
g3s = sep(g1, g2)
a, b = pin3(g2, g3s)
print("   lambda_3 = 0       : %.3e   %.3e   |g3sep - g3| = %.4f"
      % (a, b, np.abs(g3s - g3).max()))
a, b = pin3(g2, np.zeros((M,) * 6, dtype=complex))
print("   gamma_3 = 0        : %.3e   %.3e" % (a, b))
print()

print("=== 3. consequence: ndot at the saturated orbital, from gamma_2-dot ===")
print("%6s | %12s %12s %12s" % ("H", "exact", "lambda_3=0", "gamma_3=0"))
res = {"exact": [], "l3": [], "g3": []}
for tr in range(6):
    rr = np.random.default_rng(2200 + tr)
    A_ = rr.normal(size=(M, M)) + 1j * rr.normal(size=(M, M))
    h = 0.25 * (A_ + A_.conj().T)
    Bt = rr.normal(size=(M,) * 4) + 1j * rr.normal(size=(M,) * 4)
    V = Bt - np.transpose(Bt, (1, 0, 2, 3))
    V = V - np.transpose(V, (0, 1, 3, 2))
    V = 0.25 * (V + np.transpose(V.conj(), (2, 3, 0, 1)))
    row = []
    for key, g3u in (("exact", g3), ("l3", g3s), ("g3", np.zeros((M,) * 6, dtype=complex))):
        d2 = g2dot(h, V, g1, g2, g3u)
        d1 = contract_to_1rdm(d2, N)
        d1 = 0.5 * (d1 + d1.conj().T)
        v = float(np.real(u.conj() @ d1 @ u))
        res[key].append(v); row.append(v)
    print("%6d | %+12.3e %+12.3e %+12.3e" % (tr, row[0], row[1], row[2]))
print()
for k, nm in (("exact", "exact"), ("l3", "lambda_3=0"), ("g3", "gamma_3=0")):
    a = np.array(res[k])
    print("   %-11s max |ndot| %.3e   %d of %d push n_p ABOVE 1"
          % (nm, np.abs(a).max(), int((a > 1e-9).sum()), len(a)))
