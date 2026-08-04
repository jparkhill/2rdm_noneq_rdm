#!/usr/bin/env python3
"""Does the blocking survive under closure at states NEAR saturation?

At exact saturation every term vanishes, for the exact gamma_3 and under lambda_3 = 0
alike. That is the boundary itself. Forward invariance needs slightly more: approaching
the boundary from inside the feasible set, the flow must not point outward.

So interpolate psi(t) = cos(t) * AGP + sin(t) * random, which sweeps the top geminal
down from the ceiling through physical, N-representable states -- every one of them a
genuine wavefunction, unlike the deformed gamma_2 used earlier. At each, measure ndot
of the top geminal under the exact gamma_3, under gamma_3 = 0, and under lambda_3 = 0,
over several Hamiltonians.

Blocking under a closure means ndot <= 0 as the ceiling is approached. A positive ndot
near saturation is the closure pushing the state out of the representable set.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from qdiss.symbolic import mkno
from tools.eom2_coherent import ONE, TWO

M, N = 10, 4
CEIL = N * (M - N + 2) / M
dim = 1 << M
PAIRS = [(p, v) for p in range(M) for v in range(p + 1, M)]
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
    g1 = P1.conj() @ P1.T
    P2 = np.array([aa(v, aa(p, psi)) for p in range(M) for v in range(M)])
    g2 = (P2.conj() @ P2.T).reshape(M, M, M, M)
    P3 = np.array([aa(x, aa(v, aa(p, psi)))
                   for p in range(M) for v in range(M) for x in range(M)])
    g3 = (P3.conj() @ P3.T).reshape((M,) * 6)
    return g1, g2, g3


def to_pair(T):
    return np.array([[2.0 * T[p, v, q, w] for (q, w) in PAIRS] for (p, v) in PAIRS])


def sep(g1, g2):
    l2 = g2 - (np.einsum("pq,vw->pvqw", g1, g1) - np.einsum("pw,vq->pvqw", g1, g1))
    t = {("G", 1): g1, ("L", 2): l2}
    out = np.zeros((M,) * 6, dtype=complex)
    for s, sub, keys in SEP3:
        out = out + s * np.einsum(sub, *[t[k] for k in keys], optimize=True)
    return out


def ndot(h, V, g1, g2, g3, gv):
    eye = np.eye(M, dtype=complex)
    gam = {1: g1, 2: g2, 3: g3}
    tot = 0.0
    for terms, Hop, pref, hsub in ((ONE, h, 1.0, "rs"), (TWO, V, 0.25, "rtsu")):
        for sign, facs in terms:
            subs = [hsub]; ops = [Hop]
            for n, r, u, l in facs:
                subs.append("".join(u) + "".join(l))
                ops.append(eye if n == "D" else gam[r])
            T = sign * pref * np.einsum(",".join(subs) + "->pvqw", *ops, optimize=True)
            Tp = to_pair(1j * T); Tp = 0.5 * (Tp + Tp.conj().T)
            tot += float(np.real(gv.conj() @ Tp @ gv))
    return tot


vac = np.zeros(dim, dtype=complex); vac[0] = 1.0
agp = vac
for _ in range(N // 2):
    acc = np.zeros_like(agp)
    for k in range(M // 2):
        acc = acc + ad(2 * k, ad(2 * k + 1, agp))
    agp = acc
agp = agp / np.linalg.norm(agp)
popc = np.array([bin(s).count("1") for s in range(dim)])

HS = []
for tr in range(5):
    r = np.random.default_rng(77 + tr)
    A_ = r.normal(size=(M, M)) + 1j * r.normal(size=(M, M))
    h = 0.25 * (A_ + A_.conj().T)
    Bt = r.normal(size=(M,) * 4) + 1j * r.normal(size=(M,) * 4)
    V = Bt - np.transpose(Bt, (1, 0, 2, 3))
    V = V - np.transpose(V, (0, 1, 3, 2))
    HS.append((h, 0.25 * (V + np.transpose(V.conj(), (2, 3, 0, 1)))))

rr = np.random.default_rng(31337)
mix = np.where(popc == N, rr.normal(size=dim) + 1j * rr.normal(size=dim), 0.0)
mix = mix / np.linalg.norm(mix)

print("M=%d N=%d  ceiling %.4f   (%d Hamiltonians per point)" % (M, N, CEIL, len(HS)))
print()
print("%8s %10s %10s | %22s %22s"
      % ("t", "n_top", "ceil-n", "exact: max ndot", "lambda_3=0: max ndot"))
worst = {"exact": -9e9, "l3": -9e9}
for t in (0.0, 0.02, 0.05, 0.10, 0.20, 0.35, 0.55, 0.80):
    psi = np.cos(t) * agp + np.sin(t) * mix
    psi = psi / np.linalg.norm(psi)
    g1, g2, g3 = rdms(psi)
    rho = 0.5 * (to_pair(g2) + to_pair(g2).conj().T)
    w, U = np.linalg.eigh(rho)
    occ, gv = w[-1], U[:, -1]
    g3s = sep(g1, g2)
    ve = [ndot(h, V, g1, g2, g3, gv) for h, V in HS]
    vl = [ndot(h, V, g1, g2, g3s, gv) for h, V in HS]
    worst["exact"] = max(worst["exact"], max(ve))
    worst["l3"] = max(worst["l3"], max(vl))
    print("%8.2f %10.6f %10.2e | %+10.4f (%d/%d out) %+10.4f (%d/%d out)"
          % (t, occ, CEIL - occ, max(ve), sum(1 for x in ve if x > 1e-9), len(ve),
             max(vl), sum(1 for x in vl if x > 1e-9), len(vl)))
print()
print("worst outward flow at the TOP geminal:  exact %+0.4f   lambda_3=0 %+0.4f"
      % (worst["exact"], worst["l3"]))
print()
print("note: away from saturation the top geminal is interior, so ndot > 0 is")
print("allowed and expected. What matters is the trend as ceil-n -> 0.")
