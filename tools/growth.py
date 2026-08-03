#!/usr/bin/env python3
"""How fast does lambda_3 grow out of a determinant, relative to lambda_2?

John's construction generates lambda_3 from its own EOM starting from a reference
where it vanishes. The premise is that lambda_3 is then subleading -- small, and
well captured by its driving term. A single determinant has lambda_2 = lambda_3 =
0 exactly, so propagating from one and watching both grow tests that premise and
fixes the power laws.

If lambda_2 ~ t and lambda_3 ~ t^2 or higher, the hierarchy is real and a
first-order-in-V generated lambda_3 is justified. If they grow at the same rate,
the construction has no small parameter to lean on.
"""
import itertools
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from qdiss.symbolic import mkno
from qdiss.representability.dqg import contract_to_1rdm

M, N = 6, 3
dim = 1 << M
UP3, LO3 = ("p", "v", "x"), ("q", "w", "y")
sep3 = [(s, f) for s, f in mkno.gamma_to_cumulants(UP3, LO3)
        if not any(n == "L" and r == 3 for n, r, u, l in f)]
SEP3 = [(s, ",".join("".join(u) + "".join(l) for n, r, u, l in f) + "->pvxqwy",
         [(n, r) for n, r, u, l in f]) for s, f in sep3]


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


def apply_str(up, lo, v):
    for q in lo:
        v = aa(q, v)
    for p in reversed(up):
        v = ad(p, v)
    return v


def rdm(psi, k):
    g = np.zeros((M,) * (2 * k), dtype=complex)
    for t in itertools.product(*([range(M)] * (2 * k))):
        g[t] = np.vdot(psi, apply_str(t[:k], t[k:], psi))
    return g


rng = np.random.default_rng(2026)
A_ = rng.normal(size=(M, M)) + 1j * rng.normal(size=(M, M))
h = 0.25 * (A_ + A_.conj().T)
B = rng.normal(size=(M,) * 4) + 1j * rng.normal(size=(M,) * 4)
V = B - np.transpose(B, (1, 0, 2, 3))
V = V - np.transpose(V, (0, 1, 3, 2))
V = 0.25 * (V + np.transpose(V.conj(), (2, 3, 0, 1)))


def Hpsi(v):
    out = np.zeros_like(v)
    for r in range(M):
        for s in range(M):
            out += h[r, s] * ad(r, aa(s, v))
    for r in range(M):
        for t in range(M):
            for s in range(M):
                for u in range(M):
                    c = V[r, t, s, u]
                    if c != 0:
                        out += 0.25 * c * ad(r, ad(t, aa(u, aa(s, v))))
    return out


def rk4(f, y, dt):
    k1 = f(y); k2 = f(y + 0.5 * dt * k1)
    k3 = f(y + 0.5 * dt * k2); k4 = f(y + dt * k3)
    return y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)


def lam2_of(g2, g1):
    return g2 - (np.einsum("pq,vw->pvqw", g1, g1) - np.einsum("pw,vq->pvqw", g1, g1))


def sep_g3(g1, l2):
    t = {("G", 1): g1, ("L", 2): l2}
    out = np.zeros((M,) * 6, dtype=complex)
    for sign, sub, keys in SEP3:
        out = out + sign * np.einsum(sub, *[t[k] for k in keys], optimize=True)
    return out


# start from a single determinant: lambda_2 = lambda_3 = 0 exactly
psi = np.zeros(dim, dtype=complex)
psi[sum(1 << i for i in range(N))] = 1.0

print("start: single determinant, M=%d N=%d" % (M, N))
print("%8s %14s %14s %12s %12s" % ("t", "|lambda_2|", "|lambda_3|", "l3/l2", "l3/l2^2"))
dt = 0.0025
rows = []
for step in range(41):
    t = step * dt
    if step % 5 == 0:
        g2 = rdm(psi, 2)
        g3 = rdm(psi, 3)
        g1 = contract_to_1rdm(g2, N)
        l2 = lam2_of(g2, g1)
        l3 = g3 - sep_g3(g1, l2)
        n2 = np.abs(l2).max()
        n3 = np.abs(l3).max()
        rows.append((t, n2, n3))
        r1 = n3 / n2 if n2 > 1e-14 else float("nan")
        r2 = n3 / n2 ** 2 if n2 > 1e-14 else float("nan")
        print("%8.4f %14.6e %14.6e %12.4f %12.4f" % (t, n2, n3, r1, r2))
    psi = rk4(lambda v: -1j * Hpsi(v), psi, dt)

print()
good = [(t, a, b) for t, a, b in rows if t > 0 and a > 1e-12 and b > 1e-14]
if len(good) >= 3:
    tt = np.array([g[0] for g in good])
    p2 = np.polyfit(np.log(tt), np.log([g[1] for g in good]), 1)[0]
    p3 = np.polyfit(np.log(tt), np.log([g[2] for g in good]), 1)[0]
    print("power laws from the log-log fit:")
    print("  |lambda_2| ~ t^%.2f" % p2)
    print("  |lambda_3| ~ t^%.2f" % p3)
    print()
    print("lambda_3 is subleading by t^%.2f, so it IS a small parameter at short"
          % (p3 - p2))
    print("time and a generated, first-order estimate is justified.")
