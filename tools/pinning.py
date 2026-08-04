#!/usr/bin/env python3
"""Is the vanishing at saturation caused by PINNING of the higher RDM?

Rank 1 and rank 2 both protect their ceiling at second order, so no blocking factor
exists at either. The alternative mechanism: when an orbital saturates, the state
factorizes as psi = a+_p |chi> and gamma_2 is no longer free. Its blocks carrying the
index p are forced,

    gamma_2[p v, p w] = gamma_1[v, w]        for v, w != p
    gamma_2[p v, q w] = 0                    when p is in the bra pair but not the ket

so the terms of gamma_1-dot vanish because their gamma_2 factor is pinned, not because
anything multiplies them.

Two checks. First, measure the pinning residual at the saturated state and confirm it
grows with the gap. Second, and this is the real test, break the pinning by hand while
holding n_p = 1 fixed: perturb only the pinned blocks of gamma_2. If the terms come
back nonzero while the occupation is still exactly saturated, then pinning -- a
constraint on the ancilla -- is what does the protecting.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from tools.eom1 import one as ONE1, two as TWO1

M, N = 8, 4
dim = 1 << M


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
    return P1.conj() @ P1.T, (P2.conj() @ P2.T).reshape(M, M, M, M)


def terms(h, V, g1, g2, u):
    eye = np.eye(M, dtype=complex)
    gam = {1: g1, 2: g2}
    out = []
    for tl, Hop, pref, hsub in ((ONE1, h, 1.0, "rs"), (TWO1, V, 0.25, "rtsu")):
        for sign, facs in tl:
            subs = [hsub]; ops = [Hop]
            for n, r, up, lo in facs:
                subs.append("".join(up) + "".join(lo))
                ops.append(eye if n == "D" else gam[r])
            T = 1j * sign * pref * np.einsum(",".join(subs) + "->pq", *ops, optimize=True)
            T = 0.5 * (T + T.conj().T)
            out.append(float(np.real(u.conj() @ T @ u)))
    return np.array(out)


r = np.random.default_rng(808)
A_ = r.normal(size=(M, M)) + 1j * r.normal(size=(M, M))
h = 0.25 * (A_ + A_.conj().T)
Bt = r.normal(size=(M,) * 4) + 1j * r.normal(size=(M,) * 4)
V = Bt - np.transpose(Bt, (1, 0, 2, 3))
V = V - np.transpose(V, (0, 1, 3, 2))
V = 0.25 * (V + np.transpose(V.conj(), (2, 3, 0, 1)))

popc = np.array([bin(s).count("1") for s in range(dim)])
chi = np.where((popc == N - 1) & (np.arange(dim) & 1 == 0),
               r.normal(size=dim) + 1j * r.normal(size=dim), 0.0)
chi = chi / np.linalg.norm(chi)
sat = ad(0, chi) / np.linalg.norm(ad(0, chi))
mix = np.where(popc == N, r.normal(size=dim) + 1j * r.normal(size=dim), 0.0)
mix = mix / np.linalg.norm(mix)


def pin_residual(g1, g2, p=0):
    """how far gamma_2 is from the form forced by a saturated orbital p"""
    idx = [q for q in range(M) if q != p]
    a = np.array([[g2[p, v, p, w] - g1[v, w] for w in idx] for v in idx])
    b = np.array([[[g2[p, v, q, w] for w in idx] for q in idx] for v in idx])
    return np.abs(a).max(), np.abs(b).max()


print("=== pinning residual versus the gap, along the physical displacement ===")
print("%8s %11s | %14s %14s %12s" % ("t", "1-n_p", "|g2[pv,pw]-g1|", "|g2[pv,qw]|", "sum ndot"))
for t in (0.0, 0.01, 0.02, 0.04, 0.08):
    psi = np.cos(t) * sat + np.sin(t) * mix
    psi = psi / np.linalg.norm(psi)
    g1, g2 = rdms(psi)
    w, U = np.linalg.eigh(g1)
    ra, rb = pin_residual(g1, g2)
    print("%8.3f %11.2e | %14.3e %14.3e %12.3e"
          % (t, 1 - w[-1], ra, rb, terms(h, V, g1, g2, U[:, -1]).sum()))

print()
print("=== break the pinning by hand, holding n_p = 1 exactly ===")
print("perturb only the gamma_2 blocks carrying the saturated index; gamma_1 untouched")
print("%8s %11s | %14s %12s" % ("eps", "1-n_p", "pin residual", "sum ndot"))
g1, g2 = rdms(sat)
w, U = np.linalg.eigh(g1)
u = U[:, -1]
rp = np.random.default_rng(5)
Wp = rp.normal(size=(M, M)) + 1j * rp.normal(size=(M, M))
for eps in (0.0, 1e-3, 1e-2, 1e-1):
    g2m = g2.copy()
    for v in range(1, M):
        for q in range(1, M):
            d = eps * Wp[v, q]
            g2m[0, v, q, 0] += d; g2m[v, 0, 0, q] += np.conj(d)
            g2m[v, 0, q, 0] -= d; g2m[0, v, 0, q] -= np.conj(d)
    ra, rb = pin_residual(g1, g2m)
    print("%8.0e %11.2e | %14.3e %12.3e"
          % (eps, 1 - w[-1], max(ra, rb), terms(h, V, g1, g2m, u).sum()))
print()
print("n_p is exactly 1 in every row above, so any nonzero ndot is the ancilla")
print("constraint being violated, not the occupation moving off the boundary.")
