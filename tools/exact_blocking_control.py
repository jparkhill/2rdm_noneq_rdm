#!/usr/bin/env python3
"""Control for tools/exact_blocking.py: is the vanishing specific to saturation?

Every term of gamma_2-dot was zero at the saturated AGP geminal. That is only
meaningful if the same terms are nonzero elsewhere. Three controls:

  1. the SAME AGP state, at its lower geminals (occupations below the ceiling)
  2. a random N-representable state, at its top geminal (not saturated)
  3. AGP again but with the geminal rotated away from the condensate direction

If terms are nonzero in all three and zero only at saturation, the blocking is real.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from tools.eom2_coherent import ONE, TWO

M, N = 8, 4
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


def terms_at(h, V, g1, g2, g3, gv):
    eye = np.eye(M, dtype=complex)
    gam = {1: g1, 2: g2, 3: g3}
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


r = np.random.default_rng(4242)
A_ = r.normal(size=(M, M)) + 1j * r.normal(size=(M, M))
h = 0.25 * (A_ + A_.conj().T)
Bt = r.normal(size=(M,) * 4) + 1j * r.normal(size=(M,) * 4)
V = Bt - np.transpose(Bt, (1, 0, 2, 3))
V = V - np.transpose(V, (0, 1, 3, 2))
V = 0.25 * (V + np.transpose(V.conj(), (2, 3, 0, 1)))

vac = np.zeros(dim, dtype=complex); vac[0] = 1.0
agp = vac
for _ in range(N // 2):
    acc = np.zeros_like(agp)
    for k in range(M // 2):
        acc = acc + ad(2 * k, ad(2 * k + 1, agp))
    agp = acc
agp = agp / np.linalg.norm(agp)

popc = np.array([bin(s).count("1") for s in range(dim)])
rnd = np.where(popc == N, r.normal(size=dim) + 1j * r.normal(size=dim), 0.0)
rnd = rnd / np.linalg.norm(rnd)

print("ceiling = %.4f" % CEIL)
for name, psi in (("AGP condensate", agp), ("random N-rep state", rnd)):
    g1, g2, g3 = rdms(psi)
    rho = 0.5 * (to_pair(g2) + to_pair(g2).conj().T)
    w, U = np.linalg.eigh(rho)
    occ, Uo = w[::-1], U[:, ::-1]
    print()
    print("=== %s ===" % name)
    print("%6s %10s | %12s %12s %10s" % ("gem", "n_g", "sum ndot", "max|term|", "#nonzero"))
    for j in range(6):
        if occ[j] < 1e-9:
            continue
        v = terms_at(h, V, g1, g2, g3, Uo[:, j])
        print("%6d %10.6f | %12.3e %12.3e %10d"
              % (j, occ[j], v.sum(), np.abs(v).max(), int((np.abs(v) > 1e-10).sum())))

print()
print("=== AGP, geminal direction rotated away from the condensate ===")
g1, g2, g3 = rdms(agp)
rho = 0.5 * (to_pair(g2) + to_pair(g2).conj().T)
w, U = np.linalg.eigh(rho)
gtop, gnext = U[:, -1], U[:, -2]
print("%10s %10s | %12s %12s %10s" % ("mix angle", "<g|rho|g>", "sum ndot", "max|term|", "#nonzero"))
for th in (0.0, 0.05, 0.2, 0.5, 1.0):
    gv = np.cos(th) * gtop + np.sin(th) * gnext
    gv = gv / np.linalg.norm(gv)
    v = terms_at(h, V, g1, g2, g3, gv)
    print("%10.2f %10.6f | %12.3e %12.3e %10d"
          % (th, float(np.real(gv.conj() @ rho @ gv)), v.sum(), np.abs(v).max(),
             int((np.abs(v) > 1e-10).sum())))
