#!/usr/bin/env python3
"""Nagumo-correct test: can any fixed-coefficient gamma_3 block at the ceiling?

tools/fit_gamma3.py asked whether some (a,b) makes ndot vanish AT the ceiling. That
was the wrong demand twice over. Its residual was reported relative to a base term
that is identically zero, so the printed verdict was 0/0; and forward invariance does
not need ndot = 0, only ndot <= 0. A geminal is allowed to be pushed inward at its
ceiling, just not outward.

The base result stands and is the reason the question sharpens. Writing

    ndot(ceiling) = F0(H) + a*F1(H) + b*F2(H),

F0 is the gamma_3-free part, and it is zero to machine precision for every H. So
gamma_3 is the ENTIRE driver of geminal occupation change; the gamma_2-only terms are
purely rotational.

The real question is then a half-plane feasibility problem. Blocking for all H means

    a*F1(H) + b*F2(H) <= 0   for every H,

i.e. all the direction vectors (F1(H), F2(H)) lie in a common half-plane through the
origin. In two dimensions that is decidable exactly: sort the directions by angle and
look at the largest angular gap. A gap exceeding pi means a valid (a,b) exists; no gap
exceeding pi means the directions positively span the plane and no fixed-coefficient
reconstruction can block.
"""
import itertools
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from qdiss.symbolic import mkno
from qdiss.representability.dqg import contract_to_1rdm
from tools.eom2_coherent import ONE, TWO

M, N = 6, 3
CEIL = N * (M - N + 2) / M
dim = 1 << M
PAIRS = [(p, v) for p in range(M) for v in range(p + 1, M)]
UP3, LO3 = ("p", "v", "x"), ("q", "w", "y")
ORB = {"g1g1g1": [], "g1l2": []}
for sign, facs in mkno.gamma_to_cumulants(UP3, LO3):
    if any(n == "L" and r == 3 for n, r, u, l in facs):
        continue
    key = "g1l2" if any(n == "L" and r == 2 for n, r, u, l in facs) else "g1g1g1"
    ORB[key].append((sign, ",".join("".join(u) + "".join(l) for n, r, u, l in facs)
                     + "->pvxqwy", [(n, r) for n, r, u, l in facs]))


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


def orbit(key, g1, g2):
    l2 = g2 - (np.einsum("pq,vw->pvqw", g1, g1) - np.einsum("pw,vq->pvqw", g1, g1))
    t = {("G", 1): g1, ("L", 2): l2}
    out = np.zeros((M,) * 6, dtype=complex)
    for sign, sub, keys in ORB[key]:
        out = out + sign * np.einsum(sub, *[t[k] for k in keys], optimize=True)
    return out


def to_pair(T):
    return np.array([[2.0 * T[p, v, q, w] for (q, w) in PAIRS] for (p, v) in PAIRS])


def from_pair(A):
    T = np.zeros((M,) * 4, dtype=complex)
    for i, (p, v) in enumerate(PAIRS):
        for j, (q, w) in enumerate(PAIRS):
            x = 0.5 * A[i, j]
            T[p, v, q, w] = x; T[v, p, q, w] = -x
            T[p, v, w, q] = -x; T[v, p, w, q] = x
    return T


def pieces(h, V, g1, g2, g3, gvec):
    eye = np.eye(M, dtype=complex)
    gam = {1: g1, 2: g2, 3: g3}
    base = 0.0; g3part = 0.0
    for terms, Hop, pref, hsub in ((ONE, h, 1.0, "rs"), (TWO, V, 0.25, "rtsu")):
        for sign, facs in terms:
            subs = [hsub]; ops = [Hop]; has3 = False
            for n, r, u, l in facs:
                subs.append("".join(u) + "".join(l))
                ops.append(eye if n == "D" else gam[r])
                if n == "G" and r == 3:
                    has3 = True
            T = sign * pref * np.einsum(",".join(subs) + "->pvqw", *ops, optimize=True)
            Tp = to_pair(1j * T); Tp = 0.5 * (Tp + Tp.conj().T)
            val = float(np.real(gvec.conj() @ Tp @ gvec))
            if has3: g3part += val
            else: base += val
    return base, g3part


popc = np.array([bin(s).count("1") for s in range(dim)])
F0, F1, F2 = [], [], []
NT = 24
for tr in range(NT):
    r = np.random.default_rng(500 + tr)
    A_ = r.normal(size=(M, M)) + 1j * r.normal(size=(M, M))
    h = 0.25 * (A_ + A_.conj().T)
    Bt = r.normal(size=(M,) * 4) + 1j * r.normal(size=(M,) * 4)
    V = Bt - np.transpose(Bt, (1, 0, 2, 3))
    V = V - np.transpose(V, (0, 1, 3, 2))
    V = 0.25 * (V + np.transpose(V.conj(), (2, 3, 0, 1)))
    psi = np.where(popc == N, r.normal(size=dim) + 1j * r.normal(size=dim), 0.0)
    psi = psi / np.linalg.norm(psi)
    g2 = np.zeros((M,) * 4, dtype=complex)
    for t in itertools.product(*([range(M)] * 4)):
        g2[t] = np.vdot(psi, ad(t[0], ad(t[1], aa(t[3], aa(t[2], psi)))))
    rho = 0.5 * (to_pair(g2) + to_pair(g2).conj().T)
    w, U = np.linalg.eigh(rho)
    for gi in (1, 2):
        occ, gv = w[::-1][gi], U[:, ::-1][:, gi]
        rho_c = rho + (CEIL - occ) * np.outer(gv, gv.conj())
        g2c = from_pair(0.5 * (rho_c + rho_c.conj().T))
        g1c = contract_to_1rdm(g2c, N)
        b0, _ = pieces(h, V, g1c, g2c, np.zeros((M,) * 6, dtype=complex), gv)
        _, cA = pieces(h, V, g1c, g2c, orbit("g1g1g1", g1c, g2c), gv)
        _, cB = pieces(h, V, g1c, g2c, orbit("g1l2", g1c, g2c), gv)
        F0.append(b0); F1.append(cA); F2.append(cB)

F0, F1, F2 = np.array(F0), np.array(F1), np.array(F2)
print("samples: %d  (2 geminals x %d Hamiltonians), all deformed to n_g = %.3f" % (len(F0), NT, CEIL))
print()
print("gamma_3-free part F0 at the ceiling : max |F0| = %.3e" % np.abs(F0).max())
print("   -> the gamma_2-only terms contribute NOTHING to geminal occupation change.")
print("      gamma_3 is the entire driver.")
print()
print("orbit contributions at the ceiling : |F1| mean %.4f    |F2| mean %.4f"
      % (np.abs(F1).mean(), np.abs(F2).mean()))
print("lambda_3=0 value (a=b=1)           : mean %.4f, of which %d/%d are OUTWARD (>0)"
      % (np.abs(F1 + F2).mean(), int((F1 + F2 > 0).sum()), len(F0)))
print()

ang = np.sort(np.arctan2(F2, F1))
gaps = np.diff(np.concatenate([ang, [ang[0] + 2 * np.pi]]))
gmax = gaps.max()
print("=== half-plane feasibility: does some (a,b) give a*F1 + b*F2 <= 0 for ALL H? ===")
print("   largest angular gap between direction vectors : %.4f rad  (need > pi = %.4f)"
      % (gmax, np.pi))
if gmax > np.pi + 1e-9:
    k = int(np.argmax(gaps))
    mid = ang[k] + gmax / 2
    a, b = -np.cos(mid), -np.sin(mid)
    vals = a * F1 + b * F2
    print("   FEASIBLE. (a,b) = (%.4f, %.4f) gives max ndot(ceiling) = %.4e"
          % (a, b, vals.max()))
else:
    print("   INFEASIBLE. The directions positively span the plane, so for any (a,b)")
    print("   there is a Hamiltonian driving the geminal OUT through its ceiling.")
    for a, b in ((1, 1), (1, 0), (0, 1), (-1, -1)):
        v = a * F1 + b * F2
        print("      (a,b)=(%2d,%2d): %d of %d samples outward, worst %+0.4f"
              % (a, b, int((v > 0).sum()), len(v), v.max()))
