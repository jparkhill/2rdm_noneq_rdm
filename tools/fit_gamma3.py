#!/usr/bin/env python3
"""Can a reconstruction gamma_3[gamma_1, gamma_2] carry Coleman's ceiling?

The separable expansion splits into two antisymmetry-preserving orbits,

    gamma_3 = a * (gamma_1 wedge gamma_1 wedge gamma_1) + b * (gamma_1 wedge lambda_2)

with a = b = 1 the lambda_3 = 0 reconstruction. Coefficients may differ between
orbits without breaking the antisymmetry of gamma_3, but not within an orbit, so
these two numbers are the whole freedom at fixed (M, N).

Two independent demands on (a, b):

  1. Exactness at N = 2. A two-electron state has gamma_3 = 0 identically, so the
     reconstruction must vanish there. lambda_3 = 0 does not -- that was measured at
     3.16e-2 earlier. Is there a ratio a/b that does?

  2. The ceiling. gamma_3 enters gamma_2-dot linearly, so

         ndot(c) = base(H) + a * A(H) + b * B(H)

     is affine in (a, b). Demanding it vanish at s = c for EVERY Hamiltonian is a
     least-squares problem over sampled H. A small residual means the family can host
     the blocking factor; a residual comparable to the base scale means it cannot.

Then: are the two demands compatible?
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
ORB = {"g1g1g1": [], "g1l2": []}
for sign, facs in mkno.gamma_to_cumulants(UP3, LO3):
    if any(n == "L" and r == 3 for n, r, u, l in facs):
        continue
    key = "g1l2" if any(n == "L" and r == 2 for n, r, u, l in facs) else "g1g1g1"
    ORB[key].append((sign, ",".join("".join(u) + "".join(l) for n, r, u, l in facs)
                     + "->pvxqwy", [(n, r) for n, r, u, l in facs]))
print("orbit sizes: gamma_1^3 = %d terms, gamma_1*lambda_2 = %d terms"
      % (len(ORB["g1g1g1"]), len(ORB["g1l2"])))


def build(M):
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
    return ad, aa


def orbit(key, g1, g2, M):
    l2 = g2 - (np.einsum("pq,vw->pvqw", g1, g1) - np.einsum("pw,vq->pvqw", g1, g1))
    t = {("G", 1): g1, ("L", 2): l2}
    out = np.zeros((M,) * 6, dtype=complex)
    for sign, sub, keys in ORB[key]:
        out = out + sign * np.einsum(sub, *[t[k] for k in keys], optimize=True)
    return out


print()
print("=== 1. exactness at N = 2, where gamma_3 vanishes identically ===")
M2, N2 = 6, 2
ad, aa = build(M2)
dim2 = 1 << M2
popc2 = np.array([bin(s).count("1") for s in range(dim2)])
rows = []
for tr in range(6):
    r = np.random.default_rng(70 + tr)
    psi = np.where(popc2 == N2, r.normal(size=dim2) + 1j * r.normal(size=dim2), 0.0)
    psi = psi / np.linalg.norm(psi)
    g1 = np.array([[np.vdot(psi, ad(p, aa(q, psi))) for q in range(M2)] for p in range(M2)])
    g2 = np.zeros((M2,) * 4, dtype=complex)
    for t in itertools.product(*([range(M2)] * 4)):
        g2[t] = np.vdot(psi, ad(t[0], ad(t[1], aa(t[3], aa(t[2], psi)))))
    g3 = np.zeros((M2,) * 6, dtype=complex)
    for t in itertools.product(*([range(M2)] * 6)):
        g3[t] = np.vdot(psi, ad(t[0], ad(t[1], ad(t[2], aa(t[5], aa(t[4], aa(t[3], psi)))))))
    OA, OB = orbit("g1g1g1", g1, g2, M2), orbit("g1l2", g1, g2, M2)
    rows.append((np.abs(g3).max(), OA.ravel(), OB.ravel(), np.abs(OA + OB).max()))
    print("   trial %d: |exact gamma_3| = %.2e   |a=b=1 reconstruction| = %.3e"
          % (tr, rows[-1][0], rows[-1][3]))
Amat = np.concatenate([np.stack([r[1], r[2]], axis=1) for r in rows])
Amat = np.concatenate([Amat.real, Amat.imag])
u, sv, vt = np.linalg.svd(Amat, full_matrices=False)
print("   singular values of [orbitA | orbitB] stacked over trials: %s" % np.array2string(sv, precision=4))
if sv[-1] / sv[0] < 1e-8:
    ab = vt[-1]
    print("   NULL DIRECTION a/b = %.6f  -> a reconstruction vanishing at N=2 exists" % (ab[0] / ab[1]))
else:
    print("   ratio sv_min/sv_max = %.3e : orbits are linearly INDEPENDENT," % (sv[-1] / sv[0]))
    print("   so no (a,b) other than (0,0) makes gamma_3 vanish at N = 2.")

print()
print("=== 2. can (a,b) put a root at the ceiling for every H? ===")
M, N = 6, 3
CEIL = N * (M - N + 2) / M
dim = 1 << M
PAIRS = [(p, v) for p in range(M) for v in range(p + 1, M)]
ad, aa = build(M)
popc = np.array([bin(s).count("1") for s in range(dim)])


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
    """returns (base, gamma_3 contribution) of <g|gamma_2-dot|g>"""
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


F0, F1, F2 = [], [], []
for tr in range(16):
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
    occ, gv = w[::-1][1], U[:, ::-1][:, 1]
    rho_c = rho + (CEIL - occ) * np.outer(gv, gv.conj())
    g2c = from_pair(0.5 * (rho_c + rho_c.conj().T))
    g1c = contract_to_1rdm(g2c, N)
    z3 = np.zeros((M,) * 6, dtype=complex)
    b0, _ = pieces(h, V, g1c, g2c, z3, gv)
    _, cA = pieces(h, V, g1c, g2c, orbit("g1g1g1", g1c, g2c, M), gv)
    _, cB = pieces(h, V, g1c, g2c, orbit("g1l2", g1c, g2c, M), gv)
    F0.append(b0); F1.append(cA); F2.append(cB)

F0, F1, F2 = np.array(F0), np.array(F1), np.array(F2)
Amat = np.stack([F1, F2], axis=1)
sol, *_ = np.linalg.lstsq(Amat, -F0, rcond=None)
res = Amat @ sol + F0
print("   scale of base term at the ceiling : %.4f" % np.abs(F0).mean())
print("   best (a, b) = (%.4f, %.4f)" % (sol[0], sol[1]))
print("   residual |ndot(ceiling)| mean %.4f  max %.4f   (lambda_3=0 gives mean %.4f)"
      % (np.abs(res).mean(), np.abs(res).max(), np.abs(F0 + F1 + F2).mean()))
print("   relative residual = %.3f" % (np.abs(res).mean() / np.abs(F0).mean()))
if np.abs(res).mean() / np.abs(F0).mean() < 0.05:
    print("   -> the family CAN host the ceiling root")
else:
    print("   -> the family CANNOT host the ceiling root; no (a,b) works across H")
