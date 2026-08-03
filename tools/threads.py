#!/usr/bin/env python3
"""Three threads on the N=3 result.

1. Does gamma_3 = 0 preserve D across many random (psi0, H)? If it never fails,
   that is evidence for a structural theorem rather than a lucky trajectory.
2. Nagumo, evaluated for real: at t=0 the exact 2-RDM sits ON the Q boundary
   (min eig ~ 0), so for every null vector v the flux <v|Qdot|v> must be >= 0
   under exact dynamics. Compute it for exact and for the closure.
3. Contraction consistency: does the closure drive the natural occupations of
   gamma_1 = contract(gamma_2)/(N-1) above 1? Q carries hole factors, so
   occupations above 1 would break it -- the same failure mode as docs/01's
   2.81 electrons, reappearing at rank 2.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from qdiss.symbolic import mkno
from qdiss.dynamics.coherent2 import compile_terms, evaluate, rk4
from qdiss.representability.dqg import contract_to_1rdm, build_Q, build_G
from tools.eom2_coherent import ONE, TWO

M, N = 6, 3
dim = 1 << M
ONE_C = compile_terms(ONE, "rs")
TWO_C = compile_terms(TWO, "rtsu")

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


def rdm2(psi):
    g = np.zeros((M,) * 4, dtype=complex)
    for p in range(M):
        for v in range(M):
            for q in range(M):
                for w in range(M):
                    g[p, v, q, w] = np.vdot(psi, ad(p, ad(v, aa(w, aa(q, psi)))))
    return g


def sep_gamma3(g1, l2):
    t = {("G", 1): g1, ("L", 2): l2}
    out = np.zeros((M,) * 6, dtype=complex)
    for sign, sub, keys in SEP3:
        out = out + sign * np.einsum(sub, *[t[k] for k in keys], optimize=True)
    return out


def lam2_of(g2, g1):
    return g2 - (np.einsum("pq,vw->pvqw", g1, g1) - np.einsum("pw,vq->pvqw", g1, g1))


def make(h, V, alpha):
    def rhs(g2):
        g1 = contract_to_1rdm(g2, N)
        g3 = (alpha * sep_gamma3(g1, lam2_of(g2, g1)) if alpha
              else np.zeros((M,) * 6, dtype=complex))
        gam = {2: g2, 3: g3}
        return 1j * (evaluate(ONE_C, h, gam, M, out_rank=4)
                     + 0.25 * evaluate(TWO_C, V, gam, M, out_rank=4))
    return rhs


def mat(T):
    A = T.reshape(M * M, M * M)
    return 0.5 * (A + A.conj().T)


def random_system(rng):
    A = rng.normal(size=(M, M)) + 1j * rng.normal(size=(M, M))
    h = 0.5 * (A + A.conj().T)
    B = rng.normal(size=(M,) * 4) + 1j * rng.normal(size=(M,) * 4)
    V = B - np.transpose(B, (1, 0, 2, 3))
    V = V - np.transpose(V, (0, 1, 3, 2))
    V = 0.5 * (V + np.transpose(V.conj(), (2, 3, 0, 1)))
    psi = np.where(np.array([bin(s).count("1") for s in range(dim)]) == N,
                   rng.normal(size=dim) + 1j * rng.normal(size=dim), 0.0)
    return 0.5 * h, 0.5 * V, psi / np.linalg.norm(psi)


# ================= thread 1 + 3 : scan =================
print("=== threads 1 & 3: scan over random (psi0, H), closure gamma_3 = 0 ===")
print("%4s | %11s %11s | %9s %9s" % ("run", "min D", "min Q", "max occ", "tr g1"))
dt, T = 0.004, 0.8
nstep = int(T / dt)
worstD = worstQ = 0.0
maxocc_all = 0.0
for run in range(12):
    rng = np.random.default_rng(1000 + run)
    h, V, psi = random_system(rng)
    g2 = rdm2(psi)
    rhs = make(h, V, 0.0)
    mD = mQ = 1e9
    mo = 0.0
    tr = 0.0
    for step in range(nstep + 1):
        g1 = contract_to_1rdm(g2, N)
        mD = min(mD, float(np.linalg.eigvalsh(mat(g2)).min()))
        mQ = min(mQ, float(np.linalg.eigvalsh(mat(build_Q(g2, g1))).min()))
        occ = np.linalg.eigvalsh(0.5 * (g1 + g1.conj().T))
        mo = max(mo, float(occ.max()))
        tr = float(np.trace(g1).real)
        g2 = rk4(rhs, g2, dt)
    worstD = min(worstD, mD)
    worstQ = min(worstQ, mQ)
    maxocc_all = max(maxocc_all, mo)
    print("%4d | %11.3e %11.3e | %9.6f %9.6f" % (run, mD, mQ, mo, tr))
print()
print("  worst min D over all runs : %.3e" % worstD)
print("  worst min Q over all runs : %.3e" % worstQ)
print("  largest occupation seen   : %.6f  (must be <= 1 for Q)" % maxocc_all)
print()

# ================= thread 2 : Nagumo flux at t=0 =================
print("=== thread 2: boundary flux <v|Mdot|v> at t=0 on the exact 2-RDM ===")
rng = np.random.default_rng(2026)
h, V, psi = random_system(rng)
g2 = rdm2(psi)
g1 = contract_to_1rdm(g2, N)


def flux(alpha):
    rhs = make(h, V, alpha)
    dg2 = rhs(g2)
    dg1 = contract_to_1rdm(dg2, N)
    dQ = (- np.einsum("ki,lj->ijkl", np.eye(M), dg1)
          + np.einsum("kj,li->ijkl", np.eye(M), dg1)
          + np.einsum("li,kj->ijkl", np.eye(M), dg1)
          - np.einsum("lj,ki->ijkl", np.eye(M), dg1)
          + np.einsum("klij->ijkl", dg2))
    dG = (np.einsum("jl,ik->ijkl", np.eye(M), dg1) - dg2.transpose(0, 3, 2, 1))
    return dg2, dQ, dG


for name, T_ in (("D", g2), ("Q", build_Q(g2, g1)), ("G", build_G(g2, g1))):
    w, U = np.linalg.eigh(mat(T_))
    null = U[:, w < 1e-9]
    print("  %s : %d null directions of %d" % (name, null.shape[1], len(w)))
    for alpha, lbl in ((0.0, "gamma_3=0"), (1.0, "lambda_3=0")):
        dD, dQ, dG = flux(alpha)
        dot = {"D": dD, "Q": dQ, "G": dG}[name]
        vals = np.real(np.einsum("ai,ab,bi->i", null.conj(), mat(dot), null))
        print("      %-11s min <v|%sdot|v> = %+.6e   %s"
              % (lbl, name, vals.min(), "INWARD" if vals.min() > -1e-9 else "OUTWARD"))
