#!/usr/bin/env python3
"""Rank-1 Nagumo, on exactly-sampled edges of 0 <= gamma <= 1.

The rank-1 representable set is the matrix interval 0 <= gamma <= 1, which has TWO
boundaries. Nagumo there is:

    n = 0 face, null vector v of gamma      :  <v|gamma-dot|v>  >= 0
    n = 1 face, null vector v of (1-gamma)  :  <v|gamma-dot|v>  <= 0

docs/01 section 4.2 says the exact rank-1 mechanism makes BOTH vanish identically,
because gain terms carry gamma and loss terms carry (1-gamma) -- the particle side
and the hole side. So this is the rung where a two-sided CP structure is known to
exist, and the question is which side, if either, a generated lambda_2 respects.

Edges are sampled exactly, not approached: project gamma so a chosen direction has
eigenvalue exactly 0 (particle face) or exactly 1 (hole face), then evaluate the
flux there. If the flux is inward at 0 and outward at 1, the hole side is the
failure mode -- which would match hole-CP being unreachable at rank 2.
"""
import itertools
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from qdiss.symbolic import mkno
from qdiss.dynamics.coherent2 import compile_terms, evaluate
from tools.eom1 import one as ONE1, two as TWO1
from tools.eom2_coherent import ONE as ONE2, TWO as TWO2

M, N = 6, 3
dim = 1 << M
ONE1_C = compile_terms(ONE1, "rs", out="pq")
TWO1_C = compile_terms(TWO1, "rtsu", out="pq")
ONE2_C = compile_terms(ONE2, "rs")
TWO2_C = compile_terms(TWO2, "rtsu")
UP3, LO3 = ("p", "v", "x"), ("q", "w", "y")
sep3 = [(s, f) for s, f in mkno.gamma_to_cumulants(UP3, LO3)
        if not any(n == "L" and r == 3 for n, r, u, l in f)]
SEP3 = [(s, ",".join("".join(u) + "".join(l) for n, r, u, l in f) + "->pvxqwy",
         [(n, r) for n, r, u, l in f]) for s, f in sep3]

rng = np.random.default_rng(2026)
A_ = rng.normal(size=(M, M)) + 1j * rng.normal(size=(M, M))
h = 0.25 * (A_ + A_.conj().T)
B = rng.normal(size=(M,) * 4) + 1j * rng.normal(size=(M,) * 4)
V = B - np.transpose(B, (1, 0, 2, 3))
V = V - np.transpose(V, (0, 1, 3, 2))
V = 0.25 * (V + np.transpose(V.conj(), (2, 3, 0, 1)))


def wedge2(g1):
    return np.einsum("pq,vw->pvqw", g1, g1) - np.einsum("pw,vq->pvqw", g1, g1)


def sep_g3(g1, l2):
    t = {("G", 1): g1, ("L", 2): l2}
    out = np.zeros((M,) * 6, dtype=complex)
    for sign, sub, keys in SEP3:
        out = out + sign * np.einsum(sub, *[t[k] for k in keys], optimize=True)
    return out


def g1dot(g1, g2):
    gam = {1: g1, 2: g2}
    return 1j * (evaluate(ONE1_C, h, gam, M, out_rank=2)
                 + 0.25 * evaluate(TWO1_C, V, gam, M, out_rank=2))


def herm(A):
    return 0.5 * (A + A.conj().T)


def make_edge(g, v, face):
    """Force eigenvalue exactly 0 (face=0) or exactly 1 (face=1) along v."""
    P = np.eye(M, dtype=complex) - np.outer(v, v.conj())
    if face == 0:
        return herm(P @ g @ P)
    return herm(P @ g @ P) + np.outer(v, v.conj())


def occ_ok(g):
    w = np.linalg.eigvalsh(herm(g))
    return w.min() > -1e-9 and w.max() < 1 + 1e-9, w


# a correlated gamma_1 to project from: take one from a short exact propagation
def ad(p, vv):
    out = np.zeros_like(vv); bit = 1 << p
    for s in np.nonzero(vv)[0]:
        if s & bit: continue
        out[s | bit] += ((-1.0) ** bin(s & (bit - 1)).count("1")) * vv[s]
    return out


def aa(p, vv):
    out = np.zeros_like(vv); bit = 1 << p
    for s in np.nonzero(vv)[0]:
        if not (s & bit): continue
        out[s ^ bit] += ((-1.0) ** bin(s & (bit - 1)).count("1")) * vv[s]
    return out


def apply_str(up, lo, vv):
    for q in lo:
        vv = aa(q, vv)
    for p in reversed(up):
        vv = ad(p, vv)
    return vv


def rdm(psi, k):
    g = np.zeros((M,) * (2 * k), dtype=complex)
    for t in itertools.product(*([range(M)] * (2 * k))):
        g[t] = np.vdot(psi, apply_str(t[:k], t[k:], psi))
    return g


def Hpsi(vv):
    out = np.zeros_like(vv)
    for r in range(M):
        for s in range(M):
            out += h[r, s] * ad(r, aa(s, vv))
    for r in range(M):
        for t in range(M):
            for s in range(M):
                for u in range(M):
                    c = V[r, t, s, u]
                    if c != 0:
                        out += 0.25 * c * ad(r, ad(t, aa(u, aa(s, vv))))
    return out


psi = np.zeros(dim, dtype=complex)
psi[sum(1 << i for i in range(N))] = 1.0
dt = 0.0025
for _ in range(40):
    k1 = -1j * Hpsi(psi); k2 = -1j * Hpsi(psi + 0.5 * dt * k1)
    k3 = -1j * Hpsi(psi + 0.5 * dt * k2); k4 = -1j * Hpsi(psi + dt * k3)
    psi = psi + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
g1_ref = herm(rdm(psi, 1))
l2_ref = rdm(psi, 2) - wedge2(g1_ref)
print("reference gamma_1 occupations: %s"
      % np.array2string(np.linalg.eigvalsh(g1_ref), precision=4))
print()

print("flux at exactly-sampled edges; inward means >=0 at n=0 and <=0 at n=1")
print()
print("%5s | %-30s | %-30s" % ("", "n = 0 face (particle)", "n = 1 face (hole)"))
print("%5s | %14s %14s | %14s %14s"
      % ("trial", "TDHF", "generated", "TDHF", "generated"))
bad = {"t0": 0, "g0": 0, "t1": 0, "g1": 0}
ntr = 16
for t_ in range(ntr):
    r = np.random.default_rng(6100 + t_)
    v = r.normal(size=M) + 1j * r.normal(size=M)
    v /= np.linalg.norm(v)
    row = ["%5d |" % t_]
    for face in (0, 1):
        ge = make_edge(g1_ref, v, face)
        ok, w = occ_ok(ge)
        for tag, l2 in (("t", np.zeros((M,) * 4, dtype=complex)), ("g", l2_ref)):
            g2 = wedge2(ge) + l2
            d1 = herm(g1dot(ge, g2))
            f = float(np.real(v.conj() @ d1 @ v))
            key = tag + str(face)
            if (face == 0 and f < -1e-9) or (face == 1 and f > 1e-9):
                bad[key] += 1
            row.append(" %14.4e" % f)
        row.append(" |")
    print("".join(row))

print()
print("escapes:  n=0 face  TDHF %d/%d  generated %d/%d"
      % (bad["t0"], ntr, bad["g0"], ntr))
print("          n=1 face  TDHF %d/%d  generated %d/%d"
      % (bad["t1"], ntr, bad["g1"], ntr))
