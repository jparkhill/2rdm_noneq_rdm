#!/usr/bin/env python3
"""Rank-1 program: keep gamma, generate lambda_2, watch occupations and positivity.

TDHF is i gamma-dot = [f(gamma), gamma], a commutator, so it is ISOSPECTRAL --
occupations are frozen, exactly as gamma_3 = 0 is frozen at rank 2. Correcting the
one-particle dynamics means moving occupations, and the question is whether that
can be done while keeping 0 <= n <= 1.

Three dynamics compared against exact propagation:

  TDHF        lambda_2 = 0
  generated   lambda_2 integrated along the trajectory from its own EOM with
              lambda_3 = 0, so gamma_3 = separable(gamma_1, lambda_2)
  exact       full wavefunction

Reports max |n - n_exact|, and the occupation window, since leaving [0,1] is the
rank-1 form of leaving the cone.
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


def g2dot(g1, g2, g3):
    gam = {1: g1, 2: g2, 3: g3}
    return 1j * (evaluate(ONE2_C, h, gam, M, out_rank=4)
                 + 0.25 * evaluate(TWO2_C, V, gam, M, out_rank=4))


def rhs_tdhf(y):
    g1 = y[0]
    return (g1dot(g1, wedge2(g1)), np.zeros((M,) * 4, dtype=complex))


def rhs_gen(y):
    g1, l2 = y
    g2 = wedge2(g1) + l2
    d1 = g1dot(g1, g2)
    d2 = g2dot(g1, g2, sep_g3(g1, l2))
    dl2 = d2 - (np.einsum("pq,vw->pvqw", d1, g1) - np.einsum("pw,vq->pvqw", d1, g1)
                + np.einsum("pq,vw->pvqw", g1, d1) - np.einsum("pw,vq->pvqw", g1, d1))
    return (d1, dl2)


def axpy(y, a, k):
    return (y[0] + a * k[0], y[1] + a * k[1])


def step(rhs, y, dt):
    k1 = rhs(y)
    k2 = rhs(axpy(y, 0.5 * dt, k1))
    k3 = rhs(axpy(y, 0.5 * dt, k2))
    k4 = rhs(axpy(y, dt, k3))
    return (y[0] + (dt / 6.0) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]),
            y[1] + (dt / 6.0) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]))


psi0 = np.zeros(dim, dtype=complex)
psi0[sum(1 << i for i in range(N))] = 1.0
g1_0 = rdm(psi0, 1)

yT = (g1_0.copy(), np.zeros((M,) * 4, dtype=complex))
yG = (g1_0.copy(), np.zeros((M,) * 4, dtype=complex))
psi = psi0.copy()

dt = 0.0025
print("start: determinant, M=%d N=%d.  occupations start at exactly 0 or 1." % (M, N))
print()
print("%7s | %-26s | %-26s | %s" % ("t", "TDHF", "generated lambda_2", "exact"))
print("%7s | %9s %7s %7s | %9s %7s %7s | %7s %7s"
      % ("", "err(n)", "n_min", "n_max", "err(n)", "n_min", "n_max", "n_min", "n_max"))
for stp in range(81):
    if stp % 10 == 0:
        g1e = rdm(psi, 1)
        ne = np.sort(np.linalg.eigvalsh(0.5 * (g1e + g1e.conj().T)))
        row = ["%7.4f |" % (stp * dt)]
        for y in (yT, yG):
            g = y[0]
            nn = np.sort(np.linalg.eigvalsh(0.5 * (g + g.conj().T)))
            row.append(" %9.2e %7.4f %7.4f |"
                       % (np.abs(nn - ne).max(), nn.min(), nn.max()))
        row.append(" %7.4f %7.4f" % (ne.min(), ne.max()))
        print("".join(row))
    yT = step(rhs_tdhf, yT, dt)
    yG = step(rhs_gen, yG, dt)
    psi = rk4(lambda v: -1j * Hpsi(v), psi, dt)
