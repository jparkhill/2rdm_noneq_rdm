#!/usr/bin/env python3
"""Does a DYNAMICALLY GENERATED lambda_2 track the true one?

The same construction one rung down, where the answer is known. From a
determinant, lambda_2 = lambda_3 = 0 exactly, so

    lambda_2-dot = gamma_2-dot - gamma_1-dot ^ gamma_1 - gamma_1 ^ gamma_1-dot

is computable from gamma_1 and the Hamiltonian alone: gamma_2 = gamma_1 ^ gamma_1
and gamma_3 is the Wick reconstruction, both exact at a determinant. So the
leading generated estimate is

    lambda_2(t)  ~=  t * lambda_2-dot(0)

Compare that against the true lambda_2(t) from exact propagation. If it tracks,
the construction works at rank 1, where the closure is cheap (propagate the
2-index gamma_1 and reconstruct the 4-index lambda_2) and where the positivity
condition is only 0 <= n <= 1 rather than DQG.

Baseline to beat: lambda_2 = 0, which is TDHF.
"""
import itertools
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from qdiss.dynamics.coherent2 import compile_terms, evaluate
from qdiss.representability.dqg import contract_to_1rdm
from tools.eom1 import one as ONE1, two as TWO1
from tools.eom2_coherent import ONE as ONE2, TWO as TWO2

M, N = 6, 3
dim = 1 << M
ONE1_C = compile_terms(ONE1, "rs", out="pq")
TWO1_C = compile_terms(TWO1, "rtsu", out="pq")
ONE2_C = compile_terms(ONE2, "rs")
TWO2_C = compile_terms(TWO2, "rtsu")


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


def g1dot(g1, g2):
    gam = {1: g1, 2: g2}
    return 1j * (evaluate(ONE1_C, h, gam, M, out_rank=2)
                 + 0.25 * evaluate(TWO1_C, V, gam, M, out_rank=2))


def g2dot(g1, g2, g3):
    gam = {1: g1, 2: g2, 3: g3}
    return 1j * (evaluate(ONE2_C, h, gam, M, out_rank=4)
                 + 0.25 * evaluate(TWO2_C, V, gam, M, out_rank=4))


# ---- determinant reference ----
psi0 = np.zeros(dim, dtype=complex)
psi0[sum(1 << i for i in range(N))] = 1.0
g1_0 = rdm(psi0, 1)
g2_0 = rdm(psi0, 2)
g3_0 = rdm(psi0, 3)
print("determinant check: |gamma_2 - gamma_1^gamma_1| = %.3e"
      % np.abs(g2_0 - wedge2(g1_0)).max())

d1 = g1dot(g1_0, g2_0)
d2 = g2dot(g1_0, g2_0, g3_0)
l2dot = d2 - (np.einsum("pq,vw->pvqw", d1, g1_0)
              - np.einsum("pw,vq->pvqw", d1, g1_0)
              + np.einsum("pq,vw->pvqw", g1_0, d1)
              - np.einsum("pw,vq->pvqw", g1_0, d1))
print("|lambda_2-dot(0)| = %.6e" % np.abs(l2dot).max())
print()

print("%8s %14s %14s %12s %12s" % ("t", "|l2 true|", "|l2 gen|", "rel err", "vs TDHF"))
dt = 0.0025
psi = psi0.copy()
for step in range(41):
    t = step * dt
    if step % 5 == 0 and step > 0:
        g2 = rdm(psi, 2)
        g1 = contract_to_1rdm(g2, N)
        l2_true = g2 - wedge2(g1)
        l2_gen = t * l2dot
        nt = np.abs(l2_true).max()
        ng = np.abs(l2_gen).max()
        err = np.abs(l2_gen - l2_true).max() / nt
        tdhf = 1.0            # lambda_2 = 0 has relative error 1 by definition
        print("%8.4f %14.6e %14.6e %12.4f %12s"
              % (t, nt, ng, err, "%.2fx better" % (tdhf / err) if err > 0 else "-"))
    psi = rk4(lambda v: -1j * Hpsi(v), psi, dt)
