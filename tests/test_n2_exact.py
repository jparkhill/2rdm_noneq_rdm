#!/usr/bin/env python3
"""N=2 positive control: does the coherent 2-RDM EOM reproduce exact dynamics?

At N = 2, gamma_3 vanishes identically, so closing on gamma_3 = 0 is not an
approximation -- the EOM is the exact dynamics. Two things must then hold:

  1. gamma_2(t) from the EOM equals gamma_2(t) from propagating the wavefunction
     in the full Fock space, to integrator tolerance.
  2. the D condition (gamma_2 PSD) holds for all t, automatically.

Random Hermitian h and a random antisymmetrized V, which is a harder test than
any particular molecule.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from qdiss.dynamics.coherent2 import Coherent2, closure_gamma3_zero, rk4
from tools.eom2_coherent import ONE, TWO

M, N = 4, 2
rng = np.random.default_rng(5)
dim = 1 << M


def occ_bits(s):
    return [i for i in range(M) if (s >> i) & 1]


def ad(p, v):
    out = np.zeros_like(v)
    bit = 1 << p
    for s in np.nonzero(v)[0]:
        if s & bit:
            continue
        sgn = (-1.0) ** bin(s & (bit - 1)).count("1")
        out[s | bit] += sgn * v[s]
    return out


def aa(p, v):
    out = np.zeros_like(v)
    bit = 1 << p
    for s in np.nonzero(v)[0]:
        if not (s & bit):
            continue
        sgn = (-1.0) ** bin(s & (bit - 1)).count("1")
        out[s ^ bit] += sgn * v[s]
    return out


# --- Hamiltonian: h hermitian, V with the antisymmetry of a^{rt}_{su} ---
A = rng.normal(size=(M, M)) + 1j * rng.normal(size=(M, M))
h = (A + A.conj().T) / 2

B = rng.normal(size=(M,) * 4) + 1j * rng.normal(size=(M,) * 4)
V = B - np.transpose(B, (1, 0, 2, 3))
V = V - np.transpose(V, (0, 1, 3, 2))
V = (V + np.transpose(V.conj(), (2, 3, 0, 1))) / 2


def H_apply(v):
    out = np.zeros_like(v)
    for r in range(M):
        for s in range(M):
            if h[r, s] != 0:
                out += h[r, s] * ad(r, aa(s, v))
    for r in range(M):
        for t in range(M):
            for s in range(M):
                for u in range(M):
                    c = V[r, t, s, u]
                    if c != 0:
                        out += 0.5 * c * ad(r, ad(t, aa(u, aa(s, v))))
    return out


def gamma2_of(psi):
    g = np.zeros((M,) * 4, dtype=complex)
    for p in range(M):
        for v_ in range(M):
            ket = ad(p, ad(v_, np.zeros_like(psi)))  # placeholder, unused
    for p in range(M):
        for v_ in range(M):
            for q in range(M):
                for w in range(M):
                    ket = ad(p, ad(v_, aa(w, aa(q, psi))))
                    g[p, v_, q, w] = np.vdot(psi, ket)
    return g


# --- initial state: a random correlated 2-electron state ---
popc = np.array([bin(s).count("1") for s in range(dim)])
psi = np.where(popc == N, rng.normal(size=dim) + 1j * rng.normal(size=dim), 0.0)
psi = psi / np.linalg.norm(psi)

g2 = gamma2_of(psi)
eom = Coherent2(ONE, TWO, h, V, closure_gamma3_zero(M))

import os
dt = float(os.environ.get("DT",0.002))
nstep = int(round(0.5/dt))
print("M=%d N=%d  dt=%g  steps=%d" % (M, N, dt, nstep))
print()
print("%6s  %14s  %14s  %14s" % ("t", "|g2_eom-g2_exact|", "min eig D", "tr g2"))

worst = 0.0
mineig = 1e9
for step in range(nstep + 1):
    if step % 50 == 0:
        gx = gamma2_of(psi)
        err = np.abs(g2 - gx).max()
        worst = max(worst, err)
        D = g2.reshape(M * M, M * M)
        D = (D + D.conj().T) / 2
        ev = np.linalg.eigvalsh(D).min()
        mineig = min(mineig, ev)
        tr = np.einsum("pqpq->", g2).real
        print("%6.3f  %14.3e  %14.6f  %14.6f" % (step * dt, err, ev, tr))
    # advance both
    g2 = rk4(eom.ddt, g2, dt)
    psi = rk4(lambda v: -1j * H_apply(v), psi, dt)

print()
print("max |EOM - exact| over trajectory : %.3e" % worst)
print("min eigenvalue of D seen          : %.3e" % mineig)
print("VERDICT:", "EOM IS EXACT AT N=2 AND STAYS PSD"
      if worst < 1e-6 and mineig > -1e-6 else "PROBLEM")
