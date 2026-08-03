#!/usr/bin/env python3
"""Is the gamma_3 = 0 closure isospectral in gamma_2?

Conjecture: dropping the gamma_3 terms leaves gamma_2-dot built only from the
four gamma_2 terms, which carry commutator structure -- an effective pair
Hamiltonian acting on the pair space. If gamma_2-dot = -i[K, gamma_2] then

  * <v|gamma_2-dot|v> = 0 identically for v in ker gamma_2, so D is preserved
    exactly rather than approximately -- explaining -3.8e-16 over 30 null
    directions, and
  * the eigenvalues of gamma_2 are constants of motion.

The second is the decisive test. It would also be a real limitation: an
isospectral pair flow is a unitary rotation of gamma_2 and cannot relax the pair
occupation distribution at all.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from qdiss.symbolic import mkno
from qdiss.dynamics.coherent2 import compile_terms, evaluate, rk4
from qdiss.representability.dqg import contract_to_1rdm
from tools.eom2_coherent import ONE, TWO

M, N = 6, 3
dim = 1 << M
ONE_C = compile_terms(ONE, "rs")
TWO_C = compile_terms(TWO, "rtsu")


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


rng = np.random.default_rng(2026)
A = rng.normal(size=(M, M)) + 1j * rng.normal(size=(M, M))
h = 0.25 * (A + A.conj().T)
B = rng.normal(size=(M,) * 4) + 1j * rng.normal(size=(M,) * 4)
V = B - np.transpose(B, (1, 0, 2, 3))
V = V - np.transpose(V, (0, 1, 3, 2))
V = 0.25 * (V + np.transpose(V.conj(), (2, 3, 0, 1)))
psi = np.where(np.array([bin(s).count("1") for s in range(dim)]) == N,
               rng.normal(size=dim) + 1j * rng.normal(size=dim), 0.0)
psi = psi / np.linalg.norm(psi)
g2 = rdm2(psi)


def rhs0(g2):
    gam = {2: g2, 3: np.zeros((M,) * 6, dtype=complex)}
    return 1j * (evaluate(ONE_C, h, gam, M, out_rank=4)
                 + 0.25 * evaluate(TWO_C, V, gam, M, out_rank=4))


def spec(g):
    A = g.reshape(M * M, M * M)
    return np.linalg.eigvalsh(0.5 * (A + A.conj().T))


w0 = spec(g2)
print("initial gamma_2 spectrum (nonzero part):")
print("  ", np.array2string(w0[w0 > 1e-10], precision=6))
print()

# --- test 1: is gamma_2-dot anti-Hermitian-commutator-like, i.e. tr(gamma_2^k) conserved? ---
dt = 0.002
nstep = 250
print("%6s | %12s %12s %12s | %12s" %
      ("t", "|dspec|max", "tr g2", "tr g2^2", "min eig"))
g = g2.copy()
for step in range(nstep + 1):
    if step % 50 == 0:
        w = spec(g)
        A = g.reshape(M * M, M * M)
        print("%6.3f | %12.3e %12.6f %12.6f | %12.3e"
              % (step * dt, np.abs(np.sort(w) - np.sort(w0)).max(),
                 np.trace(A).real, np.trace(A @ A).real, w.min()))
    g = rk4(rhs0, g, dt)

print()
# --- test 2: extract K and check gamma_2-dot = -i[K,gamma_2] directly ---
dg = rhs0(g2)
Ag = g2.reshape(M * M, M * M)
Adg = dg.reshape(M * M, M * M)
print("is gamma_2-dot anti-Hermitian in pair space? |dot + dot^H| = %.3e"
      % np.abs(Adg + Adg.conj().T).max())
print("  (a commutator -i[K,X] with K,X Hermitian IS anti-Hermitian)")
comm_resid = np.abs(Adg @ Ag - Ag @ Adg)
print("does it commute with gamma_2? (it need not) |[dot,g2]| = %.3e" % comm_resid.max())
