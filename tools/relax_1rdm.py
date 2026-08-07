#!/usr/bin/env python3
"""Secular Redfield relaxation on the 1-RDM EOM: does it reach a stationary state?

APPROXIMATION. Propagated variable gamma_1 alone. Coherent part is the MKNO 1-RDM EOM,
qdiss.eom.gwt.terms_cumulant(1), a function of gamma_1 and lambda_2 only, validated against
qdiss.eom.fullci at 1e-15. Two closures for the coherent part:
    lam2_0   two-particle cumulant TRUNCATED, lambda_2 := 0. This is TDHF, i[(h+W)^T, gamma_1],
             and is the control: its stationary state is known to be finite-temperature
             Hartree-Fock, going to the HF determinant as T -> 0.
    blocked  ours, lambda_2 reconstructed to satisfy the partial-trace identity exactly
             (qdiss.closures, A^f + f^A). No stationary state is known for it -- that is the
             open question this tool exists to probe.

DISSIPATOR, secular Redfield with detailed balance. In the instantaneous eigenbasis of the
generator A = (h + W)^T with eigenvalues eps_i:
    populations  ndot_i = sum_{j!=i} [ k_ij n_j (1-n_i) - k_ji n_i (1-n_j) ]
    coherences   gtilde_ij -> -Gamma * gtilde_ij  for i != j
    rates        k_ij = gamma_0 / (1 + exp(beta (eps_i - eps_j)))
which satisfies k_ij/k_ji = exp(-beta(eps_i - eps_j)), so the population fixed point is
Fermi-Dirac and the n_j(1-n_i) factors are the usual Pauli blocking. Flat spectral density and a
constant pure-dephasing rate are DEFAULTS chosen here, not derived; the spectral density is the
piece to replace with the intended bath.

beta is ramped adiabatically from hot to cold. Reported: |gamma_1-dot| including the dissipator
(stationarity), |[A, gamma_1]| (commutation with the generator), idempotency |g^2 - g|, occupation
extremes, and the distance to the Fermi-Dirac profile of the final generator.
"""
import sys
import numpy as np
sys.path.insert(0, "."); sys.path.insert(0, "tools")
from molecular_2rdm import molecule
from qdiss import closures
from qdiss.eom import gwt
from qdiss.eom.evaluate import assemble

CUM1 = gwt.terms_cumulant(1)
BASE = np.array([[0.00, 0.00, 0.0], [0.80, 0.00, 0.0],
                 [1.05, 0.72, 0.0], [-0.15, 0.72, 0.0]])
geom = BASE + 0.05 * np.random.default_rng(20260805).normal(size=BASE.shape)
atom = "; ".join("H %.6f %.6f %.6f" % tuple(r) for r in geom)
mol, mf, M, h, V, enuc = molecule(atom, "sto-3g", 0, 0)
N = mol.nelectron
EYE = np.eye(M, dtype=complex)
Z4 = np.zeros((M,) * 4, dtype=complex)
GAMMA0, DEPH = 0.05, 0.05


def gen(g1):
    return (h + np.einsum("prqs,rs->pq", V, g1)).conj()


def coherent(g1, closure):
    l2 = Z4 if closure == "lam2_0" else closures.lambda2(g1, "blocked")
    return assemble(CUM1[0], CUM1[1], h, V,
                    {("D", 1): EYE, ("G", 1): g1, ("L", 2): l2}, "pq", M)


def dissipator(g1, beta):
    A = gen(g1)
    eps, U = np.linalg.eigh(0.5 * (A + A.conj().T))
    gt = U.conj().T @ g1 @ U
    nn = np.real(np.diag(gt))
    d = eps[:, None] - eps[None, :]
    K = GAMMA0 / (1.0 + np.exp(np.clip(beta * d, -60, 60)))
    np.fill_diagonal(K, 0.0)
    gain = K * (nn[None, :] * (1.0 - nn[:, None]))
    loss = K.T * (nn[:, None] * (1.0 - nn[None, :]))
    Dt = np.zeros((M, M), dtype=complex)
    np.fill_diagonal(Dt, (gain - loss).sum(axis=1))
    off = -DEPH * gt
    np.fill_diagonal(off, 0.0)
    return U @ (Dt + off) @ U.conj().T


def rhs(g1, closure, beta):
    return coherent(g1, closure) + dissipator(g1, beta)


def rk4(f, y, dt):
    k1 = f(y); k2 = f(y + dt / 2 * k1); k3 = f(y + dt / 2 * k2); k4 = f(y + dt * k3)
    return y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def fermi(eps, nel, beta):
    lo, hi = eps.min() - 20, eps.max() + 20
    for _ in range(200):
        mu = 0.5 * (lo + hi)
        s = (1.0 / (1.0 + np.exp(np.clip(beta * (eps - mu), -60, 60)))).sum()
        if s > nel:
            hi = mu
        else:
            lo = mu
    return 1.0 / (1.0 + np.exp(np.clip(beta * (eps - mu), -60, 60)))


rng = np.random.default_rng(3)
A0 = rng.normal(size=(M, M)) + 1j * rng.normal(size=(M, M))
A0 = 0.5 * (A0 + A0.conj().T)
w0, U0 = np.linalg.eigh(A0)
occ0 = np.zeros(M); occ0[:N] = 1.0
g_start = U0 @ np.diag(occ0.astype(complex)) @ U0.conj().T

dt = 0.02
SCHED = [(2.0, 400), (5.0, 400), (10.0, 400), (25.0, 600), (60.0, 800), (150.0, 1000)]
for closure in ("lam2_0", "blocked"):
    g = g_start.copy()
    print("=== closure %s ===" % closure)
    print("   %8s %5s | %11s %11s %11s | %8s %8s %11s"
          % ("beta", "T", "|g1dot|", "|[A,g1]|", "|g^2-g|", "n_min", "n_max", "|n-FermiDirac|"))
    for beta, ns in SCHED:
        for _ in range(ns):
            g = rk4(lambda y: rhs(y, closure, beta), g, dt)
            g = 0.5 * (g + g.conj().T)
        A = gen(g)
        eps, U = np.linalg.eigh(0.5 * (A + A.conj().T))
        nn = np.sort(np.real(np.diag(U.conj().T @ g @ U)))
        fd = np.sort(fermi(eps, N, beta))
        print("   %8.1f %5.3f | %11.3e %11.3e %11.3e | %8.5f %8.5f %11.3e"
              % (beta, 1.0 / beta, np.abs(rhs(g, closure, beta)).max(),
                 np.abs(A @ g - g @ A).max(), np.abs(g @ g - g).max(),
                 nn.min(), nn.max(), np.abs(nn - fd).max()))
    print()
