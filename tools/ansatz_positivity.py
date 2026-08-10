#!/usr/bin/env python3
"""Positivity, trace and accuracy of each rank-3 ansatz IN PROPAGATION, read against the exact.

THE APPROXIMATION, identical for every row. The propagated equation is the rank-2 Heisenberg EOM
    d/dt gamma_2[pv,qw] = +i <[H, a+_p a+_v a_w a_q]>
with term list from qdiss.eom.wick, term-identical to the MKNO route qdiss.eom.gwt, validated
against exact Fock-space linear algebra qdiss.eom.fullci at 1e-15. The PROPAGATED STATE VARIABLE
is gamma_2 alone, M^4 complex numbers; gamma_1 is NOT propagated but recovered at every RK4 stage
by contraction, gamma_1 = Tr_2(gamma_2)/(N-1). This is the Lackner / Mueller scheme. The ONLY
approximation is the reconstruction of the 3-cumulant lambda_3, from which
    gamma_3 = separable(gamma_1, lambda_2) + lambda_3
with `separable` the MKNO cumulant expansion of gamma_3 with lambda_3 deleted (qdiss.ansatz).
Reconstructions, defined in qdiss.ansatz and scored statically in tools/ansatz_hierarchy.py:
    null          lambda_3 = 0
    scaled        Mueller arXiv:2405.05093 Eq. 14, closed-system limit, no Lackner correction
    trace         minimum-norm lambda_3 with Tr_3 lambda_3 = anti(E, lambda_2), E = gamma_1 - I/2
    scaled+trace  Mueller Eq. 14 WITH the Lackner Eq. 43/44/45 contraction correction

INITIAL STATES. Per CLAUDE.md the wavefunction is perturbed, never the RDM: a direct kick to
gamma_2 is not N-representable and would dominate every positivity metric. Here
    pure:    |Psi(0)> = exp(-i eps mu_z) |Psi_FCI>,   eps = 0.05
    mixed:   rho(0)   = K rho_th K^dag,  K = exp(-i eps mu_z),  rho_th propto exp(-H/kT), kT = 0.5
with mu_z = sum_pq d[p,q] a+_p a_q the z dipole in the MO basis, built in the N-electron sector.

METRICS, all read AGAINST THE EXACT VALUE at the same time, never against zero.
    D[(p,v),(q,w)] = 2 gamma_2[p,v,q,w] over unordered pairs p<v: the pair matrix, positive
                     semidefinite for any N-representable gamma_2. min eig D is the CLEARANCE.
                     The exact clearance is small and can sit on the boundary for a pure state,
                     so the reported figure is the closure's clearance next to the exact one.
    n_min, n_max     extremal eigenvalues of gamma_1, which must lie in [0,1] (the Pauli bound).
    |tr g1 - N|, |tr g2 - N(N-1)|   trace errors.
    |g2 - exact|     max element deviation of gamma_2.

System: H4 asymmetric trapezoid / sto-3g, M = 8, N = 4. No spatial symmetry, N >= 3.
"""
import sys
import time
import numpy as np
import scipy.linalg as sla
import torch
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from ansatz_hierarchy import H4_TRAP, rank3_set
from molecular_2rdm import build_H_sector, molecule, sector
from pyscf import ao2mo
from qdiss import ansatz as az, gpu
from qdiss.eom import fullci, wick

ONE2, TWO2 = wick.terms(2)
dev = gpu.device()
N, EPS, KT = 4, 0.05, 0.5
DT, NS, EVERY = 0.005, 200, 50
NAMES = ("null", "scaled", "trace", "scaled+trace")

mol, mf, M, h, V, enuc = molecule(H4_TRAP, "sto-3g", 0, 0)
dets = sector(M, N)
Hs = build_H_sector(h, V, M, dets)
es, vec = np.linalg.eigh(Hs)

# z dipole in the spin-orbital MO basis, then in the N-electron sector
dz_ao = mol.intor("int1e_r")[2]
dz_mo = mf.mo_coeff.T @ dz_ao @ mf.mo_coeff
nsp = dz_mo.shape[0]
d = np.zeros((M, M), dtype=complex)
for p in range(nsp):
    for q in range(nsp):
        d[2 * p, 2 * q] = dz_mo[p, q]
        d[2 * p + 1, 2 * q + 1] = dz_mo[p, q]
idx = {dd: i for i, dd in enumerate(dets)}
Mu = np.zeros((len(dets),) * 2, dtype=complex)
for j, dd in enumerate(dets):
    b = np.zeros(1 << M, dtype=complex)
    b[dd] = 1.0
    col = np.zeros(1 << M, dtype=complex)
    for p in range(M):
        for q in range(M):
            if abs(d[p, q]) > 1e-14:
                col += d[p, q] * fullci.ad(p, fullci.aa(q, b))
    for d2, i in idx.items():
        if col[d2] != 0:
            Mu[i, j] = col[d2]
Mu = 0.5 * (Mu + Mu.conj().T)
K = sla.expm(-1j * EPS * Mu)

PAIRS = [(p, v) for p in range(M) for v in range(p + 1, M)]
hg, Vg, Ig = gpu.up(h, dev), gpu.up(V, dev), gpu.up(np.eye(M), dev)


def rdm12(dm):
    """gamma_1, gamma_2 of a sector density matrix, dropping weights below 1e-9."""
    w, U = np.linalg.eigh(dm)
    g1 = np.zeros((M,) * 2, dtype=complex)
    g2 = np.zeros((M,) * 4, dtype=complex)
    for i, wi in enumerate(w):
        if wi < 1e-9:
            continue
        p = np.zeros(1 << M, dtype=complex)
        p[dets] = U[:, i]
        g1 += wi * fullci.rdm(p, M, 1)
        g2 += wi * fullci.rdm(p, M, 2)
    return g1, g2


def pair_matrix(g2):
    A = np.array([[2.0 * g2[p, v, q, w] for (q, w) in PAIRS] for (p, v) in PAIRS])
    return 0.5 * (A + A.conj().T)


def recon3(nm, g1, g2):
    l2 = g2 - az.wedge(g1, g1)
    return az.separable(3, {("G", 1): g1, ("L", 2): l2}) + rank3_set(g1, l2, N, M)[nm]


def rhs(x, nm):
    g2 = gpu.down(x)
    g1 = np.einsum("pxqx->pq", g2) / (N - 1)
    g3 = gpu.up(recon3(nm, g1, g2), dev)
    return gpu.assemble(ONE2, TWO2, hg, Vg,
                        {("D", 1): Ig, ("G", 1): gpu.up(g1, dev), ("G", 2): x,
                         ("G", 3): g3}, "pvqw")


def rk4(f, y, dt):
    k1 = f(y); k2 = f(y + dt / 2 * k1); k3 = f(y + dt / 2 * k2); k4 = f(y + dt * k3)
    return y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


dm_pure = np.zeros((len(dets),) * 2, dtype=complex)
dm_pure[0, 0] = 1.0
w = np.exp(-(es - es[0]).real / KT)
w /= w.sum()
STATES = (("pure FCI ground, dipole-kicked eps=%.2f" % EPS, K @ dm_pure @ K.conj().T),
          ("thermal kT=%.1f, dipole-kicked eps=%.2f" % (KT, EPS),
           K @ np.diag(w.astype(complex)) @ K.conj().T))

print("H4 asymmetric trapezoid / sto-3g   M=%d N=%d   sector %d   E_FCI %.9f"
      % (M, N, len(dets), es[0].real + enuc))
print("propagated: gamma_2 alone, gamma_1 = Tr_2(gamma_2)/(N-1); only lambda_3 approximated")
print()

for sname, dm0 in STATES:
    dm0e = vec @ dm0 @ vec.conj().T          # sector eigenbasis -> determinant basis
    g1_0, g2_0 = rdm12(dm0e)
    print("=== initial state: %s ===" % sname)
    print("    exact at t=0:  min eig D %+.6f   occupations [%.4f .. %.4f]"
          % (np.linalg.eigvalsh(pair_matrix(g2_0)).min(),
             *np.linalg.eigvalsh(0.5 * (g1_0 + g1_0.conj().T))[[0, -1]].real))
    st = {nm: gpu.up(g2_0, dev) for nm in NAMES}
    print("    %-13s %5s | %9s | %10s %10s | %11s %11s | %9s %9s"
          % ("ansatz", "t", "|g2-exact|", "|tr g1-N|", "|tr g2-NN|",
             "min eig D", "exact D", "n_min", "n_max"))
    t0 = time.perf_counter()
    for step in range(NS + 1):
        if step % EVERY == 0:
            t = step * DT
            ph = np.exp(-1j * es.real * t)
            dmt = vec @ (ph[:, None] * dm0 * ph.conj()[None, :]) @ vec.conj().T
            g1e, g2e = rdm12(dmt)
            de = np.linalg.eigvalsh(pair_matrix(g2e)).min()
            for nm in NAMES:
                g2 = gpu.down(st[nm])
                if not np.isfinite(g2).all():
                    print("    %-13s %5.2f | %s" % (nm, t, "diverged")); continue
                g1 = np.einsum("pxqx->pq", g2) / (N - 1)
                occ = np.linalg.eigvalsh(0.5 * (g1 + g1.conj().T))
                print("    %-13s %5.2f | %9.2e | %10.2e %10.2e | %+11.6f %+11.6f | %9.6f %9.6f"
                      % (nm, t, np.abs(g2 - g2e).max(), abs(np.trace(g1).real - N),
                         abs(np.einsum("pvpv->", g2).real - N * (N - 1)),
                         np.linalg.eigvalsh(pair_matrix(g2)).min(), de,
                         occ.min(), occ.max()))
            print("    %-13s %5.2f | %9s | %10s %10s | %+11.6f %+11.6f | %9.6f %9.6f"
                  % ("exact", t, "--", "--", "--", de, de,
                     *np.linalg.eigvalsh(0.5 * (g1e + g1e.conj().T))[[0, -1]].real))
            print()
        if step < NS:
            for nm in NAMES:
                if torch.isfinite(st[nm]).all():
                    st[nm] = rk4(lambda y, k=nm: rhs(y, k), st[nm], DT)
    print("    wall time %.1f s\n" % (time.perf_counter() - t0))
