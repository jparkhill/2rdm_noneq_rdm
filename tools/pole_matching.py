#!/usr/bin/env python3
"""Assign closure poles to FCI states by response fingerprint, not by proximity.

Sorting Jacobian eigenvalues and reading off the lowest is wrong: the gamma_2 Jacobian has many
more modes than there are physical states, and the extras are redundant directions (occupied-
occupied and virtual-virtual rotations) carrying no transition density. Matching each FCI energy
to its nearest pole is equally wrong, since with hundreds of candidates something is always close.

Scheme, per John: fingerprint every mode by its response to several random Hermitian one-body
observables, then assign correspondence by maximum cosine between fingerprint vectors. Three
refinements applied:
  (a) magnitude gate first. Cosine is scale-free, so a noise mode with a random direction can
      score well; modes whose total response is below a fraction of the largest are discarded
      before any matching. This is what actually removes the redundant rotations.
  (b) phase-insensitive cosine |<a,b>|/(|a||b|), since eigenvectors carry an arbitrary phase.
  (c) degenerate FCI levels are matched as blocks. Within a multiplet the eigenvectors are only
      defined up to a unitary, so individual fingerprints are meaningless there.

Fingerprint of a Jacobian mode: contract its delta-gamma_2 to delta-gamma_1 and take
sum_pq O[p,q] delta-gamma_1[p,q]. Fingerprint of FCI state k: <0|O_hat|k> in the sector.

Geometry: the asymmetric trapezoid with 0.05 * normal jitter on every coordinate. Note this
breaks SPATIAL degeneracy only; spin multiplets survive any spatial perturbation because H
commutes with S^2.
"""
import sys
import numpy as np
import torch
sys.path.insert(0, "."); sys.path.insert(0, "tools")
from pyscf import gto
from molecular_2rdm import ORB_B, build_H_sector, sector, separable_gpu
from qdiss import gpu
from qdiss.eom import fullci, gwt, wick
from qdiss.eom.evaluate import assemble as np_assemble
from pyscf import tdscf
from unrestricted_spectra import spin_orbital_uhf, uhf_stable

ONE2, TWO2 = wick.terms(2)
dev = gpu.device()
BASE = np.array([[0.00, 0.00, 0.0], [0.80, 0.00, 0.0],
                 [1.05, 0.72, 0.0], [-0.15, 0.72, 0.0]])
rj = np.random.default_rng(20260805)
geom = BASE + 0.05 * rj.normal(size=BASE.shape)
atom = "; ".join("H %.6f %.6f %.6f" % tuple(r) for r in geom)
mol = gto.M(atom=atom, basis="sto-3g", verbose=0)
N = mol.nelectron
mfu = uhf_stable(mol)
M, h, V = spin_orbital_uhf(mol, mfu)[0:3]
dets = sector(M, N)
es, vec = np.linalg.eigh(build_H_sector(h, V, M, dets))
exc = (es - es[0]).real
print("H4 trapezoid + 0.05*normal jitter, UHF reference, M=%d N=%d" % (M, N))
print("FCI excitations (Ha): %s" % np.array2string(np.sort(exc)[1:11], precision=4))

NOBS = 6
ro = np.random.default_rng(11)
OBS = []
for _ in range(NOBS):
    A = ro.normal(size=(M, M)) + 1j * ro.normal(size=(M, M))
    OBS.append(0.5 * (A + A.conj().T))

idx = {d: i for i, d in enumerate(dets)}
nd = len(dets)
Osec = []
for O in OBS:
    Om = np.zeros((nd, nd), dtype=complex)
    b = np.zeros(1 << M, dtype=complex)
    for j, d in enumerate(dets):
        b[:] = 0; b[d] = 1.0
        col = np.zeros(1 << M, dtype=complex)
        for p in range(M):
            for q in range(M):
                if abs(O[p, q]) > 1e-14:
                    col += O[p, q] * fullci.ad(p, fullci.aa(q, b))
        for d2, i in idx.items():
            Om[i, j] = col[d2]
    Osec.append(Om)
fci_fp = np.array([[np.vdot(vec[:, 0], Om @ vec[:, k]) for Om in Osec] for k in range(nd)])
fci_mag = np.linalg.norm(fci_fp, axis=1)
bright = [k for k in range(1, nd) if fci_mag[k] > 0.02 * fci_mag[1:].max()]
groups, used = [], set()
for k in bright:
    if k in used:
        continue
    grp = [j for j in bright if abs(exc[j] - exc[k]) < 1e-6]
    used.update(grp)
    groups.append(grp)
groups = sorted(groups, key=lambda g: exc[g[0]])[:6]
print("bright FCI levels (mag > 2%% of max): %d, grouped into %d degenerate blocks"
      % (len(bright), len(groups)))

hg, Vg, Ig = gpu.up(h, dev), gpu.up(V, dev), gpu.up(np.eye(M), dev)
PAIRS = [(p, v) for p in range(M) for v in range(p + 1, M)]
NP = len(PAIRS)


def p2t(Z):
    T = np.zeros((M,) * 4, dtype=complex)
    for i, (p, v) in enumerate(PAIRS):
        for j, (q, w) in enumerate(PAIRS):
            x = 0.5 * Z[i, j]
            T[p, v, q, w] = x; T[v, p, q, w] = -x
            T[p, v, w, q] = -x; T[v, p, w, q] = x
    return T


B2 = []
for a in range(NP):
    for b in range(a, NP):
        for ph in (1.0, 1j):
            Z = np.zeros((NP, NP), dtype=complex); Z[a, b] += ph; Z[b, a] += np.conj(ph)
            if np.abs(Z).max() < 1e-12:
                continue
            B2.append(p2t(Z))
n2 = len(B2)
Bg = [gpu.up(b, dev) for b in B2]
nrm = [float(np.real(np.tensordot(b.conj(), b, axes=4))) for b in B2]
psi0 = np.zeros(1 << M, dtype=complex); psi0[dets] = vec[:, 0]
g2_0 = fullci.rdm(psi0, M, 2)


def recon(nm, g1, g2):
    if nm == "muller":
        return gpu.muller_gamma3(g1, g2, N, M, True, ORB_B)
    if nm == "muller_raw":
        return gpu.muller_gamma3(g1, g2, N, M, False, ORB_B)
    return separable_gpu(g1, g2, M)


def rhs(x, nm):
    g1 = torch.einsum("pxqx->pq", x) / (N - 1)
    return gpu.assemble(ONE2, TWO2, hg, Vg,
                        {("D", 1): Ig, ("G", 1): g1, ("G", 2): x,
                         ("G", 3): recon(nm, g1, x)}, "pvqw")


eps = 1e-6


def match(oms, fps, label, ntot):
    """Magnitude-gate, then phase-insensitive cosine against each degenerate FCI block."""
    mag = np.linalg.norm(fps, axis=1)
    live = mag > 0.05 * mag.max()
    row = []
    for g in groups:
        tgt = fci_fp[g].sum(axis=0)
        tn = np.linalg.norm(tgt)
        best, bw = -1.0, np.nan
        for i in np.where(live)[0]:
            cos = abs(np.vdot(tgt, fps[i])) / (tn * mag[i] + 1e-30)
            if cos > best:
                best, bw = cos, oms[i]
        row.append(bw)
    err = np.mean([abs(r - exc[g[0]]) for r, g in zip(row, groups)])
    print("%-12s | %s   mean|err| %.4f  (%d of %d modes bright)"
          % (label, "  ".join("%8.4f" % r for r in row), err, int(live.sum()), ntot))


print()
print("%-12s | %s" % ("method", "  ".join("%8.4f" % exc[g[0]] for g in groups)))
print("%-12s | %s" % ("FCI", "  ".join("%8s" % "--" for g in groups)))

td = tdscf.TDHF(mfu); td.nstates = 14; td.verbose = 0; td.kernel()
tdp = np.sort(np.asarray(td.e).real)
print("%-12s | %s   (reference list, not fingerprint-matched)"
      % ("TDHF(pyscf)", "  ".join("%8.4f" % x for x in tdp[:len(groups)])))

CUM1 = gwt.terms_cumulant(1)
EYE = np.eye(M, dtype=complex)
ZERO4 = np.zeros((M,) * 4, dtype=complex)


def d1_lam2zero(x):
    return np_assemble(CUM1[0], CUM1[1], h, V,
                       {("D", 1): EYE, ("G", 1): x, ("L", 2): ZERO4}, "pq", M)


def fock(g):
    return h + np.einsum("prqs,rs->pq", V, g)


gfp = np.zeros((M, M), dtype=complex)
for _ in range(800):
    _, vv = np.linalg.eigh(0.5 * (fock(gfp).conj() + fock(gfp).conj().conj().T))
    nw = vv[:, :N] @ vv[:, :N].conj().T
    if np.abs(nw - gfp).max() < 1e-13:
        gfp = nw; break
    gfp = 0.6 * nw + 0.4 * gfp
f0 = d1_lam2zero(gfp)
B1 = []
for i in range(M):
    E = np.zeros((M, M), dtype=complex); E[i, i] = 1.0; B1.append(E)
for i in range(M):
    for j in range(i + 1, M):
        E = np.zeros((M, M), dtype=complex); E[i, j] = 1.0; E[j, i] = 1.0; B1.append(E)
        E = np.zeros((M, M), dtype=complex); E[i, j] = 1j; E[j, i] = -1j; B1.append(E)
n1 = len(B1); J1 = np.zeros((n1, n1))
for k in range(n1):
    dp = (d1_lam2zero(gfp + eps * B1[k]) - f0) / eps
    for m in range(n1):
        J1[m, k] = float(np.real(np.trace(B1[m].conj().T @ dp))) / float(
            np.real(np.trace(B1[m].conj().T @ B1[m])))
w1, U1 = np.linalg.eig(J1)
keep1 = np.where(w1.imag > 1e-6)[0]
fp1, om1 = [], []
for c in keep1:
    dg1 = sum(U1[m, c] * B1[m] for m in range(n1))
    fp1.append([np.einsum("pq,pq->", O, dg1) for O in OBS])
    om1.append(w1[c].imag)
print("   tdhf(claude) fixed point |gamma_1-dot| = %.2e" % np.abs(f0).max())
match(np.array(om1), np.array(fp1), "tdhf(claude)", len(om1))
for nm in ("lam3_0", "muller", "muller_raw"):
    f0 = rhs(gpu.up(g2_0, dev), nm)
    J = np.zeros((n2, n2))
    for k in range(n2):
        dn = gpu.down((rhs(gpu.up(g2_0, dev) + eps * Bg[k], nm) - f0) / eps)
        for m in range(n2):
            J[m, k] = float(np.real(np.tensordot(B2[m].conj(), dn, axes=4))) / nrm[m]
    w, Uj = np.linalg.eig(J)
    keep = np.where(w.imag > 1e-6)[0]
    fps, oms = [], []
    for c in keep:
        dg2 = sum(Uj[m, c] * B2[m] for m in range(n2))
        dg1 = np.einsum("pxqx->pq", dg2) / (N - 1)
        fps.append([np.einsum("pq,pq->", O, dg1) for O in OBS])
        oms.append(w[c].imag)
    match(np.array(oms), np.array(fps), nm, len(oms))
