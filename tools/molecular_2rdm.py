#!/usr/bin/env python3
"""2-RDM truncated comparison on molecular systems, GPU, against exact propagation.

Derivation: qdiss.eom.wick term lists (identical to qdiss.eom.gwt), validated against
qdiss.eom.fullci at 1e-15; contracted on device by qdiss.gpu, itself validated against the numpy
path at complex128 round-off in tools/gpu_validate.py. Propagated state: gamma_2 ALONE with
gamma_1 by contraction, which is the Lackner/Mueller scheme.

Integrals from PySCF in the convention of qdiss.eom.hamiltonian, each Hamiltonian checked against
the PySCF FCI total energy. The exact reference is built in the N-electron SECTOR rather than the
full Fock space, so the dense matrix is C(M,N) square instead of 2^M square; that is what makes
M = 12 reachable.

Closures compared, all producing gamma_3 from gamma_1 and gamma_2:
    muller      Mueller Eq. 14 plus the Lackner contraction correction
    muller_raw  the same without the correction
    lam3_0      ours, gamma_3 = the separable part (lambda_3 = 0)
"""
import itertools
import sys
import time
import numpy as np
import torch
sys.path.insert(0, ".")
from pyscf import ao2mo, fci, gto, scf
from qdiss import gpu
from qdiss.eom import fullci, wick
from qdiss.symbolic import mkno

UP3, LO3 = ("p", "v", "x"), ("q", "w", "y")
ORB_B = [(s, ",".join("".join(u) + "".join(l) for n, r, u, l in f) + "->pvxqwy",
          [(n, r) for n, r, u, l in f])
         for s, f in mkno.gamma_to_cumulants(UP3, LO3)
         if not any(n == "L" and r == 3 for n, r, u, l in f)
         and any(n == "L" and r == 2 for n, r, u, l in f)]
SEP3 = [(s, ",".join("".join(u) + "".join(l) for n, r, u, l in f) + "->pvxqwy",
         [(n, r) for n, r, u, l in f])
        for s, f in mkno.gamma_to_cumulants(UP3, LO3)
        if not any(n == "L" and r == 3 for n, r, u, l in f)]
ONE2, TWO2 = wick.terms(2)
dev = gpu.device()


def molecule(atom, basis, charge=0, spin=0):
    mol = gto.M(atom=atom, basis=basis, charge=charge, spin=spin, verbose=0)
    mf = (scf.ROHF(mol) if spin else scf.RHF(mol)).run()
    C = mf.mo_coeff
    n = C.shape[1]
    hc = C.T @ mf.get_hcore() @ C
    eri = ao2mo.restore(1, ao2mo.kernel(mol, C), n)
    M = 2 * n
    h = np.zeros((M, M), dtype=complex)
    for p in range(n):
        for q in range(n):
            h[2 * p, 2 * q] = hc[p, q]
            h[2 * p + 1, 2 * q + 1] = hc[p, q]
    phys = np.zeros((M,) * 4)
    for p in range(n):
        for q in range(n):
            for r in range(n):
                for s in range(n):
                    v = eri[p, r, q, s]
                    for sp in (0, 1):
                        for sq in (0, 1):
                            phys[2 * p + sp, 2 * q + sq, 2 * r + sp, 2 * s + sq] += v
    V = (phys - np.transpose(phys, (0, 1, 3, 2))).astype(complex)
    return mol, mf, M, h, V, float(mol.energy_nuc())


def sector(M, N):
    return [d for d in range(1 << M) if bin(d).count("1") == N]


def build_H_sector(h, V, M, dets):
    """H restricted to the N-electron sector, dim C(M,N)."""
    idx = {d: i for i, d in enumerate(dets)}
    n = len(dets)
    H = np.zeros((n, n), dtype=complex)
    basis = np.zeros(1 << M, dtype=complex)
    hnz = [(p, q) for p in range(M) for q in range(M) if abs(h[p, q]) > 1e-14]
    vnz = [tuple(x) for x in np.argwhere(np.abs(V) > 1e-14)]
    for j, d in enumerate(dets):
        basis[:] = 0.0
        basis[d] = 1.0
        col = np.zeros(1 << M, dtype=complex)
        for p, q in hnz:
            col += h[p, q] * fullci.ad(p, fullci.aa(q, basis))
        for p, t, q, u in vnz:
            col += 0.25 * V[p, t, q, u] * fullci.ad(p, fullci.ad(t, fullci.aa(u, fullci.aa(q, basis))))
        for d2, i in idx.items():
            if col[d2] != 0:
                H[i, j] = col[d2]
    return 0.5 * (H + H.conj().T)


def wedge_np(a, b):
    return np.einsum("pq,vw->pvqw", a, b) - np.einsum("pw,vq->pvqw", a, b)


def separable_gpu(g1, g2, M):
    l2 = g2 - gpu.wedge(g1, g1)
    t = {("G", 1): g1, ("L", 2): l2}
    out = None
    for s, sub, keys in SEP3:
        piece = s * torch.einsum(sub, *[t[k] for k in keys])
        out = piece if out is None else out + piece
    return out


NAMES = ("muller2405", "muller2405_raw", "lam3_0")
_DRIVER_SYS = (("H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0", "sto-3g", 4, 0, 0, "H4 chain / sto-3g"),
       ("Li 0 0 0; H 0 0 1.6", "sto-3g", 4, 0, 0, "LiH / sto-3g"))

for atom, basis, N, chg, spin, label in (_DRIVER_SYS if __name__ == "__main__" else ()):
    mol, mf, M, h, V, enuc = molecule(atom, basis, chg, spin)
    dets = sector(M, N)
    t0 = time.perf_counter()
    Hs = build_H_sector(h, V, M, dets)
    es, vec = np.linalg.eigh(Hs)
    efci = fci.FCI(mf).kernel()[0]
    print("=== %s   M=%d N=%d   sector dim %d ===" % (label, M, N, len(dets)))
    print("   ground E %.9f vs PySCF FCI %.9f  (diff %.2e)   H build %.1fs"
          % (es[0].real + enuc, efci, abs(es[0].real + enuc - efci),
             time.perf_counter() - t0))
    rr = np.random.default_rng(7)
    c = vec[:, 0] + 0.3 * (rr.normal(size=len(dets)) + 1j * rr.normal(size=len(dets)))
    c /= np.linalg.norm(c)
    psi0 = np.zeros(1 << M, dtype=complex)
    psi0[dets] = c
    g2_0 = fullci.rdm(psi0, M, 2)
    hg, Vg, Ig = gpu.up(h, dev), gpu.up(V, dev), gpu.up(np.eye(M), dev)

    def recon(nm, g1, g2):
        if nm == "muller2405":
            return gpu.muller2405_gamma3(g1, g2, N, M, True, ORB_B)
        if nm == "muller2405_raw":
            return gpu.muller2405_gamma3(g1, g2, N, M, False, ORB_B)
        if nm == "lam3_0":
            return separable_gpu(g1, g2, M)
        raise ValueError("unknown reconstruction %r; moment negation is excluded" % nm)

    def rhs(x, nm):
        g1 = torch.einsum("pxqx->pq", x) / (N - 1)
        g3 = recon(nm, g1, x)
        return gpu.assemble(ONE2, TWO2, hg, Vg,
                            {("D", 1): Ig, ("G", 1): g1, ("G", 2): x, ("G", 3): g3}, "pvqw")

    def rk4(f, y, dt):
        k1 = f(y); k2 = f(y + dt / 2 * k1); k3 = f(y + dt / 2 * k2); k4 = f(y + dt * k3)
        return y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

    PAIRS = [(p, v) for p in range(M) for v in range(p + 1, M)]

    def to_pair(T):
        A = np.array([[2.0 * T[p, v, q, w] for (q, w) in PAIRS] for (p, v) in PAIRS])
        return 0.5 * (A + A.conj().T)

    st = {nm: gpu.up(g2_0, dev) for nm in NAMES}
    dt, ns, ev = 0.005, 200, 50
    tstart = time.perf_counter()
    print("   |gamma_2-exact| / min eig D / n_min(gamma_1)")
    print("   %5s | %s" % ("t", " | ".join("%-22s" % nm for nm in NAMES)))
    for step in range(ns + 1):
        if step % ev == 0:
            t = step * dt
            ct = vec @ (np.exp(-1j * es.real * t) * (vec.conj().T @ c))
            pt = np.zeros(1 << M, dtype=complex)
            pt[dets] = ct
            ex = fullci.rdm(pt, M, 2)
            cells = []
            for nm in NAMES:
                g2 = gpu.down(st[nm])
                if not np.isfinite(g2).all():
                    cells.append("%-22s" % "diverged"); continue
                g1 = np.einsum("pxqx->pq", g2) / (N - 1)
                cells.append("%7.4f %7.4f %6.3f"
                             % (np.abs(g2 - ex).max(),
                                np.linalg.eigvalsh(to_pair(g2)).min(),
                                np.linalg.eigvalsh(0.5 * (g1 + g1.conj().T)).min()))
            print("   %5.2f | %s" % (t, " | ".join(cells)))
        if step < ns:
            for nm in NAMES:
                if torch.isfinite(st[nm]).all():
                    st[nm] = rk4(lambda y, k=nm: rhs(y, k), st[nm], dt)
    if dev.type == "cuda":
        torch.cuda.synchronize()
    print("   propagation wall time %.1f s for %d steps x %d closures"
          % (time.perf_counter() - tstart, ns, len(NAMES)))
    print()
