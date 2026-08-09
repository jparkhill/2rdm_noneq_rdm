#!/usr/bin/env python3
"""Full excitation-energy table from an UNRESTRICTED, spin-broken reference.

Reference orbitals. Previous runs used RHF, so the one-particle basis and the ADC reference were
both restricted. Here the reference is UHF started from a spin-mixed guess and then driven to an
internally stable solution by repeated stability analysis. Only FCI is invariant to this choice;
TDHF, ADC(2) and every closure see a different one-particle basis.

Spin-orbital integrals are built from the separate alpha and beta coefficient matrices, so the
aa, ab and bb two-electron blocks are transformed independently before antisymmetrisation.

Methods, all energies in Hartree:
    FCI      exact diagonalisation in the N-electron sector (all S_z), basis-invariant
    TDHF     Jacobian of the rank-1 EOM with lambda_2 = 0 at its own UHF fixed point
    ADC(2)   pyscf.adc UADC, method adc(2), method_type ee
    lam3_0   ours: MKNO cumulant expansion with lambda_3 = 0, gamma_2 propagated alone
    muller   Mueller Eq. 14 with the Lackner contraction correction, gamma_2 propagated alone
    muller_raw   the same without the correction

The closure poles come from response diagonalisation of the gamma_2 flow at the FCI ground-state
gamma_2. The FCI ground state is a fixed point of the exact flow but not of any closure, so an
inhomogeneous term is present; it shifts no eigenvalue.
"""
import sys
import time
import numpy as np
import torch
sys.path.insert(0, "."); sys.path.insert(0, "tools")
from pyscf import ao2mo, gto, scf
from molecular_2rdm import ORB_B, build_H_sector, sector, separable_gpu
from qdiss import gpu
from qdiss.eom import fullci, wick
from qdiss.eom.evaluate import assemble as np_assemble

ONE1, TWO1 = wick.terms(1)
ONE2, TWO2 = wick.terms(2)
dev = gpu.device()
H4_TRAP = "H 0.00 0.00 0.0; H 0.80 0.00 0.0; H 1.05 0.72 0.0; H -0.15 0.72 0.0"


def uhf_stable(mol, seed=0.3):
    mf = scf.UHF(mol)
    dm = mf.get_init_guess()
    dm = (dm[0].copy(), dm[1].copy())
    dm[0][0, 0] += seed
    dm[1][0, 0] -= seed
    mf.kernel(dm0=dm)
    for _ in range(6):
        mo = mf.stability()[0]
        if isinstance(mo, (list, tuple)):
            same = all(np.allclose(a, b) for a, b in zip(mo, mf.mo_coeff))
        else:
            same = np.allclose(mo, mf.mo_coeff)
        if same:
            break
        mf = mf.run(mf.make_rdm1(mo, mf.mo_occ))
    return mf


def spin_orbital_uhf(mol, mf):
    Ca, Cb = mf.mo_coeff
    n = Ca.shape[1]
    M = 2 * n
    hao = mf.get_hcore()
    ha, hb = Ca.T @ hao @ Ca, Cb.T @ hao @ Cb
    h = np.zeros((M, M), dtype=complex)
    for p in range(n):
        for q in range(n):
            h[2 * p, 2 * q] = ha[p, q]
            h[2 * p + 1, 2 * q + 1] = hb[p, q]
    eaa = ao2mo.restore(1, ao2mo.general(mol, (Ca, Ca, Ca, Ca), compact=False), n)
    ebb = ao2mo.restore(1, ao2mo.general(mol, (Cb, Cb, Cb, Cb), compact=False), n)
    eab = ao2mo.general(mol, (Ca, Ca, Cb, Cb), compact=False).reshape(n, n, n, n)
    phys = np.zeros((M,) * 4)
    blocks = {(0, 0): eaa, (1, 1): ebb, (0, 1): eab,
              (1, 0): np.transpose(eab, (2, 3, 0, 1))}
    for p in range(n):
        for q in range(n):
            for r in range(n):
                for s in range(n):
                    for sp in (0, 1):
                        for sq in (0, 1):
                            phys[2 * p + sp, 2 * q + sq, 2 * r + sp, 2 * s + sq] += \
                                blocks[(sp, sq)][p, r, q, s]
    V = (phys - np.transpose(phys, (0, 1, 3, 2))).astype(complex)
    return M, h, V


mol = gto.M(atom=H4_TRAP, basis="sto-3g", verbose=0)
N = mol.nelectron
mfr = scf.RHF(mol).run()
mfu = uhf_stable(mol)
print("RHF  E = %.9f" % mfr.e_tot)
print("UHF  E = %.9f   <S^2> = %.4f   (spin broken by %.2e)"
      % (mfu.e_tot, mfu.spin_square()[0], abs(mfu.e_tot - mfr.e_tot)))
M, h, V = spin_orbital_uhf(mol, mfu)
enuc = float(mol.energy_nuc())
dets = sector(M, N)
t0 = time.perf_counter()
es, vec = np.linalg.eigh(build_H_sector(h, V, M, dets))
exc = np.sort((es - es[0]).real)
exc = exc[exc > 1e-8][:10]
print("FCI  E = %.9f   (basis-invariant check)   sector %d, %.0fs"
      % (es[0].real + enuc, len(dets), time.perf_counter() - t0))
print()
res = {"FCI": exc}

from pyscf import adc
try:
    myadc = adc.ADC(mfu)
    myadc.method = "adc(2)"
    myadc.method_type = "ee"
    myadc.verbose = 0
    out = myadc.kernel(nroots=10)
    res["ADC(2)"] = np.sort(np.atleast_1d(out[0]))[:10]
except Exception as ex:
    print("ADC(2) failed:", type(ex).__name__, str(ex)[:160])

EYE = np.eye(M, dtype=complex)


def herm_basis(dim):
    B = []
    for i in range(dim):
        E = np.zeros((dim, dim), dtype=complex); E[i, i] = 1.0; B.append(E)
    for i in range(dim):
        for j in range(i + 1, dim):
            E = np.zeros((dim, dim), dtype=complex); E[i, j] = 1.0; E[j, i] = 1.0; B.append(E)
            E = np.zeros((dim, dim), dtype=complex); E[i, j] = 1j; E[j, i] = -1j; B.append(E)
    return B


def fock(g):
    return h + np.einsum("prqs,rs->pq", V, g)


g = np.zeros((M, M), dtype=complex)
for _ in range(600):
    _, vv = np.linalg.eigh(0.5 * (fock(g).conj() + fock(g).conj().conj().T))
    new = vv[:, :N] @ vv[:, :N].conj().T
    if np.abs(new - g).max() < 1e-13:
        g = new; break
    g = 0.6 * new + 0.4 * g


def d1(x):
    g2 = np.einsum("pq,vw->pvqw", x, x) - np.einsum("pw,vq->pvqw", x, x)
    return np_assemble(ONE1, TWO1, h, V,
                       {("D", 1): EYE, ("G", 1): x, ("G", 2): g2}, "pq", M)


B1 = herm_basis(M); n1 = len(B1); J = np.zeros((n1, n1)); f0 = d1(g); eps = 1e-6
for k in range(n1):
    dp = (d1(g + eps * B1[k]) - f0) / eps
    for m in range(n1):
        J[m, k] = float(np.real(np.trace(B1[m].conj().T @ dp))) / float(
            np.real(np.trace(B1[m].conj().T @ B1[m])))
w = np.linalg.eigvals(J)
res["TDHF"] = np.sort(np.abs(w.imag[w.imag > 1e-6]))[:10]
print("TDHF at its UHF fixed point: |gamma_1-dot| = %.1e" % np.abs(f0).max())

psi0 = np.zeros(1 << M, dtype=complex); psi0[dets] = vec[:, 0]
g2_0 = fullci.rdm(psi0, M, 2)
hg, Vg, Ig = gpu.up(h, dev), gpu.up(V, dev), gpu.up(np.eye(M), dev)
PAIRS = [(p, v) for p in range(M) for v in range(p + 1, M)]
NP = len(PAIRS)


def pair_to_tensor(Z):
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
            B2.append(pair_to_tensor(Z))
n2 = len(B2)
Bg = [gpu.up(b, dev) for b in B2]
nrm = [float(np.real(np.tensordot(b.conj(), b, axes=4))) for b in B2]


def recon(nm, g1, g2):
    if nm == "muller2405":
        return gpu.muller2405_gamma3(g1, g2, N, M, True, ORB_B)
    if nm == "muller2405_raw":
        return gpu.muller2405_gamma3(g1, g2, N, M, False, ORB_B)
    return separable_gpu(g1, g2, M)


def rhs(x, nm):
    g1 = torch.einsum("pxqx->pq", x) / (N - 1)
    return gpu.assemble(ONE2, TWO2, hg, Vg,
                        {("D", 1): Ig, ("G", 1): g1, ("G", 2): x,
                         ("G", 3): recon(nm, g1, x)}, "pvqw")


for nm in ("lam3_0", "muller2405", "muller2405_raw"):
    f0 = rhs(gpu.up(g2_0, dev), nm)
    Jm = np.zeros((n2, n2))
    for k in range(n2):
        dn = gpu.down((rhs(gpu.up(g2_0, dev) + eps * Bg[k], nm) - f0) / eps)
        for m in range(n2):
            Jm[m, k] = float(np.real(np.tensordot(B2[m].conj(), dn, axes=4))) / nrm[m]
    ww = np.linalg.eigvals(Jm)
    res[nm] = np.sort(np.abs(ww.imag[ww.imag > 1e-6]))[:10]

print()
order = ["FCI", "TDHF", "ADC(2)", "lam3_0", "muller2405", "muller2405_raw"]
print("Excitation energies, Hartree, sorted ascending, UHF spin-broken reference")
print("%-11s | %s" % ("method", "  ".join("%7d" % (i + 1) for i in range(10))))
for k in order:
    if k not in res:
        continue
    a = np.asarray(res[k])[:10]
    print("%-11s | %s" % (k, "  ".join("%7.4f" % x for x in a)))
print()
print("%-11s | %9s" % ("method", "mean |err|"))
for k in order[1:]:
    if k not in res:
        continue
    a = np.asarray(res[k])
    m = [a[np.argmin(np.abs(a - e))] for e in exc]
    print("%-11s | %9.4f" % (k, np.mean(np.abs(np.array(m) - exc))))
