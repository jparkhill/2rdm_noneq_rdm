#!/usr/bin/env python3
"""Singlet-restricted pole table with Hungarian assignment.

Three fixes over tools/pole_matching.py, which failed its own acceptance test:
  (1) SINGLETS ONLY. FCI is diagonalised in the S_z = 0 sector and filtered by <S^2> ~ 0.
      tdhf(claude) uses a spin-ADAPTED perturbation basis, B[2p,2q] = B[2p+1,2q+1] = b[p,q],
      which for a closed-shell reference spans exactly the singlet excitations, so triplets never
      enter the Jacobian. The gamma_2 closure modes are filtered post hoc by the spin symmetry of
      their contracted delta-gamma_1: singlet modes have the alpha-alpha and beta-beta blocks
      equal, M_S = 0 triplets have them opposite.
  (2) HUNGARIAN ASSIGNMENT. One-to-one matching by scipy linear_sum_assignment on the cosine
      matrix, replacing six independent argmaxes, which previously assigned one mode to two
      different FCI levels.
  (3) ACCEPTANCE TEST. tdhf(claude) must reproduce TDHF(pyscf) singlets. It is the only row whose
      correct answer is independently known, so no closure number is reported unless it passes.

Fingerprint: response of several random Hermitian one-body observables. For a mode, contract
delta-gamma_2 to delta-gamma_1 and take sum_pq O[p,q] delta-gamma_1[p,q]; for FCI state k,
<0|O_hat|k>. Cosine is phase-insensitive, |<a,b>|/(|a||b|).

Geometry: asymmetric trapezoid plus 0.05 * normal jitter. RHF reference (UHF is identical here).
Energies in Hartree.
"""
import sys
import numpy as np
import torch
from scipy.optimize import linear_sum_assignment
sys.path.insert(0, "."); sys.path.insert(0, "tools")
from pyscf import gto, scf, tdscf
from molecular_2rdm import ORB_B, molecule, separable_gpu
from qdiss import gpu
from qdiss.eom import fullci, gwt
from qdiss.eom.evaluate import assemble as np_assemble

CUM1 = gwt.terms_cumulant(1)
from qdiss.eom import wick
ONE2, TWO2 = wick.terms(2)
dev = gpu.device()
BASE = np.array([[0.00, 0.00, 0.0], [0.80, 0.00, 0.0],
                 [1.05, 0.72, 0.0], [-0.15, 0.72, 0.0]])
geom = BASE + 0.05 * np.random.default_rng(20260805).normal(size=BASE.shape)
atom = "; ".join("H %.6f %.6f %.6f" % tuple(r) for r in geom)
mol, mf, M, h, V, enuc = molecule(atom, "sto-3g", 0, 0)
N = mol.nelectron
n = M // 2

sz0 = [d for d in range(1 << M)
       if bin(d & 0x5555555555).count("1") == N // 2
       and bin(d & 0xAAAAAAAAAA).count("1") == N // 2]
idx = {d: i for i, d in enumerate(sz0)}
nd = len(sz0)


def op_matrix(apply_fn):
    Om = np.zeros((nd, nd), dtype=complex)
    b = np.zeros(1 << M, dtype=complex)
    for j, d in enumerate(sz0):
        b[:] = 0; b[d] = 1.0
        col = apply_fn(b)
        for d2, i in idx.items():
            Om[i, j] = col[d2]
    return Om


def Hfn(b):
    col = np.zeros_like(b)
    for p in range(M):
        for q in range(M):
            if abs(h[p, q]) > 1e-14:
                col += h[p, q] * fullci.ad(p, fullci.aa(q, b))
    for p, t, q, u in np.argwhere(np.abs(V) > 1e-14):
        col += 0.25 * V[p, t, q, u] * fullci.ad(p, fullci.ad(t, fullci.aa(u, fullci.aa(q, b))))
    return col


def S2fn(b):
    sp = np.zeros_like(b); sm = np.zeros_like(b)
    for p in range(n):
        sp += fullci.ad(2 * p, fullci.aa(2 * p + 1, b))
    out = np.zeros_like(b)
    for p in range(n):
        out += fullci.ad(2 * p + 1, fullci.aa(2 * p, sp))
    return out


Hs = op_matrix(Hfn); Hs = 0.5 * (Hs + Hs.conj().T)
S2s = op_matrix(S2fn); S2s = 0.5 * (S2s + S2s.conj().T)
es, vec = np.linalg.eigh(Hs)
s2 = np.real(np.einsum("ik,ij,jk->k", vec.conj(), S2s, vec))
sing = np.where(np.abs(s2) < 1e-6)[0]
print("S_z=0 sector dim %d; singlet states %d; E0 = %.9f" % (nd, len(sing), es[sing[0]] + enuc))
exc_all = (es - es[sing[0]]).real
sing_exc = np.array([exc_all[k] for k in sing[1:]])
keep = np.argsort(sing_exc)[:6]
targets = sing_exc[keep]
tstate = [sing[1:][k] for k in keep]
print("FCI singlet excitations (Ha): %s" % np.array2string(targets, precision=4))

NOBS = 6
ro = np.random.default_rng(11)
OBS = []
for _ in range(NOBS):
    A = ro.normal(size=(n, n)) + 1j * ro.normal(size=(n, n))
    a = 0.5 * (A + A.conj().T)
    O = np.zeros((M, M), dtype=complex)
    for p in range(n):
        for q in range(n):
            O[2 * p, 2 * q] = a[p, q]; O[2 * p + 1, 2 * q + 1] = a[p, q]
    OBS.append(O)
Osec = [op_matrix(lambda b, O=O: sum(O[p, q] * fullci.ad(p, fullci.aa(q, b))
                                    for p in range(M) for q in range(M)
                                    if abs(O[p, q]) > 1e-14)) for O in OBS]
# Full transition density matrices T_k[p,q] = <0|a+_p a_q|k>, all M^2 components, which is
# strictly more discriminating than projecting onto a handful of random observables.
def embed(c):
    f = np.zeros(1 << M, dtype=complex); f[sz0] = c; return f


def ladder(psi, k):
    import itertools as _it
    rows = []
    for ix in _it.product(*([range(M)] * k)):
        v = psi
        for pp in ix:
            v = fullci.aa(pp, v)
        rows.append(v)
    return np.array(rows)


p0 = embed(vec[:, sing[0]])
L01, L02 = ladder(p0, 1), ladder(p0, 2)
fci_T1, fci_T2 = [], []
for k in tstate:
    pk = embed(vec[:, k])
    fci_T1.append((L01.conj() @ ladder(pk, 1).T).reshape(M, M))
    fci_T2.append((L02.conj() @ ladder(pk, 2).T).reshape(M, M, M, M))
fci_T1 = np.array(fci_T1); fci_T2 = np.array(fci_T2)

td = tdscf.TDHF(mf); td.singlet = True; td.nstates = 10; td.verbose = 0; td.kernel()
tdp = np.sort(np.asarray(td.e).real)


def assign_bio(J, basisB, targetsH, label):
    """Biorthogonal projection. J is real and non-normal, so the component of a target along
    mode i is l_i^dag x with l_i^dag the i-th ROW of R^{-1}, not r_i^dag x. Score is the norm of
    the projection of the target onto mode i, |c_i| ||r_i|| / ||x||. Hungarian on the scores.
    """
    w, R = np.linalg.eig(J)
    Rinv = np.linalg.inv(R)
    keep = np.where(w.imag > 1e-6)[0]
    rn = np.array([np.linalg.norm(R[:, i]) for i in range(R.shape[1])])
    nb = len(basisB)
    bn = np.array([float(np.real(np.tensordot(b.conj(), b, axes=b.ndim))) for b in basisB])
    S = np.zeros((len(targets), len(keep)))
    for a in range(len(targets)):
        H = targetsH[a]
        x = np.array([float(np.real(np.tensordot(basisB[m].conj(), H, axes=H.ndim))) / bn[m]
                      for m in range(nb)])
        nx = np.linalg.norm(x)
        c = Rinv @ x
        for b, i in enumerate(keep):
            S[a, b] = abs(c[i]) * rn[i] / (nx + 1e-30)
    r, cc = linear_sum_assignment(-S)
    row = np.full(len(targets), np.nan); sc = np.zeros(len(targets))
    for a, b in zip(r, cc):
        row[a] = w[keep[b]].imag; sc[a] = S[a, b]
    err = np.nanmean(np.abs(row - targets))
    print("%-13s | %s  mean|err| %.4f  (%d modes)"
          % (label, "  ".join("%8.4f" % x for x in row), err, len(keep)))
    print("%-13s | %s  <- biorthogonal projection weight"
          % ("", "  ".join("%8.3f" % x for x in sc)))
    return err


EYE = np.eye(M, dtype=complex)
Z4 = np.zeros((M,) * 4, dtype=complex)


def d1_lam2zero(x):
    return np_assemble(CUM1[0], CUM1[1], h, V,
                       {("D", 1): EYE, ("G", 1): x, ("L", 2): Z4}, "pq", M)


def fock(g):
    return h + np.einsum("prqs,rs->pq", V, g)


g = np.zeros((M, M), dtype=complex)
for _ in range(800):
    _, vv = np.linalg.eigh(0.5 * (fock(g).conj() + fock(g).conj().conj().T))
    nw = vv[:, :N] @ vv[:, :N].conj().T
    if np.abs(nw - g).max() < 1e-13:
        g = nw; break
    g = 0.6 * nw + 0.4 * g
f0 = d1_lam2zero(g)
B1 = []
for i in range(n):
    b = np.zeros((n, n), dtype=complex); b[i, i] = 1.0; B1.append(b)
for i in range(n):
    for j in range(i + 1, n):
        b = np.zeros((n, n), dtype=complex); b[i, j] = 1.0; b[j, i] = 1.0; B1.append(b)
        b = np.zeros((n, n), dtype=complex); b[i, j] = 1j; b[j, i] = -1j; B1.append(b)
Bsa = []
for b in B1:
    E = np.zeros((M, M), dtype=complex)
    for p in range(n):
        for q in range(n):
            E[2 * p, 2 * q] = b[p, q]; E[2 * p + 1, 2 * q + 1] = b[p, q]
    Bsa.append(E)
n1 = len(Bsa); J1 = np.zeros((n1, n1)); eps = 1e-6
for k in range(n1):
    dp = (d1_lam2zero(g + eps * Bsa[k]) - f0) / eps
    for m in range(n1):
        J1[m, k] = float(np.real(np.trace(Bsa[m].conj().T @ dp))) / float(
            np.real(np.trace(Bsa[m].conj().T @ Bsa[m])))
H1 = [0.5 * (T + T.conj().T) for T in fci_T1]

print()
print("%-13s | %s" % ("method", "  ".join("%8.4f" % x for x in targets)))
print("%-13s | %s" % ("FCI singlet", "  ".join("%8s" % "--" for _ in targets)))
print("%-13s | %s  (reference)" % ("TDHF(pyscf)",
                                   "  ".join("%8.4f" % x for x in tdp[:len(targets)])))
assign_bio(J1, Bsa, H1, "tdhf(claude)")
w1c, _ = np.linalg.eig(J1)
om1 = np.sort(w1c.imag[w1c.imag > 1e-6])
nn = min(len(om1), len(tdp))
print("   mode set: %s" % np.array2string(om1, precision=4))

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
n2 = len(B2); Bg = [gpu.up(b, dev) for b in B2]
nrm = [float(np.real(np.tensordot(b.conj(), b, axes=4))) for b in B2]
psi0 = np.zeros(1 << M, dtype=complex); psi0[sz0] = vec[:, sing[0]]
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


for nm in ("lam3_0", "muller", "muller_raw"):
    f0b = rhs(gpu.up(g2_0, dev), nm)
    J = np.zeros((n2, n2))
    for k in range(n2):
        dn = gpu.down((rhs(gpu.up(g2_0, dev) + eps * Bg[k], nm) - f0b) / eps)
        for m in range(n2):
            J[m, k] = float(np.real(np.tensordot(B2[m].conj(), dn, axes=4))) / nrm[m]
    H2 = [0.5 * (T + np.transpose(T.conj(), (2, 3, 0, 1))) for T in fci_T2]
    assign_bio(J, B2, H2, nm)
