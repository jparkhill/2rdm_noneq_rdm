#!/usr/bin/env python3
r"""Is the rank-2 flow a COMMUTATOR flow?  Where does the spectral motion of Mat_2 live?

QUESTION (route F). At rank 1 the mean-field flow is a commutator, gamma_1-dot = i[f^T, gamma_1],
so the spectrum of gamma_1 is conserved for ANY energy functional and 0 <= gamma_1 <= 1 is
structural. The rank-2 analogue would be a closure for which

    Xdot = i [K, X] ,      X = Mat_2 ,  K Hermitian on pair space,

so that the spectrum of Mat_2 is conserved and gamma_2 >= 0 is structural. The exact flow is
Xdot = i(S^dag X - X S) with S = K - i Gamma/2 in general NON-Hermitian (XR-POSITIVITY.md), and the
part of the flow that is NOT a commutator is exactly the SPECTRAL motion: in the eigenbasis of X

    Xdot_ij = i K_ij (x_j - x_i) - (1/2) Gamma_ij (x_i + x_j) ,

so the block-diagonal of Xdot over the eigenspaces of X (x_i = x_j) cannot come from any K. Hence:

    a Hermitian-S closure exists  <=>  blockdiag_eig(Xdot) = 0,
    and |blockdiag_eig(Xdot)| / |Xdot| measures how far a flow is from EVERY such closure.

STATES. Near-equilibrium (dipole-kicked ground state, pure and thermal) and far-from-equilibrium
(the Hartree-Fock determinant at t = 0, where all cumulants vanish and V scatters into them; and an
interaction quench, ground state of h + 0.3 V propagated exactly under h + V, sampled at t = 0, 1, 2;
also the thermal version of that quench). RDMs of a mixed state are ensemble averages, exact.

FOR EACH STATE AND FLOW: spectral fraction (block-diagonal with degeneracy tolerance DEG_TOL),
the eigenvalue-rate range, the rate at near-kernel eigenvalues (positivity-violation rate; exact must
be ~0), and max |x_i-dot| / x_i over interior eigenvalues (the minimal diagonal dissipation).
Flows: exact; the gamma_2-carrying DRIFT (gamma_3 := 0 in the term list); the collision part
Phi[gamma_3] = total - drift with gamma_3 exact, separable (lambda_3 = 0), Mueller-2405 +/- Lackner.
Rank 1 control: exact gamma_1-dot and TDHF (must be spectral-free).

Symbols: Mat_2[(p,v),(q,w)] = 2 gamma_2[p,v,q,w] over sorted pairs p<v, q<w, gamma_2[pvqw] =
<a+_p a+_v a_w a_q>. Ground truth qdiss.eom.fullci; term list qdiss.eom.wick.
"""
import sys
import numpy as np
import scipy.linalg as sla
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from molecular_2rdm import build_H_sector, molecule, sector
from qdiss import ansatz as az
from qdiss import muller2405
from qdiss.eom import fullci, wick

H3 = ("H 0.00 0.00 0.0; H 1.00 0.00 0.0; H 0.35 0.90 0.0", 3, 1)
H4 = ("H 0.00 0.00 0.0; H 0.80 0.00 0.0; H 1.05 0.72 0.0; H -0.15 0.72 0.0", 4, 0)
KICK = 0.3
QUENCH_LAMBDA = 0.3
KER_TOL = 1e-8
INT_TOL = 1e-3
DEG_TOL = 1e-6


def sector_operator(op1, M, dets):
    idx = {d: i for i, d in enumerate(dets)}
    O = np.zeros((len(dets),) * 2, dtype=complex)
    for j, d in enumerate(dets):
        b = np.zeros(1 << M, dtype=complex)
        b[d] = 1.0
        col = np.zeros(1 << M, dtype=complex)
        for p in range(M):
            for q in range(M):
                if abs(op1[p, q]) > 1e-14:
                    col += op1[p, q] * fullci.ad(p, fullci.aa(q, b))
        for d2, i in idx.items():
            if col[d2] != 0:
                O[i, j] = col[d2]
    return 0.5 * (O + O.conj().T)


def dipole(mol, mf, M):
    dx = mf.mo_coeff.T @ mol.intor("int1e_r")[0] @ mf.mo_coeff
    assert np.abs(dx).max() > 1e-8, "null dipole component"
    d = np.zeros((M, M), dtype=complex)
    for p in range(dx.shape[0]):
        for q in range(dx.shape[0]):
            d[2 * p, 2 * q] = d[2 * p + 1, 2 * q + 1] = dx[p, q]
    return d


def pairs(M):
    return [(p, v) for p in range(M) for v in range(p + 1, M)]


def to_mat2(T, M):
    P = pairs(M)
    X = np.zeros((len(P), len(P)), dtype=complex)
    for I, (p, v) in enumerate(P):
        for J, (q, w) in enumerate(P):
            X[I, J] = 2.0 * T[p, v, q, w]
    return X


def ensemble_rdms(dm, dets, M, H):
    """Exact RDMs (k=1,2,3) and rates (k=1,2) of a sector density matrix dm, by eigen-ensemble."""
    w, U = np.linalg.eigh(0.5 * (dm + dm.conj().T))
    out = [0] * 5
    for i in np.where(w > 1e-12)[0]:
        psi = np.zeros(1 << M, dtype=complex)
        psi[dets] = U[:, i]
        vals = (fullci.rdm(psi, M, 1), fullci.rdm(psi, M, 2), fullci.rdm(psi, M, 3),
                fullci.rdm_dot(psi, H, M, 1), fullci.rdm_dot(psi, H, M, 2))
        out = [o + w[i] * v for o, v in zip(out, vals)]
    return out, int((w > 1e-3).sum())


def analyse(X, Xd, label, boundary=(0.0,)):
    X = 0.5 * (X + X.conj().T)
    Xd = 0.5 * (Xd + Xd.conj().T)
    x, U = np.linalg.eigh(X)
    D = U.conj().T @ Xd @ U
    # block-diagonal over degenerate eigenspaces
    same = np.abs(x[:, None] - x[None, :]) < DEG_TOL
    B = np.where(same, D, 0.0)
    xdot = np.real(np.diag(D))
    nrm = np.linalg.norm(Xd)
    spec = np.linalg.norm(B) / nrm if nrm > 0 else 0.0
    dist = np.min([np.abs(x - b) for b in boundary], axis=0)
    ker = dist < KER_TOL
    interior = dist > INT_TOL
    ker_rate = np.abs(xdot[ker]).max() if ker.any() else 0.0
    gam = np.abs(xdot[interior] / dist[interior]).max() if interior.any() else 0.0
    # entropy rate of the spectrum: fermionic -x ln x - (1-x) ln(1-x) at rank 1, -x ln x at rank 2
    xs = np.clip(x, 1e-300, None)
    if len(boundary) == 2:
        dS = -np.sum(xdot * (np.log(xs) - np.log(np.clip(1 - x, 1e-300, None))))
    else:
        dS = -np.sum(xdot * (np.log(xs) + 1.0))
    print("    %-36s |Xd|=%.4f  spectral=%.4f  n_ker=%2d  |xdot|ker=%.1e  "
          "max|xdot/dist|=%.3f  xdot in [%+.1e,%+.1e]  dS/dt=%+.2e"
          % (label, nrm, spec, ker.sum(), ker_rate, gam, xdot.min(), xdot.max(), dS))
    return spec


def run(label, dm, dets, M, N, h, V):
    H = fullci.build_H(h, V, M)
    (g1, g2, g3, d1, d2), rk = ensemble_rdms(dm, dets, M, H)
    l2 = az.cumulant(2, g2, {("G", 1): g1})
    n1 = np.linalg.eigvalsh(g1)
    EYE = np.eye(M, dtype=complex)
    print("=== %s   M=%d N=%d  rank(rho)~%d  occupations [%.4f, %.4f]  |lambda_2|=%.4f ==="
          % (label, M, N, rk, n1.min(), n1.max(), np.linalg.norm(l2)))
    print("  RANK 1  (X = gamma_1, boundaries 0 and 1)")
    analyse(g1, d1, "exact", boundary=(0.0, 1.0))
    tab1 = {("D", 1): EYE, ("G", 1): g1, ("G", 2): az.wedge(g1, g1)}
    analyse(g1, wick.gamma_dot(1, h, V, tab1, M), "TDHF (lambda_2 = 0)", boundary=(0.0, 1.0))
    tab1[("G", 2)] = g2
    chk1 = np.abs(wick.gamma_dot(1, h, V, tab1, M) - d1).max()

    print("  RANK 2  (X = Mat_2, boundary 0)")
    X = to_mat2(g2, M)
    x = np.linalg.eigvalsh(0.5 * (X + X.conj().T))
    print("    Mat_2 dim %d, eigenvalues in [%.1e, %.4f]" % (X.shape[0], x.min(), x.max()))
    Xd_exact = to_mat2(d2, M)
    tab2 = {("D", 1): EYE, ("G", 1): g1, ("G", 2): g2, ("G", 3): np.zeros((M,) * 6, dtype=complex)}
    Xd_drift = to_mat2(wick.gamma_dot(2, h, V, tab2, M), M)
    tab2[("G", 3)] = g3
    chk2 = np.abs(wick.gamma_dot(2, h, V, tab2, M) - d2).max()
    print("    [check] wick vs fullci: rank 1 %.1e, rank 2 %.1e" % (chk1, chk2))
    analyse(X, Xd_exact, "exact total")
    analyse(X, Xd_drift, "drift (gamma_3 := 0 in term list)")
    analyse(X, Xd_exact - Xd_drift, "Phi[gamma_3 exact]")
    recon = {
        "lambda_3 = 0": az.separable(3, {("G", 1): g1, ("L", 2): l2}),
        "muller2405 + Lackner": muller2405.gamma3(g1, g2, N, True),
        "muller2405 raw": muller2405.gamma3(g1, g2, N, False),
    }
    for name, g3r in recon.items():
        tab2[("G", 3)] = g3r
        Xd_r = to_mat2(wick.gamma_dot(2, h, V, tab2, M), M)
        analyse(X, Xd_r - Xd_drift, "Phi[%s]" % name)
        print("    %-36s |Xd_closure - Xd_exact|/|Xd_exact| = %.4f"
              % ("", np.linalg.norm(Xd_r - Xd_exact) / np.linalg.norm(Xd_exact)))
    print()


def propagate(dm, es, vec, t):
    ph = np.exp(-1j * es.real * t)
    return vec @ (ph[:, None] * (vec.conj().T @ dm @ vec) * ph.conj()[None, :]) @ vec.conj().T


for atom, N, spin in (H3, H4):
    mol, mf, M, h, V, enuc = molecule(atom, "sto-3g", 0, spin)
    dets = sector(M, N)
    tag = "H3 bent" if N == 3 else "H4 trapezoid"
    es, vec = np.linalg.eigh(build_H_sector(h, V, M, dets))
    Uk = sla.expm(-1j * KICK * sector_operator(dipole(mol, mf, M), M, dets))

    # near-equilibrium: kicked ground state (pure), kicked thermal state
    dm = np.outer(Uk @ vec[:, 0], (Uk @ vec[:, 0]).conj())
    run(tag + " | kicked ground state", dm, dets, M, N, h, V)
    if N == 4:
        w = np.exp(-(es - es[0]) / 0.5)
        w /= w.sum()
        dm = Uk @ (vec * w) @ vec.conj().T @ Uk.conj().T
        run(tag + " | kicked thermal kT=0.5", dm, dets, M, N, h, V)

    # far from equilibrium: HF determinant at t = 0 (all cumulants zero)
    hf_det = sum(1 << p for p in range(N))     # MO basis: lowest N spin orbitals
    i0 = dets.index(hf_det)
    dm = np.zeros((len(dets),) * 2, dtype=complex)
    dm[i0, i0] = 1.0
    run(tag + " | HF determinant, t=0", dm, dets, M, N, h, V)

    # interaction quench: ground state of h + lam V, evolved under h + V
    esA, vecA = np.linalg.eigh(build_H_sector(h, QUENCH_LAMBDA * V, M, dets))
    dm0 = np.outer(vecA[:, 0], vecA[:, 0].conj())
    for t in (0.0, 1.0, 2.0):
        run(tag + " | V-quench %.1f->1, pure, t=%.0f" % (QUENCH_LAMBDA, t),
            propagate(dm0, es, vec, t), dets, M, N, h, V)
    if N == 4:
        w = np.exp(-(esA - esA[0]) / 0.5)
        w /= w.sum()
        dmT = (vecA * w) @ vecA.conj().T
        run(tag + " | V-quench, thermal kT=0.5, t=1", propagate(dmT, es, vec, 1.0), dets, M, N, h, V)
