#!/usr/bin/env python3
r"""Convex reference: the D + contraction variational 2-RDM bound, by ADMM.

WHY. tools/ladder_closure.py tried to optimise the ladder directly and failed -- constraint
residuals 1e-1 to 3e-1, and an oracle error that INCREASED with chi', which is impossible for a
nested family. That indicts the optimiser, not the ansatz, so a trustworthy reference point is
needed before any more gradient descent.

THE SIMPLIFICATION. The convex reference does not need the ladder at all. EVERY ladder, at every
chi', produces a Mat_2 that is positive semidefinite and contraction-consistent with gamma_1. So

    minimise  E   subject to   Mat_2 >= 0 ,   Tr_2 gamma_2 = (N-1) gamma_1

is a LOWER BOUND on what any chi' can achieve, and it is a convex semidefinite program on a
C(M,2) x C(M,2) matrix rather than an (M r_1)-square one. This is the classical D-condition
variational 2-RDM bound with contraction consistency, and its looseness is a known and expected
feature -- the point here is to measure it for OUR systems so the ladder results can be judged
against something real.

SYMBOLS. Mat_2[(p,v),(q,w)] = 2 gamma_2[p,v,q,w] over sorted pairs p<v, q<w; the antisymmetric
extension is gamma_2[p,v,q,w] = sign * Mat_2[sort(p,v), sort(q,w)] / 2, zero when p = v or q = w.
Tr_2 gamma_2[p,q] = sum_x gamma_2[p,x,q,x]. Energy in the convention of qdiss.eom.hamiltonian,
E = sum_pq h gamma_1 + sum_{a,b sorted} V[a,b] Mat_2[a,b] / 2 + E_nuc, checked against the exact
gamma_2 before use.

ALGORITHM. ADMM on  min <C,X> s.t. A(X) = b, X >= 0, in a real orthonormal basis of Hermitian
matrices. X-step is the affine projection of (Z - U - C/rho), computed once from a cached
(A A^T)^{-1}; Z-step is the eigenvalue projection onto the PSD cone; U accumulates. Convergence is
reported as primal and dual residuals so the reader can judge, and the run is checked by confirming
that the EXACT gamma_2 is feasible (residual ~ 0) and that its energy is E_FCI.

Ground truth qdiss.eom.fullci and the PySCF FCI energy.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from molecular_2rdm import build_H_sector, molecule, sector
from pyscf import fci
from qdiss import ansatz as az
from qdiss.eom import fullci

CASES = (("H 0.00 0.00 0.0; H 0.80 0.00 0.0; H 1.05 0.72 0.0; H -0.15 0.72 0.0",
          "sto-3g", 4, 0, "H4 trapezoid"),
         ("Li 0 0 0; H 0 0 1.6", "sto-3g", 4, 0, "LiH"))


def herm_basis(n):
    """Orthonormal real basis of n x n Hermitian matrices (Hilbert-Schmidt)."""
    B = []
    for i in range(n):
        E = np.zeros((n, n), dtype=complex)
        E[i, i] = 1.0
        B.append(E)
    s = 1.0 / np.sqrt(2.0)
    for i in range(n):
        for j in range(i + 1, n):
            E = np.zeros((n, n), dtype=complex)
            E[i, j] = s
            E[j, i] = s
            B.append(E)
            E = np.zeros((n, n), dtype=complex)
            E[i, j] = 1j * s
            E[j, i] = -1j * s
            B.append(E)
    return B


def run(atom, basis, N, spin, label):
    mol, mf, M, h, V, enuc = molecule(atom, basis, 0, spin)
    dets = sector(M, N)
    es, vec = np.linalg.eigh(build_H_sector(h, V, M, dets))
    psi = np.zeros(1 << M, dtype=complex)
    psi[dets] = vec[:, 0]
    g1, g2 = fullci.rdm(psi, M, 1), fullci.rdm(psi, M, 2)
    e_fci, e_hf = float(fci.FCI(mf).kernel()[0]), float(mf.e_tot)
    PAIR = az.tuples(M, 2)
    np_ = len(PAIR)
    pos = {t: i for i, t in enumerate(PAIR)}
    Vp = np.array([[V[p, v, q, w] for (q, w) in PAIR] for (p, v) in PAIR])
    Mat2ex = az.tensor_to_coord(g2, M, 2) * 2.0
    e1 = float(np.real(np.einsum("pq,pq->", h, g1))) + enuc

    AI = np.zeros((M, M), dtype=int)
    SA = np.zeros((M, M))
    for p in range(M):
        for v in range(M):
            if p != v:
                AI[p, v] = pos[(min(p, v), max(p, v))]
                SA[p, v] = 1.0 if p < v else -1.0

    def contract(Mat2):
        g2t = 0.5 * SA[:, :, None, None] * SA[None, None, :, :] * Mat2[AI][:, :, AI]
        return np.einsum("pxqx->pq", g2t)

    # objective and constraint map in a real Hermitian basis
    B = herm_basis(np_)
    C = np.array([float(np.real(np.einsum("ab,ab->", Vp, Bk))) / 2.0 for Bk in B])
    HB = herm_basis(M)
    A = np.zeros((len(HB), len(B)))
    for k, Bk in enumerate(B):
        c = contract(Bk)
        A[:, k] = [float(np.real(np.einsum("pq,pq->", Hk.conj(), c))) for Hk in HB]
    tgt = (N - 1.0) * g1
    b = np.array([float(np.real(np.einsum("pq,pq->", Hk.conj(), tgt))) for Hk in HB])

    x_ex = np.array([float(np.real(np.einsum("ab,ab->", Bk.conj(), Mat2ex))) for Bk in B])
    print("=== %s  M=%d N=%d | E_HF %.6f  E_FCI %.6f  E_corr %.6f ==="
          % (label, M, N, e_hf, e_fci, e_fci - e_hf))
    print("    feasibility of the EXACT gamma_2: |A x_ex - b| = %.2e,  E = %.9f (dev %.2e)"
          % (np.abs(A @ x_ex - b).max(), e1 + C @ x_ex, abs(e1 + C @ x_ex - e_fci)))

    AAT_inv = np.linalg.inv(A @ A.T)

    def proj_affine(w):
        return w - A.T @ (AAT_inv @ (A @ w - b))

    # stack the basis once; a Python sum over 784 arrays per iteration is the bottleneck
    Bs = np.array(B)

    def mat(x):
        return np.einsum("k,kab->ab", x.astype(complex), Bs)

    def vecof(Z):
        return np.real(np.einsum("kab,ab->k", Bs.conj(), Z))

    rho = 1.0
    z = np.zeros(len(B))
    u = np.zeros_like(z)
    for it in range(20000):
        x = proj_affine(z - u - C / rho)
        Xm = mat(x + u)
        w, W = np.linalg.eigh(0.5 * (Xm + Xm.conj().T))
        Zm = (W * np.clip(w, 0.0, None)) @ W.conj().T
        z_new = vecof(Zm)
        u = u + x - z_new
        rp, rd = np.linalg.norm(x - z_new), rho * np.linalg.norm(z_new - z)
        z = z_new
        if rp < 1e-9 and rd < 1e-9:
            break
    Xf = mat(z)
    e_sdp = e1 + C @ z
    print("    ADMM: %d iters, primal res %.2e, dual res %.2e, min eig Mat_2 %.2e, "
          "|contraction res| %.2e" % (it + 1, rp, rd, np.linalg.eigvalsh(0.5 * (Xf + Xf.conj().T)).min(),
                                      np.abs(contract(Xf) - tgt).max()))
    print("    E(D + contraction SDP) = %.9f    E_SDP - E_FCI = %+.6f    = %.0f%% of E_corr"
          % (e_sdp, e_sdp - e_fci, 100.0 * (e_sdp - e_hf) / (e_fci - e_hf)))
    print("    -> any ladder, at any chi', is bounded below by this number.")
    print()


for a, bs, n, s, l in (CASES if __name__ == "__main__" else ()):
    run(a, bs, n, s, l)
