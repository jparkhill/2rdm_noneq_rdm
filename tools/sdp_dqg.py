#!/usr/bin/env python3
r"""How much of the 0.597 hartree does DQG recover? The same ADMM with Q and G added.

CONTEXT. tools/sdp_reference.py found that the two conditions the ladder construction makes
structural -- D >= 0 and contraction consistency -- bound the H4 energy at -2.540306 against
E_FCI = -1.943433, i.e. 1207% of E_corr. That killed the plan to fix the ladder variationally. The
question this file answers is how much is recoverable by adding the other two standard
N-representability cones, since CAR already hands the ladder the hole condition for free.

THE PROGRAM. gamma_1 is fixed to the exact FCI value; the variable is Mat_2, a C(M,2)-square
Hermitian matrix with Mat_2[(p,v),(q,w)] = 2 gamma_2[p,v,q,w] over sorted pairs.

    minimise  E = sum h gamma_1 + sum_{a,b} V[a,b] Mat_2[a,b] / 2 + E_nuc
    subject to  Tr_2 gamma_2 = (N-1) gamma_1                            (affine)
                D = Mat_2                             >= 0
                Q = two-hole matrix, antisym basis     >= 0             (affine in Mat_2)
                G = particle-hole matrix, M^2 square   >= 0             (affine in Mat_2)

Q and G come from qdiss.representability.dqg, whose build_Q is documented as checked against
explicit Fock space at 1.2e-16 and whose antisym_matrix uses the same factor-2 pair convention as
Mat_2. Both are affine in Mat_2 once gamma_1 is fixed.

ALGORITHM. ADMM with one splitting variable per cone. The x-step minimises
<C,x> + (rho/2) sum_c |M_c x + m_c - Z_c + U_c|^2 subject to A x = b, which is the KKT system
[[rho M^T M, A^T], [A, 0]] factored once; the Z-step projects each block onto the PSD cone. Hermitian
matrices are carried as real vectors through an explicit isometry (diagonal, then sqrt(2) times the
real and imaginary upper triangles) so no basis stack is stored -- the G block alone would otherwise
need 4096 matrices of size 64 x 64.

VALIDATION BEFORE USE: the exact gamma_2 must be feasible for all three cones and give E_FCI.

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
from qdiss.representability import dqg

CASES = (("H 0.00 0.00 0.0; H 0.80 0.00 0.0; H 1.05 0.72 0.0; H -0.15 0.72 0.0",
          "sto-3g", 4, 0, "H4 trapezoid"),)


def iso(n):
    """Real-vector <-> n x n Hermitian isometry (Hilbert-Schmidt preserving)."""
    iu, ju = np.triu_indices(n, 1)
    s = np.sqrt(2.0)

    def h2v(H):
        return np.concatenate([np.real(np.diag(H)), s * np.real(H[iu, ju]), s * np.imag(H[iu, ju])])

    def v2h(v):
        H = np.zeros((n, n), dtype=complex)
        H[np.arange(n), np.arange(n)] = v[:n]
        nt = len(iu)
        up = (v[n:n + nt] + 1j * v[n + nt:]) / s
        H[iu, ju] = up
        H[ju, iu] = up.conj()
        return H
    return h2v, v2h, n * n


def run(atom, basis, N, spin, label):
    mol, mf, M, h, V, enuc = molecule(atom, basis, 0, spin)
    dets = sector(M, N)
    es, vec = np.linalg.eigh(build_H_sector(h, V, M, dets))
    psi = np.zeros(1 << M, dtype=complex)
    psi[dets] = vec[:, 0]
    g1, g2 = fullci.rdm(psi, M, 1), fullci.rdm(psi, M, 2)
    e_fci, e_hf = float(fci.FCI(mf).kernel()[0]), float(mf.e_tot)
    PAIR = az.tuples(M, 2)
    npair, pos = len(PAIR), {t: i for i, t in enumerate(PAIR)}
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

    def to_g2(Mat2):
        return 0.5 * SA[:, :, None, None] * SA[None, None, :, :] * Mat2[AI][:, :, AI]

    def contract(Mat2):
        return np.einsum("pxqx->pq", to_g2(Mat2))

    hD, vD, nD = iso(npair)
    hG, vG, nG = iso(M * M)
    hA, vA, _ = iso(M)

    def cones(Mat2):
        d2 = to_g2(Mat2)
        Qm, _ = dqg.antisym_matrix(dqg.build_Q(d2, g1))
        Gm = dqg.as_matrix(dqg.build_G(d2, g1))
        return Mat2, Qm, Gm

    # affine maps M_c x + m_c, built column by column
    x_ex = hD(Mat2ex)
    z0 = cones(np.zeros((npair, npair), dtype=complex))
    m_c = [hD(z0[0]), hD(z0[1]), hG(z0[2])]
    Mcols = [[], [], []]
    Acols = []
    for k in range(nD):
        e = np.zeros(nD)
        e[k] = 1.0
        Mk = vD(e)
        cz = cones(Mk)
        Mcols[0].append(hD(cz[0]) - m_c[0])
        Mcols[1].append(hD(cz[1]) - m_c[1])
        Mcols[2].append(hG(cz[2]) - m_c[2])
        Acols.append(hA(contract(Mk)))
    Mm = [np.array(c).T for c in Mcols]
    A = np.array(Acols).T
    b = hA((N - 1.0) * g1)
    C = np.array([float(np.real(np.einsum("ab,ab->", Vp, vD(np.eye(nD)[k])))) / 2.0
                  for k in range(nD)])

    print("=== %s  M=%d N=%d | E_HF %.6f  E_FCI %.6f  E_corr %.6f ==="
          % (label, M, N, e_hf, e_fci, e_fci - e_hf))
    ce = cones(Mat2ex)
    print("    exact gamma_2 feasibility: |A x - b| %.2e | min eig D %.2e  Q %.2e  G %.2e | E dev %.2e"
          % (np.abs(A @ x_ex - b).max(),
             np.linalg.eigvalsh(0.5 * (ce[0] + ce[0].conj().T)).min(),
             np.linalg.eigvalsh(0.5 * (ce[1] + ce[1].conj().T)).min(),
             np.linalg.eigvalsh(0.5 * (ce[2] + ce[2].conj().T)).min(),
             abs(e1 + C @ x_ex - e_fci)))

    for tag, use in (("D only", [0]), ("D+Q", [0, 1]), ("D+Q+G", [0, 1, 2])):
        Ms = np.vstack([Mm[i] for i in use])
        ms = np.concatenate([m_c[i] for i in use])
        sizes = [Mm[i].shape[0] for i in use]
        rho = 1.0
        K = np.block([[rho * (Ms.T @ Ms), A.T], [A, np.zeros((A.shape[0],) * 2)]])
        Kf = np.linalg.inv(K)
        z = np.zeros(Ms.shape[0])
        u = np.zeros_like(z)
        for it in range(30000):
            rhs = np.concatenate([rho * (Ms.T @ (z - u - ms)) - C, b])
            x = (Kf @ rhs)[:nD]
            w = Ms @ x + ms + u
            zn, off = [], 0
            for i, n in zip(use, sizes):
                blk = w[off:off + n]
                Hm = (vD(blk) if i in (0, 1) else vG(blk))
                ev, W = np.linalg.eigh(0.5 * (Hm + Hm.conj().T))
                Z = (W * np.clip(ev, 0.0, None)) @ W.conj().T
                zn.append(hD(Z) if i in (0, 1) else hG(Z))
                off += n
            zn = np.concatenate(zn)
            u = u + (Ms @ x + ms) - zn
            rp = np.linalg.norm(Ms @ x + ms - zn)
            rd = rho * np.linalg.norm(zn - z)
            z = zn
            if rp < 1e-8 and rd < 1e-8:
                break
        ev = e1 + C @ x
        cz = cones(vD(x))
        print("    %-7s E = %.6f  E-E_FCI %+.6f  (%.0f%% of E_corr) | %5d its, rp %.1e rd %.1e"
              " | min eig D %.1e Q %.1e G %.1e"
              % (tag, ev, ev - e_fci, 100.0 * (ev - e_hf) / (e_fci - e_hf), it + 1, rp, rd,
                 np.linalg.eigvalsh(0.5 * (cz[0] + cz[0].conj().T)).min(),
                 np.linalg.eigvalsh(0.5 * (cz[1] + cz[1].conj().T)).min(),
                 np.linalg.eigvalsh(0.5 * (cz[2] + cz[2].conj().T)).min()))


for a, bs, n, s, l in (CASES if __name__ == "__main__" else ()):
    run(a, bs, n, s, l)
