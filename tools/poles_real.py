#!/usr/bin/env python3
"""Pole structure of cumulant-reconstruction 1-RDM EOMs on real 2- and 3-electron systems.

Not RDMFT: no energy functional is minimized and nothing is inserted as a coupling matrix.
The object under test is a reconstruction of the CUMULANT lambda_2[gamma_1], substituted into
the unapproximated MKNO 1-RDM EOM (qdiss.eom.gwt cumulant form, 2 + 12 terms in gamma_1 and
lambda_2, validated against qdiss.eom.fullci at 1e-15). The RDMFT functionals below supply the
SHAPE of the cumulant only.

Matrix elements from PySCF, spin-orbital, in the convention of qdiss.eom.hamiltonian, each
Hamiltonian checked against the PySCF FCI total energy.

Reference point. At a determinantal gamma_1 every occupation-diagonal reconstruction vanishes
identically, and the sqrt(n) ones have divergent derivatives for occupied-virtual pairs, so
linear response there is either trivial or singular. Poles are therefore taken at the exact
FCI ground-state gamma_1, which has fractional occupations and is common to every method, so
the comparison is like-for-like. TDHF at its own HF stationary point is reported as well,
since that is the textbook TDHF pole.

Poles are eigenvalues of the Jacobian of gamma_1-dot over a Hermitian basis. Where the
reference is not a fixed point of a given EOM there is an inhomogeneous term, which shifts no
eigenvalue; that caveat is stated rather than hidden.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
from pyscf import ao2mo, gto, scf
from qdiss.eom import fullci, gwt
from qdiss.eom.evaluate import assemble

CUM1 = gwt.terms_cumulant(1)
NAMES = ("hf", "blocked", "muller", "power06", "bbc3ish", "pnofish", "mp2")


def build(atom, basis, charge=0, spin=0):
    mol = gto.M(atom=atom, basis=basis, charge=charge, spin=spin, verbose=0)
    mf = (scf.ROHF(mol) if spin else scf.RHF(mol)).run()
    C = mf.mo_coeff
    nmo = C.shape[1]
    hc = C.T @ mf.get_hcore() @ C
    eri = ao2mo.restore(1, ao2mo.kernel(mol, C), nmo)
    M = 2 * nmo
    h = np.zeros((M, M), dtype=complex)
    for p in range(nmo):
        for q in range(nmo):
            h[2 * p, 2 * q] = hc[p, q]
            h[2 * p + 1, 2 * q + 1] = hc[p, q]
    phys = np.zeros((M,) * 4)
    for p in range(nmo):
        for q in range(nmo):
            for r in range(nmo):
                for s in range(nmo):
                    v = eri[p, r, q, s]
                    for sp in (0, 1):
                        for sq in (0, 1):
                            phys[2 * p + sp, 2 * q + sq, 2 * r + sp, 2 * s + sq] += v
    V = (phys - np.transpose(phys, (0, 1, 3, 2))).astype(complex)
    return M, h, V, float(mol.energy_nuc())


def wedge(a, b):
    return np.einsum("pq,vw->pvqw", a, b) - np.einsum("pw,vq->pvqw", a, b)


def from_c(c, U, M):
    D = np.zeros((M,) * 4, dtype=complex)
    for p in range(M):
        for v in range(M):
            D[p, v, v, p] += 0.5 * c[p, v]
            D[p, v, p, v] -= 0.5 * c[p, v]
    return np.einsum("ap,bv,pvqw,cq,dw->abcd", U, U, D, U.conj(), U.conj(), optimize=True)


def c_matrix(name, n):
    M = len(n)
    prod = np.outer(n, n)
    root = np.sqrt(np.clip(prod, 0, None))
    if name == "muller":
        return prod - root
    if name == "power06":
        return prod - np.clip(prod, 0, None) ** 0.6
    if name == "bbc3ish":
        X = root.copy()
        weak = np.outer(n < 0.5, n < 0.5)
        X = np.where(weak & ~np.eye(M, dtype=bool), -X, X)
        np.fill_diagonal(X, np.diag(prod))
        return prod - X
    if name == "pnofish":
        order = np.argsort(-n)
        partner = np.empty(M, dtype=int)
        for k in range(0, M - 1, 2):
            partner[order[k]] = order[k + 1]
            partner[order[k + 1]] = order[k]
        if M % 2:
            partner[order[-1]] = order[-1]
        X = prod.copy()
        for p in range(M):
            X[p, partner[p]] = -root[p, partner[p]]
        return prod - X
    if name == "blocked":
        a = n - n * n
        T = a.sum()
        den = T - 2 * a
        if np.abs(den).min() < 1e-12:
            return np.zeros((M, M))
        Q = (a / den).sum()
        f = -a / ((1.0 + Q) * den)
        return -2.0 * (np.outer(a, f) + np.outer(f, a))
    raise ValueError(name)


def lam2(g1, name, h, V, M):
    if name == "hf":
        return np.zeros((M,) * 4, dtype=complex)
    n, U = np.linalg.eigh(0.5 * (g1 + g1.conj().T))
    n = np.clip(n, 1e-10, 1.0)
    if name != "mp2":
        return from_c(c_matrix(name, n), U, M)
    F = h + np.einsum("prqs,rs->pq", V, g1)
    e = np.real(np.diag(U.conj().T @ F @ U))
    Vn = np.einsum("ap,bv,abcd,cq,dw->pvqw", U.conj(), U.conj(), V, U, U, optimize=True)
    occ = np.einsum("q,w,p,v->pvqw", n, n, 1 - n, 1 - n)
    den = (e[None, None, :, None] + e[None, None, None, :]
           - e[:, None, None, None] - e[None, :, None, None])
    den = np.where(np.abs(den) < 1e-6, np.inf, den)
    T = Vn * occ / den
    T = T - np.transpose(T, (1, 0, 2, 3))
    T = T - np.transpose(T, (0, 1, 3, 2))
    L = 0.25 * (T + np.transpose(T.conj(), (2, 3, 0, 1)))
    return np.einsum("ap,bv,pvqw,cq,dw->abcd", U, U, L, U.conj(), U.conj(), optimize=True)


def gdot(g1, name, h, V, M):
    t = {("D", 1): np.eye(M, dtype=complex), ("G", 1): g1,
         ("L", 2): lam2(g1, name, h, V, M)}
    return assemble(CUM1[0], CUM1[1], h, V, t, "pq", M)


def herm_basis(M):
    B = []
    for i in range(M):
        E = np.zeros((M, M), dtype=complex); E[i, i] = 1.0; B.append(E)
    for i in range(M):
        for j in range(i + 1, M):
            E = np.zeros((M, M), dtype=complex); E[i, j] = 1.0; E[j, i] = 1.0; B.append(E)
            E = np.zeros((M, M), dtype=complex); E[i, j] = 1j; E[j, i] = -1j; B.append(E)
    return B


def poles(g1, name, h, V, M, eps=1e-5):
    B = herm_basis(M)
    n = len(B)
    J = np.zeros((n, n))
    f0 = gdot(g1, name, h, V, M)
    for k, Bk in enumerate(B):
        dp = (gdot(g1 + eps * Bk, name, h, V, M) - f0) / eps
        for m, Bm in enumerate(B):
            J[m, k] = float(np.real(np.trace(Bm.conj().T @ dp))) / float(
                np.real(np.trace(Bm.conj().T @ Bm)))
    ev = np.linalg.eigvals(J)
    return np.sort(np.abs(ev.imag[ev.imag > 1e-6]))


SYS = (("H 0 0 0; H 0 0 0.74", "6-31g", 2, 0, 0, "H2 / 6-31g   (2e)"),
       ("H 0 0 0; H 0 0 1.6; H 0 0 3.2", "sto-3g", 3, 0, 1, "H3 linear / sto-3g   (3e)"))

for atom, basis, N, chg, spin, label in SYS:
    M, h, V, enuc = build(atom, basis, chg, spin)
    Hm = fullci.build_H(h, V, M)
    dim = 1 << M
    popc = np.array([bin(s).count("1") for s in range(dim)])
    sub = np.where(popc == N)[0]
    es, vec = np.linalg.eigh(Hm[np.ix_(sub, sub)])
    psi0 = np.zeros(dim, dtype=complex)
    psi0[sub] = vec[:, 0]
    g1_fci = fullci.rdm(psi0, M, 1)
    nfci = np.sort(np.linalg.eigvalsh(0.5 * (g1_fci + g1_fci.conj().T)))[::-1]
    exc = np.sort((es - es[0]).real)[1:7]
    print("=== %s   M=%d N=%d   E_FCI %.8f ===" % (label, M, N, es[0].real + enuc))
    print("   FCI occupations : %s" % np.array2string(nfci[:5], precision=4))
    print("   exact poles     : %s" % np.array2string(exc, precision=4))

    def fock(g):
        return h + np.einsum("prqs,rs->pq", V, g)

    g = np.zeros((M, M), dtype=complex)
    for _ in range(400):
        _, vv = np.linalg.eigh(0.5 * (fock(g).conj() + fock(g).conj().conj().T))
        new = vv[:, :N] @ vv[:, :N].conj().T
        if np.abs(new - g).max() < 1e-12:
            g = new
            break
        g = 0.6 * new + 0.4 * g
    print("   TDHF at HF ref  : %s   (|g1dot| %.1e)"
          % (np.array2string(poles(g, "hf", h, V, M)[:6], precision=4),
             np.abs(gdot(g, "hf", h, V, M)).max()))
    print()
    print("   poles at the FCI reference gamma_1")
    for name in NAMES:
        l2 = lam2(g1_fci, name, h, V, M)
        tid = np.abs(np.einsum("pxqx->pq", l2)
                     + (g1_fci - g1_fci @ g1_fci)).max()
        om = poles(g1_fci, name, h, V, M)
        print("   %-9s |l2| %7.4f trace %8.1e | %s"
              % (name, np.abs(l2).max(), tid, np.array2string(om[:6], precision=4)))
    print()
