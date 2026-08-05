#!/usr/bin/env python3
"""Poles, trace errors and positivity errors for each 2-RDM closure on molecular systems.

The three figures of merit are the ones that were asked for: pole positions against FCI and their
errors, trace errors in the 1- and 2-RDM, and positivity errors in the 1- and 2-RDM. No other
quantity is ranked.

APPROXIMATION, common to every row below. The equation propagated is the rank-2 Heisenberg EOM,
    gamma_2-dot[pv,qw] = +i <[H, a+_p a+_v a_w a_q]>
derived by plain vacuum Wick in qdiss.eom.wick, term-identical to the MKNO route qdiss.eom.gwt,
and validated against exact Fock-space linear algebra qdiss.eom.fullci at 1e-15. The PROPAGATED
STATE VARIABLE is gamma_2 alone; gamma_1 is not propagated but recovered at every step as
gamma_1 = Tr_2(gamma_2)/(N-1). This is the scheme of Lackner et al. and of Mueller et al. The only
approximation is the reconstruction of gamma_3 from gamma_1 and gamma_2:

    muller      Mueller Eq. 14 closed-system limit,
                (N-2)/N [g2 g1 + g2 g1 + g1 g2] - 2(N-1)(N-2)/N^2 g1 g1 g1, antisymmetrised,
                then made contraction-consistent by the Lackner Eq. 44/45 correction so that
                Tr_3(gamma_3) = (N-2) gamma_2 exactly.
    muller_raw  the same reconstruction WITHOUT the Lackner correction, so contraction
                consistency is violated. Included only to isolate what the correction does.
    lam3_0      ours: lambda_3 = 0, i.e. gamma_3 is the separable part of the MKNO cumulant
                expansion, built from gamma_1 and lambda_2. Contraction-inconsistent.
Moment negation (gamma_3 = 0) is deliberately absent. It is not a reconstruction but a deletion of
the 8 terms in which V reaches a third electron, leaving gamma_2-dot = -i[h(+)h+V, gamma_2] with a
gamma_2-independent generator, so its poles are differences of two-electron eigenvalues of H_2
(verified to 4.6e-05) and carry no N-body information; N never enters. It is excluded by the
project ground rules.

Integrals from PySCF in the convention of qdiss.eom.hamiltonian; each Hamiltonian is checked
against the PySCF FCI total energy before use. The exact reference is diagonalised in the
N-electron sector. Contractions run on GPU via qdiss.gpu, validated against the numpy path at
complex128 round-off in tools/gpu_validate.py.

POLES are obtained by response diagonalisation: the Jacobian of the gamma_2 flow is built by
finite differences over a Hermitian, pair-antisymmetric basis at the exact FCI ground-state
gamma_2, and its eigenvalue imaginary parts are the poles. The reference is the FCI excitation
spectrum of the same Hamiltonian in the same sector. Caveat stated rather than hidden: the FCI
ground state is a fixed point of the exact flow but not of any closure, so an inhomogeneous term
is present; it shifts no eigenvalue but means these are natural frequencies of the linearised
flow rather than textbook response poles.
"""
import sys
import time
import numpy as np
import torch
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from molecular_2rdm import (NAMES, ORB_B, SEP3, build_H_sector, molecule, sector, separable_gpu)
from qdiss import gpu
from qdiss.eom import fullci, wick

ONE2, TWO2 = wick.terms(2)
dev = gpu.device()

# Asymmetric trapezoidal H4, nearest-neighbour distances about 20 percent shorter than the
# 1.0 Angstrom linear chain, chosen to remove the threefold pole degeneracy that made the linear
# chain uninformative. Parallel sides 0.80 and 1.20; legs 0.762 and 0.735, so no mirror symmetry.
H4_TRAP = "H 0.00 0.00 0.0; H 0.80 0.00 0.0; H 1.05 0.72 0.0; H -0.15 0.72 0.0"
SYS = ((H4_TRAP, "sto-3g", 4, 0, 0, "H4 asymmetric trapezoid / sto-3g"),)

for atom, basis, N, chg, spin, label in SYS:
    mol, mf, M, h, V, enuc = molecule(atom, basis, chg, spin)
    dets = sector(M, N)
    Hs = build_H_sector(h, V, M, dets)
    es, vec = np.linalg.eigh(Hs)
    psi0 = np.zeros(1 << M, dtype=complex)
    psi0[dets] = vec[:, 0]
    g2_0 = fullci.rdm(psi0, M, 2)
    g1_0 = np.einsum("pxqx->pq", g2_0) / (N - 1)
    hg, Vg, Ig = gpu.up(h, dev), gpu.up(V, dev), gpu.up(np.eye(M), dev)
    PAIRS = [(p, v) for p in range(M) for v in range(p + 1, M)]
    NP = len(PAIRS)

    def recon(nm, g1, g2):
        if nm == "muller":
            return gpu.muller_gamma3(g1, g2, N, M, True, ORB_B)
        if nm == "muller_raw":
            return gpu.muller_gamma3(g1, g2, N, M, False, ORB_B)
        if nm == "lam3_0":
            return separable_gpu(g1, g2, M)
        raise ValueError("unknown reconstruction %r; moment negation is excluded" % nm)

    def rhs(x, nm):
        g1 = torch.einsum("pxqx->pq", x) / (N - 1)
        return gpu.assemble(ONE2, TWO2, hg, Vg,
                            {("D", 1): Ig, ("G", 1): g1, ("G", 2): x,
                             ("G", 3): recon(nm, g1, x)}, "pvqw")

    def pair_to_tensor(Z):
        T = np.zeros((M,) * 4, dtype=complex)
        for i, (p, v) in enumerate(PAIRS):
            for j, (q, w) in enumerate(PAIRS):
                x = 0.5 * Z[i, j]
                T[p, v, q, w] = x; T[v, p, q, w] = -x
                T[p, v, w, q] = -x; T[v, p, w, q] = x
        return T

    basis_t = []
    for a in range(NP):
        for b in range(a, NP):
            for ph in (1.0, 1j):
                Z = np.zeros((NP, NP), dtype=complex)
                Z[a, b] += ph; Z[b, a] += np.conj(ph)
                if np.abs(Z).max() < 1e-12:
                    continue
                basis_t.append(pair_to_tensor(Z))
    n = len(basis_t)
    Bg = [gpu.up(b, dev) for b in basis_t]
    norms = [float(np.real(np.tensordot(b.conj(), b, axes=4))) for b in basis_t]

    exc = np.sort((es - es[0]).real)
    exc = exc[exc > 1e-8][:8]
    print("=== %s   M=%d N=%d   sector dim %d ===" % (label, M, N, len(dets)))
    print("    ground E %.9f   FCI excitations %s"
          % (es[0].real + enuc, np.array2string(exc, precision=4)))
    print()
    print("TABLE 1  poles by response diagonalisation of the gamma_2 flow at the FCI ground state")
    print("         propagated variable gamma_2 alone, gamma_1 = Tr_2(gamma_2)/(N-1)")
    print("    %-11s | %s | %9s" % ("gamma_3 recon", "  ".join("%8.4f" % e for e in exc),
                                    "mean |err|"))
    print("    %-11s | %s |" % ("FCI", "  ".join("%8s" % "--" for e in exc)))
    eps = 1e-6
    for nm in NAMES:
        t0 = time.perf_counter()
        f0 = rhs(gpu.up(g2_0, dev), nm)
        J = np.zeros((n, n))
        for k in range(n):
            d = (rhs(gpu.up(g2_0, dev) + eps * Bg[k], nm) - f0) / eps
            dn = gpu.down(d)
            for m in range(n):
                J[m, k] = float(np.real(np.tensordot(basis_t[m].conj(), dn, axes=4))) / norms[m]
        w = np.linalg.eigvals(J)
        om = np.sort(np.abs(w.imag[w.imag > 1e-6]))
        matched = [om[np.argmin(np.abs(om - e))] if len(om) else np.nan for e in exc]
        err = np.mean([abs(m - e) for m, e in zip(matched, exc)])
        print("    %-11s | %s | %9.4f   (%.0fs)"
              % (nm, "  ".join("%8.4f" % m for m in matched), err,
                 time.perf_counter() - t0))
    print()

    print("TABLE 2  trace and positivity errors along propagation, same approximation per row")
    print("         tr gamma_1 must equal %d, tr gamma_2 must equal %d" % (N, N * (N - 1)))
    print("    %-11s | %5s | %10s %10s | %10s %10s %10s"
          % ("recon", "t", "|tr g1-N|", "|tr g2-NN|", "min n_i", "max n_i", "min eig D"))
    st = {nm: gpu.up(g2_0, dev) for nm in NAMES}

    def rk4(f, y, dt):
        k1 = f(y); k2 = f(y + dt / 2 * k1); k3 = f(y + dt / 2 * k2); k4 = f(y + dt * k3)
        return y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

    dt, ns, evy = 0.005, 200, 100
    rows = {nm: [] for nm in NAMES}
    for step in range(ns + 1):
        if step % evy == 0:
            for nm in NAMES:
                g2 = gpu.down(st[nm])
                g1 = np.einsum("pxqx->pq", g2) / (N - 1)
                occ = np.linalg.eigvalsh(0.5 * (g1 + g1.conj().T))
                A = np.array([[2.0 * g2[p, v, q, w] for (q, w) in PAIRS] for (p, v) in PAIRS])
                rows[nm].append((step * dt,
                                 abs(np.trace(g1).real - N),
                                 abs(np.einsum("pvpv->", g2).real - N * (N - 1)),
                                 occ.min(), occ.max(),
                                 np.linalg.eigvalsh(0.5 * (A + A.conj().T)).min()))
        if step < ns:
            for nm in NAMES:
                st[nm] = rk4(lambda y, k=nm: rhs(y, k), st[nm], dt)
    for nm in NAMES:
        for r in rows[nm]:
            print("    %-11s | %5.2f | %10.2e %10.2e | %10.6f %10.6f %10.6f"
                  % (nm, *r))
    print()
