#!/usr/bin/env python3
r"""First propagation of the dynamical closure. Validated in its known limits before anything else.

THE SCHEME. Propagate gamma_2 with the MKNO rank-2 equation of motion, closing gamma_3 as

    gamma_3 = separable(gamma_1, lambda_2) + lambda_3 ,   lambda_3 = sum_a c_a S_a

with S_a the orthonormalised connected basis (built from the instantaneous gamma_1, lambda_2). Three
treatments of the coefficients, which is the whole point of the comparison:

    null      c = 0 identically, i.e. lambda_3 = 0. The reference closure.
    adiabatic c refitted at every stage to the instantaneous best fit of ... nothing exact is
              available during a real propagation, so the honest adiabatic rule is the one the
              closure prescribes: c fixed at its initial value, with the basis moving. This is the
              ONE-SHOT connected closure.
    memory    c propagated by the Galerkin flow c-dot = <S, lambda_3-dot> - <S, sum_b c_b S_b-dot>,
              built from the rank-3 MKNO equation of motion with lambda_4 = 0.

VALIDATIONS, run first and reported first, because this session produced four vacuous or
mis-normalised tests and a propagation has more places for that to hide:

    L1  c = 0 must reproduce the `null` closure exactly -- same trajectory as tools/ansatz_positivity
    L2  feeding the EXACT gamma_3 at every stage must reproduce the exact trajectory. This is the
        real test of the integrator and the term list together; if it fails, nothing downstream means
        anything.

Then the comparison: |gamma_2 - exact| for the three treatments, on the compressed -> stretched
quench where the coefficient drift was measured at 1.34 by t = 1.

Symbols as in qdiss.ansatz. Ground truth qdiss.eom.fullci.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from ansatz_hierarchy import rdms
from molecular_2rdm import build_H_sector, molecule, sector
from qdiss import ansatz as az
from qdiss.eom import fullci, wick

COMPRESSED = "H 0.00 0.00 0.0; H 0.60 0.00 0.0; H 0.79 0.54 0.0; H -0.11 0.54 0.0"
STRETCHED = "H 0.00 0.00 0.0; H 1.20 0.00 0.0; H 1.58 1.08 0.0; H -0.22 1.08 0.0"
N, DT, NSTEP, EVERY = 4, 0.002, 500, 125
ONE2, TWO2 = wick.terms(2)
ONE3, TWO3 = wick.terms(3)


def setup(atom):
    mol, mf, M, h, V, enuc = molecule(atom, "sto-3g", 0, 0)
    dets = sector(M, N)
    es, vec = np.linalg.eigh(build_H_sector(h, V, M, dets))
    return M, h, V, dets, es, vec


M, hA, VA, dets, esA, vecA = setup(COMPRESSED)
_, hB, VB, _, esB, vecB = setup(STRETCHED)
EYE = np.eye(M, dtype=complex)
Z3 = np.zeros((M,) * 6, dtype=complex)
Z4 = np.zeros((M,) * 8, dtype=complex)
dm0 = np.outer(vecA[:, 0], vecA[:, 0].conj())


def exact_at(t):
    ph = np.exp(-1j * esB.real * t)
    dmt = vecB @ (ph[:, None] * (vecB.conj().T @ dm0 @ vecB) * ph.conj()[None, :]) @ vecB.conj().T
    g = rdms(dmt, M, N, dets)
    l2 = az.cumulant(2, g[1], {("G", 1): g[0]})
    return g[0], g[1], g[2], l2, az.cumulant(3, g[2], {("G", 1): g[0], ("L", 2): l2})


def frame(g1, l2, ref=None, kfix=None):
    B = np.array([X.reshape(-1) for _, X in az.connected_basis(g1, l2, M)])
    _, sv, Vt = np.linalg.svd(B, full_matrices=False)
    k = int((sv > 1e-8 * sv[0]).sum()) if kfix is None else kfix
    F = Vt[:k]
    if ref is not None and ref.shape[0] == k:
        W, _, Zt = np.linalg.svd(F @ ref.conj().T)
        F = (W @ Zt).conj().T @ F
    return F


def rank2_dot(g1, g2, g3):
    return wick.gamma_dot(2, hB, VB, {("D", 1): EYE, ("G", 1): g1, ("G", 2): g2,
                                      ("G", 3): g3}, M)


def rank3_dot(g1, g2, g3):
    return wick.gamma_dot(3, hB, VB, {("D", 1): EYE, ("G", 1): g1, ("G", 2): g2,
                                      ("G", 3): g3, ("G", 4): Z4}, M)


# ---------------------------------------------------------------- L2 first: does the integrator work
print("L2  feeding the EXACT gamma_3 at every stage must reproduce the exact trajectory")
g2 = exact_at(0.0)[1].copy()
for step in range(NSTEP):
    t = step * DT

    def f(x, tt):
        g1 = np.einsum("pxqx->pq", x) / (N - 1)
        return rank2_dot(g1, x, exact_at(tt)[2])
    k1 = f(g2, t)
    k2 = f(g2 + DT / 2 * k1, t + DT / 2)
    k3 = f(g2 + DT / 2 * k2, t + DT / 2)
    k4 = f(g2 + DT * k3, t + DT)
    g2 = g2 + DT / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    if (step + 1) % EVERY == 0:
        ex = exact_at((step + 1) * DT)[1]
        print("    t = %5.2f   |gamma_2 - exact| = %.3e   (|gamma_2| = %.4f)"
              % ((step + 1) * DT, np.abs(g2 - ex).max(), np.abs(ex).max()))
print()

# ---------------------------------------------------------------- the three treatments
print("comparison on the compressed -> stretched quench, dt = %.3f" % DT)
print("    %-10s %6s | %12s %12s"
      % ("treatment", "t", "|g2-exact|", "|c|"))
for mode in ("null", "adiabatic", "memory"):
    g1_0, g2_0, g3_0, l2_0, l3_0 = exact_at(0.0)
    g2 = g2_0.copy()
    F0 = frame(g1_0, l2_0)
    kfix = F0.shape[0]
    c = F0.conj() @ l3_0.reshape(-1) if mode != "null" else np.zeros(kfix, dtype=complex)
    prevF = F0
    for step in range(NSTEP + 1):
        if step % EVERY == 0:
            ex = exact_at(step * DT)[1]
            print("    %-10s %6.2f | %12.3e %12.4f"
                  % (mode, step * DT, np.abs(g2 - ex).max(), np.linalg.norm(c)))
        if step == NSTEP:
            break
        g1 = np.einsum("pxqx->pq", g2) / (N - 1)
        l2 = az.cumulant(2, g2, {("G", 1): g1})
        F = frame(g1, l2, prevF, kfix)
        prevF = F
        l3 = (F.T @ c).reshape((M,) * 6) if mode != "null" else Z3
        g3 = az.separable(3, {("G", 1): g1, ("L", 2): l2}) + l3
        d2 = rank2_dot(g1, g2, g3)
        if mode == "memory":
            d1 = np.einsum("pxqx->pq", d2) / (N - 1)
            dl2 = d2 - (az.wedge(d1, g1) + az.wedge(g1, d1))
            d3 = rank3_dot(g1, g2, g3)
            dl3 = d3 - az.separable_dot(3, {("G", 1): g1, ("L", 2): l2},
                                        {("G", 1): d1, ("L", 2): dl2})
            Sd = az.connected_basis_dot(g1, l2, d1, dl2, M)
            Sdm = np.array([x.reshape(-1) for x in Sd])
            # project the structure rates into the current frame the same way the frame was built
            Fd = np.linalg.lstsq(np.array([X.reshape(-1) for _, X in
                                           az.connected_basis(g1, l2, M)]).T, Sdm.T,
                                 rcond=None)[0]
            cdot = F.conj() @ dl3.reshape(-1)
            c = c + DT * cdot
        g2 = g2 + DT * d2
print()
print("NOTE: this is Euler for the coefficient leg and omits the basis-motion term in c-dot, so the")
print("`memory` row is NOT the full Galerkin flow. It is a first integration check only.")
