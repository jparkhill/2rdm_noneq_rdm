#!/usr/bin/env python3
"""Trace-corrected cumulant reconstructions: can pole quality and the trace identity coexist?

Not RDMFT. The object is a reconstruction of the cumulant lambda_2[gamma_1] inserted into the
unapproximated MKNO 1-RDM EOM (qdiss.eom.gwt, validated against qdiss.eom.fullci at 1e-15).
Matrix elements from PySCF, checked against FCI.

The tension found previously: the reconstruction satisfying the partial-trace identity exactly
(blocked) has mediocre poles, while the ones with good poles (Mueller family, mp2) violate it
by 0.1 to 0.5.

The fix is closed-form. For any antisymmetric W with partial trace S = sum_x W[px,qx],

    Corr(S) = (1/(M-2)) (S^I + I^S) - [tr S/((M-2)(M-1))] (I^I)

has partial trace exactly S, using sum_x (a^b)[px,qx] = a tr(b) - a b. So W - Corr(S) is
trace-free, and adding the blocked particular solution, whose partial trace is -A with
A = gamma_1 - gamma_1^2, gives

    lambda_2_tc = W - Corr(S_W) + lambda_2_blocked

which satisfies the identity exactly while keeping every coherence of W. Only the trace
channel is replaced, and it is replaced by the blocked structure rather than smeared over the
identity.

Reported for each reconstruction, raw and trace-corrected: the trace residual, the norm, and
the poles at the FCI reference against exact excitation energies.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from poles_real import (NAMES, build, c_matrix, from_c, gdot, herm_basis, lam2, poles, wedge)
from qdiss.eom import fullci
from qdiss.eom.evaluate import assemble
from qdiss.eom import gwt

CUM1 = gwt.terms_cumulant(1)


def sym_wedge(a, b):
    return wedge(a, b) + wedge(b, a)


def corr(S, M):
    I = np.eye(M, dtype=complex)
    return ((1.0 / (M - 2)) * sym_wedge(S, I)
            - (np.trace(S) / ((M - 2) * (M - 1))) * wedge(I, I))


def blocked_l2(g1, M):
    n, U = np.linalg.eigh(0.5 * (g1 + g1.conj().T))
    n = np.clip(n, 1e-10, 1.0)
    return from_c(c_matrix("blocked", n), U, M)


def lam2_tc(g1, name, h, V, M):
    W = lam2(g1, name, h, V, M)
    S = np.einsum("pxqx->pq", W)
    return W - corr(S, M) + blocked_l2(g1, M)


def gdot_tc(g1, name, h, V, M):
    t = {("D", 1): np.eye(M, dtype=complex), ("G", 1): g1,
         ("L", 2): lam2_tc(g1, name, h, V, M)}
    return assemble(CUM1[0], CUM1[1], h, V, t, "pq", M)


def poles_tc(g1, name, h, V, M, eps=1e-5):
    B = herm_basis(M)
    n = len(B)
    J = np.zeros((n, n))
    f0 = gdot_tc(g1, name, h, V, M)
    for k, Bk in enumerate(B):
        dp = (gdot_tc(g1 + eps * Bk, name, h, V, M) - f0) / eps
        for m, Bm in enumerate(B):
            J[m, k] = float(np.real(np.trace(Bm.conj().T @ dp))) / float(
                np.real(np.trace(Bm.conj().T @ Bm)))
    ev = np.linalg.eigvals(J)
    return np.sort(np.abs(ev.imag[ev.imag > 1e-6]))


SYS = (("H 0 0 0; H 0 0 0.74", "6-31g", 2, 0, 0, "H2 / 6-31g   (2e)"),
       ("H 0 0 0; H 0 0 1.6; H 0 0 3.2", "sto-3g", 3, 0, 1, "H3 linear / sto-3g   (3e)"))
TEST = ("blocked", "muller", "bbc3ish", "mp2")

for atom, basis, N, chg, spin, label in SYS:
    M, h, V, enuc = build(atom, basis, chg, spin)
    Hm = fullci.build_H(h, V, M)
    popc = np.array([bin(s).count("1") for s in range(1 << M)])
    sub = np.where(popc == N)[0]
    es, vec = np.linalg.eigh(Hm[np.ix_(sub, sub)])
    psi0 = np.zeros(1 << M, dtype=complex)
    psi0[sub] = vec[:, 0]
    g1 = fullci.rdm(psi0, M, 1)
    A = g1 - g1 @ g1
    exc = np.sort((es - es[0]).real)[1:7]
    print("=== %s ===" % label)
    print("   exact poles : %s" % np.array2string(exc, precision=4))
    print("   TDHF (l2=0) : %s" % np.array2string(poles(g1, "hf", h, V, M)[:6], precision=4))
    print()
    print("   %-9s %-4s | %8s %7s | %s" % ("recon", "tc", "trace", "|l2|", "poles"))
    for name in TEST:
        for tc in (False, True):
            L = lam2_tc(g1, name, h, V, M) if tc else lam2(g1, name, h, V, M)
            tid = np.abs(np.einsum("pxqx->pq", L) + A).max()
            om = poles_tc(g1, name, h, V, M) if tc else poles(g1, name, h, V, M)
            print("   %-9s %-4s | %8.1e %7.4f | %s"
                  % (name, "yes" if tc else "no", tid, np.abs(L).max(),
                     np.array2string(om[:6], precision=4)))
    print()
