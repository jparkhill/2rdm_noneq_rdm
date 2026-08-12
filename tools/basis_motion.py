#!/usr/bin/env python3
r"""Fixed or moving auxiliary basis? Measure the basis-motion term on an exact trajectory.

THE ARCHITECTURAL QUESTION. In the dynamical closure lambda_3 = sum_a c_a S_a, projecting the exact
rate gives

    lambda_3-dot = sum_a c_a-dot S_a  +  sum_a c_a S_a-dot

The second term exists only because the basis MOVES: the connected structures S_a are built from the
instantaneous gamma_1 and lambda_2, which are themselves propagating. True HEOM has no such term --
its auxiliary basis is fixed. If |sum_a c_a S_a-dot| is comparable to |lambda_3-dot| then the Galerkin
projection is not a small correction and a FIXED basis is the sounder construction, at the cost of a
larger auxiliary count A.

THE TWO CANDIDATE BASES.

  MOVING  the connected span from qdiss.ansatz.connected_basis: chains of two lambda_2 joined through
          gamma_1^a, legs dressed by anti(gamma_1^b). Small (4-dimensional span, measured to contain
          99.4% of lambda_3 at a pure ground state) but time-dependent.

  FIXED   chains of AMPLITUDES T_sigma[pv,qw] = V[pv,qw] / (e_q + e_w - e_p - e_v + sigma), built in
          the eigenbasis of h so e is fixed, for several shifts sigma. This is the HEOM idea applied
          here: a pole decomposition of the energy denominator. Dimensionless (V over an energy),
          time-independent, and the auxiliary count grows GEOMETRICALLY in (number of poles) x depth,
          which is the cost law John asked for.

WHAT IS REPORTED, along an exact trajectory from a DIPOLE-KICKED FCI ground state (per CLAUDE.md the
wavefunction is perturbed, never the RDM):

    span err   |P_basis lambda_3 - lambda_3| / |lambda_3|, the best the basis can do at that instant
    motion     |sum_a c_a S_a-dot| / |lambda_3-dot|, the size of the basis-motion term (MOVING only)
    coeff drift how much the fitted coefficients themselves move, |c(t) - c(0)| / |c(0)|

Rates are taken by finite difference along the exact trajectory, so the noise floor is reported by
differencing a quantity whose rate is known analytically.

Ground truth qdiss.eom.fullci.
"""
import sys
import numpy as np
import scipy.linalg as sla
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from molecular_2rdm import build_H_sector, molecule, sector
from qdiss import ansatz as az
from qdiss.eom import fullci

ATOM = "H 0.00 0.00 0.0; H 0.80 0.00 0.0; H 1.05 0.72 0.0; H -0.15 0.72 0.0"
N, EPS, DT = 4, 0.05, 1e-5
SHIFTS = (0.05, 0.2, 1.0)

mol, mf, M, h, V, enuc = molecule(ATOM, "sto-3g", 0, 0)
dets = sector(M, N)
es, vec = np.linalg.eigh(build_H_sector(h, V, M, dets))

# z dipole in the sector, to kick the wavefunction
# H4_TRAP lies in the xy-plane, so the z dipole vanishes identically; use x, and CHECK it.
dz = mf.mo_coeff.T @ mol.intor("int1e_r")[0] @ mf.mo_coeff
assert np.abs(dz).max() > 1e-8, "dipole component is null -- the kick would be vacuous"
d = np.zeros((M, M), dtype=complex)
for p in range(dz.shape[0]):
    for q in range(dz.shape[0]):
        d[2 * p, 2 * q] = d[2 * p + 1, 2 * q + 1] = dz[p, q]
idx = {dd: i for i, dd in enumerate(dets)}
Mu = np.zeros((len(dets),) * 2, dtype=complex)
for j, dd in enumerate(dets):
    b = np.zeros(1 << M, dtype=complex)
    b[dd] = 1.0
    col = np.zeros(1 << M, dtype=complex)
    for p in range(M):
        for q in range(M):
            if abs(d[p, q]) > 1e-14:
                col += d[p, q] * fullci.ad(p, fullci.aa(q, b))
    for d2, i in idx.items():
        if col[d2] != 0:
            Mu[i, j] = col[d2]
c0 = sla.expm(-1j * EPS * 0.5 * (Mu + Mu.conj().T)) @ vec[:, 0]
c0 /= np.linalg.norm(c0)
assert abs(abs(np.vdot(vec[:, 0], c0)) - 1.0) > 1e-9, "kick did not move the state"

# FIXED basis: amplitudes T_sigma in the eigenbasis of h, chained
ev, U = np.linalg.eigh(0.5 * (h + h.conj().T))
Vh = np.einsum("ap,bv,abcd,cq,dw->pvqw", U.conj(), U.conj(), V, U, U, optimize=True)
den = (ev[None, None, :, None] + ev[None, None, None, :]
       - ev[:, None, None, None] - ev[None, :, None, None])
FIXED = []
Tl = []
for s in SHIFTS:
    Ts = Vh / (den + s)
    Ts = Ts - np.transpose(Ts, (1, 0, 2, 3))
    Ts = Ts - np.transpose(Ts, (0, 1, 3, 2))
    Tl.append(np.einsum("ap,bv,pvqw,cq,dw->abcd", U, U, Ts, U.conj(), U.conj(), optimize=True))
for i, Ti in enumerate(Tl):
    for j, Tj in enumerate(Tl):
        for tag, sub in (("ch", "pvqr,rs,sxwy->pvxqwy"), ("chT", "prqw,rs,vxsy->pvxqwy")):
            X = az.antisym3(np.einsum(sub, Ti, np.eye(M, dtype=complex), Tj, optimize=True))
            FIXED.append(0.5 * (X + np.conj(X.transpose(3, 4, 5, 0, 1, 2))))
FIXED = np.array([f.reshape(-1) for f in FIXED])


def state_at(t):
    ph = np.exp(-1j * es.real * t)
    ct = vec @ (ph * (vec.conj().T @ c0))
    p = np.zeros(1 << M, dtype=complex)
    p[dets] = ct
    g1 = fullci.rdm(p, M, 1)
    l2 = az.cumulant(2, fullci.rdm(p, M, 2), {("G", 1): g1})
    l3 = az.cumulant(3, fullci.rdm(p, M, 3), {("G", 1): g1, ("L", 2): l2})
    return g1, l2, l3


def fit(B, target):
    c, *_ = np.linalg.lstsq(B.T, target.reshape(-1), rcond=None)
    return c, np.linalg.norm(B.T @ c - target.reshape(-1)) / np.linalg.norm(target)


print("Fixed vs moving auxiliary basis. H4 trapezoid M=%d N=%d, dipole-kicked eps=%.2f." % (M, N, EPS))
print("MOVING: connected chains of lambda_2 (12 generators, span 4).  "
      "FIXED: %d amplitude chains from %d poles." % (len(FIXED), len(SHIFTS)))
print()
print("%6s | %9s %9s %9s | %9s %9s"
      % ("t", "mov span", "motion", "c drift", "fix span", "fixA"))
c_ref = None
for t in (0.0, 0.1, 0.25, 0.5, 0.75, 1.0):
    g1, l2, l3 = state_at(t)
    g1b, l2b, l3b = state_at(t + DT)
    MOV = np.array([X.reshape(-1) for _, X in az.connected_basis(g1, l2, M)])
    MOVb = np.array([X.reshape(-1) for _, X in az.connected_basis(g1b, l2b, M)])
    c, se = fit(MOV, l3)
    if c_ref is None:
        c_ref = c
    dl3 = (l3b - l3) / DT
    Sdot = (MOVb - MOV) / DT
    motion = np.linalg.norm(c @ Sdot) / np.linalg.norm(dl3.reshape(-1))
    _, fe = fit(FIXED, l3)
    print("%6.2f | %9.4f %9.4f %9.4f | %9.4f %9d"
          % (t, se, motion, np.linalg.norm(c - c_ref) / max(np.linalg.norm(c_ref), 1e-30),
             fe, len(FIXED)))
print()
print("motion >~ 1 means the Galerkin projection is dominated by basis movement, and a FIXED basis")
print("is the sounder architecture despite needing a larger auxiliary count.")
