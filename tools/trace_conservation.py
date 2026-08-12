#!/usr/bin/env python3
r"""Is the trace identity a CONSERVED QUANTITY of the joint flow, and does it survive truncation?

THE REFRAMING. lambda_3 is not a functional of (gamma_1, lambda_2); it is an independent dynamical
object. So the single-trace identity

    Delta  =  Tr_3 lambda_3  -  anti(E, lambda_2)  =  0 ,     E = gamma_1 - I/2,
    Tr_3 T[pv,qw] = sum_y T[p v y, q w y],  anti(E,X) = E contracted into each of X's four slots

should be read not as a constraint on a reconstruction but as an INVARIANT of the joint flow of
(gamma_1, lambda_2, lambda_3). If the flow conserves it, then imposing it once at t = 0 and
propagating a projection that preserves a linear invariant keeps it exactly, and no reconstruction
ever has to satisfy it pointwise.

WHAT DECIDES WHETHER THAT WORKS. lambda_3-dot comes from the rank-3 equation of motion, which carries
gamma_4, i.e. lambda_4. Two possibilities:

  (i) d/dt Delta vanishes for ANY lambda_4 -- then truncating the auxiliary tower at lambda_4 = 0
      preserves the identity automatically and the conservation is free;
  (ii) d/dt Delta depends on lambda_4 -- then truncation BREAKS the conservation and the identity
      must be imposed on the projected coefficient flow as an explicit linear condition.

Either outcome is usable, but they demand different machinery, so this must be measured.

METHOD. At an exact FCI state, non-stationary (an eigenstate makes every rate vanish and the test
vacuous), take the exact rates gamma_k-dot from qdiss.eom.fullci.rdm_dot, convert them to cumulant
rates by differentiating the MKNO expansion,

    lambda_2-dot = gamma_2-dot - (gamma_1-dot wedge gamma_1 + gamma_1 wedge gamma_1-dot)
    lambda_3-dot = gamma_3-dot - d/dt separable_3(gamma_1, lambda_2)   [by finite difference of the
                                                                        separable part, to avoid
                                                                        hand-differentiating 16 terms]

and evaluate d/dt Delta three ways: with the exact gamma_4 (control, must be ~0 since the identity
holds at all times), and with the rank-3 equation of motion truncated at gamma_4 = 0 and at
lambda_4 = 0, which are the two candidate tower terminations.

Ground truth qdiss.eom.fullci. Term lists qdiss.eom.wick, validated against it at 8.7e-16.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from molecular_2rdm import build_H_sector, molecule, sector
from qdiss import ansatz as az
from qdiss.eom import fullci, wick
from qdiss.symbolic import mkno

UP4, LO4 = ("p", "v", "x", "z"), ("q", "w", "y", "t")
SEP4 = [(s, ",".join("".join(u) + "".join(l) for n, r, u, l in f) + "->pvxzqwyt",
         [(n, r) for n, r, u, l in f])
        for s, f in mkno.gamma_to_cumulants(UP4, LO4)
        if not any(n == "L" and r == 4 for n, r, u, l in f)]

CASES = (("H 0.00 0.00 0.0; H 1.00 0.00 0.0; H 0.35 0.90 0.0", "sto-3g", 3, 1, "H3 bent"),
         ("H 0.00 0.00 0.0; H 0.80 0.00 0.0; H 1.05 0.72 0.0; H -0.15 0.72 0.0",
          "sto-3g", 4, 0, "H4 trapezoid"))
DT = 1e-6


def cums(psi, M):
    g1 = fullci.rdm(psi, M, 1)
    g2 = fullci.rdm(psi, M, 2)
    g3 = fullci.rdm(psi, M, 3)
    l2 = az.cumulant(2, g2, {("G", 1): g1})
    l3 = az.cumulant(3, g3, {("G", 1): g1, ("L", 2): l2})
    return g1, g2, g3, l2, l3


def delta(g1, l2, l3):
    return np.einsum("pvyqwy->pvqw", l3) - az.trace_target(3, g1, l2)


for atom, basis, N, spin, label in CASES:
    mol, mf, M, h, V, enuc = molecule(atom, basis, 0, spin)
    dets = sector(M, N)
    es, vec = np.linalg.eigh(build_H_sector(h, V, M, dets))
    rng = np.random.default_rng(7)
    c = vec[:, 0] + 0.4 * (rng.normal(size=len(dets)) + 1j * rng.normal(size=len(dets)))
    c /= np.linalg.norm(c)
    psi = np.zeros(1 << M, dtype=complex)
    psi[dets] = c
    H = fullci.build_H(h, V, M)
    g1, g2, g3, l2, l3 = cums(psi, M)
    EYE = np.eye(M, dtype=complex)

    # exact rates
    d1 = fullci.rdm_dot(psi, H, M, 1)
    d2 = fullci.rdm_dot(psi, H, M, 2)
    d3 = fullci.rdm_dot(psi, H, M, 3)
    dl2 = d2 - (az.wedge(d1, g1) + az.wedge(g1, d1))
    # exact lambda_3-dot: finite difference of the separable part along the true trajectory
    ph = np.exp(-1j * es.real * DT)
    ct = vec @ (ph * (vec.conj().T @ c))
    pt = np.zeros(1 << M, dtype=complex)
    pt[dets] = ct
    g1b, g2b, g3b, l2b, l3b = cums(pt, M)
    dl3_fd = (l3b - l3) / DT

    print("=== %s   M=%d N=%d ===" % (label, M, N))
    print("    |Delta| at t=0 (must be ~0)                    = %.3e" % np.abs(delta(g1, l2, l3)).max())
    print("    |lambda_3-dot (fd) - consistency check| : |lambda_2-dot fd - analytic| = %.3e"
          % np.abs((l2b - l2) / DT - dl2).max())

    def ddelta(dl3):
        # d/dt [ Tr_3 lambda_3 - anti(E, lambda_2) ],  E-dot = gamma_1-dot
        return (np.einsum("pvyqwy->pvqw", dl3)
                - az.anti(d1, l2, 3) - az.anti(g1 - 0.5 * EYE, dl2, 3))

    print("    d/dt Delta with the EXACT rate                 = %.3e" % np.abs(ddelta(dl3_fd)).max())

    # the two candidate terminations, via the rank-3 EOM
    g4 = fullci.rdm(psi, M, 4)
    sep4 = None
    tab = {("D", 1): EYE, ("G", 1): g1, ("G", 2): g2, ("G", 3): g3}
    for tag, g4v in (("gamma_4 exact", g4),
                     ("gamma_4 = 0", np.zeros((M,) * 8, dtype=complex)),
                     ("lambda_4 = 0", None)):
        if g4v is None:
            t4 = {("G", 1): g1, ("L", 2): l2, ("L", 3): l3}
            sep4 = np.zeros((M,) * 8, dtype=complex)
            for s, sub, keys in SEP4:
                sep4 += s * np.einsum(sub, *[t4[k] for k in keys], optimize=True)
            g4v = sep4
        tab[("G", 4)] = g4v
        d3a = wick.gamma_dot(3, h, V, tab, M)
        # convert to a lambda_3 rate: subtract d/dt separable_3, taken from the exact trajectory
        dsep3 = (d3 - dl3_fd)          # since lambda_3-dot = gamma_3-dot - d/dt sep_3
        dl3a = d3a - dsep3
        print("    d/dt Delta with %-14s               = %.3e   (|d3| = %.4f)"
              % (tag, np.abs(ddelta(dl3a)).max(), np.abs(d3a).max()))
    print()
