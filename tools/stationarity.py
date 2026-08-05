#!/usr/bin/env python3
"""How stationary is the FCI density under each cumulant reconstruction?

The reference used for every pole table is the exact FCI ground-state gamma_1. Because
gamma_1-dot depends only on gamma_1 and gamma_2, and the FCI ground state is an eigenstate of
H, the EXACT gamma_1-dot vanishes there. So |gamma_1-dot| under a reconstruction measures how
far that reconstruction is from making the true density stationary, and how much weight the
linearised poles can carry.

Also reported: the true TDHF stationary point, where |gamma_1-dot| is zero by construction, so
its poles are textbook response poles while the others are natural frequencies of a linearised
flow at a non-fixed point.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from poles_real import build, gdot, lam2, wedge
from trace_corrected import gdot_tc, lam2_tc
from qdiss.eom import fullci
from qdiss.eom.evaluate import assemble
from qdiss.eom import gwt

CUM1 = gwt.terms_cumulant(1)
SYS = (("H 0 0 0; H 0 0 0.74", "6-31g", 2, 0, 0, "H2 / 6-31g (2e)"),
       ("H 0 0 0; H 0 0 1.6; H 0 0 3.2", "sto-3g", 3, 0, 1, "H3 linear / sto-3g (3e)"))

for atom, basis, N, chg, spin, label in SYS:
    M, h, V, enuc = build(atom, basis, chg, spin)
    Hm = fullci.build_H(h, V, M)
    popc = np.array([bin(s).count("1") for s in range(1 << M)])
    sub = np.where(popc == N)[0]
    es, vec = np.linalg.eigh(Hm[np.ix_(sub, sub)])
    psi0 = np.zeros(1 << M, dtype=complex)
    psi0[sub] = vec[:, 0]
    g1 = fullci.rdm(psi0, M, 1)
    l2ex = fullci.rdm(psi0, M, 2) - wedge(g1, g1)
    t = {("D", 1): np.eye(M, dtype=complex), ("G", 1): g1, ("L", 2): l2ex}
    dex = assemble(CUM1[0], CUM1[1], h, V, t, "pq", M)
    print("=== %s ===" % label)
    print("   exact lambda_2 at the FCI density : |gamma_1-dot| = %.3e  (must be zero)"
          % np.abs(dex).max())
    print("   %-10s | %12s %12s | %12s" % ("recon", "|g1dot| raw", "|g1dot| tc", "|l2-l2exact|"))
    for name in ("hf", "blocked", "muller", "bbc3ish", "mp2"):
        raw = np.abs(gdot(g1, name, h, V, M)).max()
        tc = np.abs(gdot_tc(g1, name, h, V, M)).max()
        err = np.abs(lam2(g1, name, h, V, M) - l2ex).max()
        print("   %-10s | %12.5f %12.5f | %12.5f" % (name, raw, tc, err))

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
    print("   TDHF at its own HF fixed point    : |gamma_1-dot| = %.2e"
          % np.abs(gdot(g, "hf", h, V, M)).max())
    print("   |gamma_1_HF - gamma_1_FCI| = %.5f" % np.abs(g - g1).max())
    print()
