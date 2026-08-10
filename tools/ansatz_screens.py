#!/usr/bin/env python3
"""Screens 1 and 2 of CLAUDE.md applied to the reconstruction hierarchy, at ranks 2 and 3.

SCREEN 2, EXTENSIVITY, is tested as pure algebra on the reconstruction map, with no chemistry, so
the answer is exact and free of any basis-matching question.

    Construction: take a system whose exact cumulants are known, here H3 bent / sto-3g at its FCI
    ground state with M_A = 6 spin orbitals, and EMBED it in a larger space of M spin orbitals by
    padding every tensor with zeros. The added orbitals are empty and non-interacting, so the
    exact lambda_k of the embedded system is exactly the padded lambda_k, and the exact partial
    trace target S is exactly the padded S. Any reconstruction that is extensive must therefore
    return the padded lambda_k unchanged.

    leak      max |lambda_k^rec| over index tuples that touch a spectator orbital, divided by
              max |lambda_k^rec| over index tuples wholly inside the fragment. Screen 2 demands
              EXACTLY zero cross-fragment elements, so any nonzero value is a failure.
    drift     max change of the within-fragment block relative to the M = 6 answer, divided by
              the magnitude of that block. Screen 2 demands the within-fragment part not change
              when a distant spectator is added.

SCREEN 1, LIPSCHITZ, is the derivative norm of the reconstruction with respect to its inputs:
    rank 2:  input gamma_1,           output lambda_2
    rank 3:  input (gamma_1, lambda_2), output lambda_3
sampled over 120 random Hermitian directions at step 1e-4, reported as
max ||d lambda_k||_F / ||d input||_F. Screen 1 rejects any closure whose value diverges as the
state approaches the physical boundary, which is where an unregularised gamma_1-dependent inverse
or a square root of an occupation blows up. Both a near-boundary state (the pure FCI ground state,
occupations reaching 0.0104 and 0.9881) and an interior state (thermal, kT = 0.5) are used,
because a finite constant at an interior point proves nothing.

Symbols as in qdiss.ansatz. Ground truth qdiss.eom.fullci in the N-electron sector.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from ansatz_hierarchy import fro, rank2_set, rank3_set, rdms
from molecular_2rdm import build_H_sector, molecule, sector
from qdiss import ansatz as az
from qdiss.eom import fullci

H3_BENT = "H 0.00 0.00 0.0; H 1.00 0.00 0.0; H 0.35 0.90 0.0"
N, MA = 3, 6


def pad(T, M):
    """Embed a tensor on MA orbitals into M orbitals by zero padding."""
    out = np.zeros((M,) * T.ndim, dtype=complex)
    out[tuple(slice(0, MA) for _ in range(T.ndim))] = T
    return out


def split_norms(T, M, k):
    """(max |T| inside the fragment, max |T| touching a spectator)."""
    A = np.abs(T)
    inside = A[tuple(slice(0, MA) for _ in range(2 * k))].max()
    B = A.copy()
    B[tuple(slice(0, MA) for _ in range(2 * k))] = 0.0
    return float(inside), float(B.max())


mol, mf, M0, h, V, enuc = molecule(H3_BENT, "sto-3g", 0, 1)
assert M0 == MA
dets = sector(MA, N)
es, vec = np.linalg.eigh(build_H_sector(h, V, MA, dets))
psi = np.zeros(1 << MA, dtype=complex)
psi[dets] = vec[:, 0]
g1A = fullci.rdm(psi, MA, 1)
g2A = fullci.rdm(psi, MA, 2)
g3A = fullci.rdm(psi, MA, 3)
l2A = az.cumulant(2, g2A, {("G", 1): g1A})
l3A = az.cumulant(3, g3A, {("G", 1): g1A, ("L", 2): l2A})

print("SCREEN 2  extensivity of the minimum-norm partial-trace reconstruction ('trace' level)")
print("          fragment H3 bent / sto-3g, MA = 6, embedded in M spin orbitals by zero padding")
print("          the exact lambda_k is the padded lambda_k, so leak and drift must both be 0")
print()
print("    %-6s %4s | %12s %12s | %12s %12s"
      % ("rank", "M", "|within|", "leak/within", "drift/within", "|exact blk|"))
for k, lam, low in ((2, l2A, g1A), (3, l3A, l2A)):
    ref = None
    for M in (6, 8, 10, 12):
        S = az.trace_target(k, pad(g1A, M), None if k == 2 else pad(l2A, M))
        # zero padding of gamma_1 makes the spectator orbitals empty, but E = gamma_1 - I/2
        # is -1/2 on them, so the rank-3 target must be rebuilt from the padded cumulants and
        # then restricted: verify it is in fact the padded target before using it.
        Sref = pad(az.trace_target(k, g1A, None if k == 2 else l2A), M)
        assert np.abs(S - Sref).max() < 1e-12, "padded target mismatch at rank %d, M %d" % (k, M)
        L = az.min_norm_trace(Sref, M, k)
        ins, out = split_norms(L, M, k)
        blk = L[tuple(slice(0, MA) for _ in range(2 * k))]
        if ref is None:
            ref = blk
        drift = float(np.abs(blk - ref).max())
        print("    %-6d %4d | %12.6f %12.3e | %12.3e %12.6f"
              % (k, M, ins, out / max(ins, 1e-30), drift / max(np.abs(ref).max(), 1e-30),
                 np.abs(lam).max()))
    print()

print("SCREEN 1  Lipschitz constant of each reconstruction, max over 120 directions at step 1e-4")
print()
SYS_STATES = []
dm = np.zeros((len(dets),) * 2, dtype=complex)
dm[0, 0] = 1.0
SYS_STATES.append(("pure ground (near boundary)", vec @ dm @ vec.conj().T))
w = np.exp(-(es - es[0]).real / 0.5)
w /= w.sum()
SYS_STATES.append(("thermal kT=0.5 (interior)", vec @ np.diag(w.astype(complex)) @ vec.conj().T))

rng = np.random.default_rng(11)
DELTA, NDIR = 1e-4, 120


def herm(M):
    A = rng.normal(size=(M, M)) + 1j * rng.normal(size=(M, M))
    A = 0.5 * (A + A.conj().T)
    return A / fro(A)


def antipair(M):
    T = rng.normal(size=(M,) * 4) + 1j * rng.normal(size=(M,) * 4)
    T = T - np.transpose(T, (1, 0, 2, 3))
    T = T - np.transpose(T, (0, 1, 3, 2))
    T = T + np.conj(np.transpose(T, (2, 3, 0, 1)))
    return T / fro(T)


for sname, dmat in SYS_STATES:
    g1, g2, g3 = rdms(dmat, MA, N, dets)
    l2 = az.cumulant(2, g2, {("G", 1): g1})
    occ = np.linalg.eigvalsh(0.5 * (g1 + g1.conj().T))
    print("  --- %s   occupations [%.4f .. %.4f] ---" % (sname, occ.min(), occ.max()))
    base2 = rank2_set(g1, N, MA)
    base3 = rank3_set(g1, l2, N, MA)
    lip2 = {nm: 0.0 for nm in base2}
    lip3 = {nm: 0.0 for nm in base3}
    for _ in range(NDIR):
        d1 = herm(MA)
        p2 = rank2_set(g1 + DELTA * d1, N, MA)
        for nm in base2:
            lip2[nm] = max(lip2[nm], fro(p2[nm] - base2[nm]) / DELTA)
        d2 = antipair(MA)
        p3 = rank3_set(g1 + DELTA * d1, l2 + DELTA * d2, N, MA)
        nrm = fro(DELTA * d1) + fro(DELTA * d2)
        for nm in base3:
            lip3[nm] = max(lip3[nm], fro(p3[nm] - base3[nm]) / nrm)
    print("    rank 2  " + "   ".join("%s %.3g" % (nm, v) for nm, v in lip2.items()))
    print("    rank 3  " + "   ".join("%s %.3g" % (nm, v) for nm, v in lip3.items()))
    print()
