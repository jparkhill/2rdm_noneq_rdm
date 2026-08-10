#!/usr/bin/env python3
"""What structure the exact cumulant has in non-interacting copies that the trace construction misses.

DEFINITIONS. For a composite of non-interacting fragments in a product state, the exact k-cumulant
lambda_k obeys two separate requirements, which this file keeps apart:
    CLAUSE 1, no leak:  lambda_k[p_1..p_k, q_1..q_k] = 0 whenever the 2k indices split across
                        fragments. This is the defining connectedness of a cumulant. It is a
                        statement about lambda, not gamma: gamma_2 = gamma_1 wedge gamma_1 +
                        lambda_2 must carry cross terms through the wedge.
    CLAUSE 2, no drift: lambda_k restricted to a fragment equals that fragment's own lambda_k,
                        with no dependence on the size, particle number or content of the rest.

    Tr_k T[p_1..p_{k-1}, q_1..q_{k-1}] = sum_y T[p_1..p_{k-1} y, q_1..q_{k-1} y]
    P = the orthogonal projector onto the row space of Tr_k, realised as
        P L = min_norm_trace(Tr_k L), verified orthogonal to L - P L at 1e-14 in
        tools/ansatz_validate.py. ||P lambda||^2 / ||lambda||^2 is exactly the fraction of the
        cumulant the single-trace identity can determine.

THREE QUESTIONS, each measured rather than argued.

Q1. IS THE TRACE PROJECTOR ITSELF EXTENSIVE? P is built from the global partial-trace map on M
    orbitals, so it need not commute with the fragment decomposition. If P lambda_exact leaks
    while lambda_exact does not, then the trace identity and extensivity cannot both be imposed
    by a clean projection, and every trace-consistent closure inherits a leak.

Q2. HOW MUCH DOES THE FRAGMENT-WISE SOLUTION BUY? Both of these satisfy the SAME global
    constraint Tr_k lambda_k = S:
        global    min_norm_trace(S, M, k)                  minimum norm over all of Mat_k
        blockwise direct sum over fragments f of min_norm_trace(S_f, M_f, k)
                                                           minimum norm within the extensive
                                                           subspace; identically the fragment
                                                           answer, embedded
    The blockwise one is extensive by construction. Comparing their errors against lambda_exact
    says whether extensivity is merely the absence of an artifact or an actual source of accuracy.

Q3. HOW DOES THE FAILURE SCALE WITH THE NUMBER OF COPIES? n non-interacting copies of the same
    fragment: the exact within-block magnitude is independent of n, while a global minimum-norm
    solution spreads the required trace over n times as many slots.

CONSTRUCTION. Fragment A = H3 bent / sto-3g, M_A = 6, N_A = 3, FCI ground state. Fragment B =
H2 / 0.9 Angstrom, M_B = 4, N_B = 2. For the composite the exact cumulants are the block direct
sums, verified against the explicitly constructed product state at 8.9e-16 in
tools/ansatz_extensivity.py, so they are used directly here. Ground truth qdiss.eom.fullci.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from ansatz_hierarchy import fro
from molecular_2rdm import build_H_sector, molecule, sector
from qdiss import ansatz as az
from qdiss.eom import fullci

H3_BENT = "H 0.00 0.00 0.0; H 1.00 0.00 0.0; H 0.35 0.90 0.0"
H2 = "H 0 0 0; H 0 0 0.9"


def frag(atom, basis, chg, spin, Nx):
    mol, mf, Mx, h, V, enuc = molecule(atom, basis, chg, spin)
    dets = sector(Mx, Nx)
    es, vec = np.linalg.eigh(build_H_sector(h, V, Mx, dets))
    psi = np.zeros(1 << Mx, dtype=complex)
    psi[dets] = vec[:, 0]
    g1 = fullci.rdm(psi, Mx, 1)
    l2 = az.cumulant(2, fullci.rdm(psi, Mx, 2), {("G", 1): g1})
    l3 = az.cumulant(3, fullci.rdm(psi, Mx, 3), {("G", 1): g1, ("L", 2): l2})
    return Mx, Nx, g1, l2, l3


def embed(T, M, off, k):
    """Place a 2k-index tensor living on orbitals [off, off+m) inside an M-orbital tensor."""
    out = np.zeros((M,) * (2 * k), dtype=complex)
    m = T.shape[0]
    out[tuple(slice(off, off + m) for _ in range(2 * k))] = T
    return out


def dsum(parts, M, k):
    """Block direct sum: parts is a list of (tensor, offset)."""
    out = np.zeros((M,) * (2 * k), dtype=complex)
    for T, off in parts:
        m = T.shape[0]
        out[tuple(slice(off, off + m) for _ in range(2 * k))] += T
    return out


def leak(T, M, k, bounds):
    """(max inside any single fragment, max on any tuple straddling fragments)."""
    A = np.abs(T)
    inside = 0.0
    keep = np.zeros((M,) * (2 * k), dtype=bool)
    for lo, hi in bounds:
        sl = tuple(slice(lo, hi) for _ in range(2 * k))
        inside = max(inside, float(A[sl].max()))
        keep[sl] = True
    return inside, float(A[~keep].max())


MA, NA, g1A, l2A, l3A = frag(H3_BENT, "sto-3g", 0, 1, 3)
MB, NB, g1B, l2B, l3B = frag(H2, "sto-3g", 0, 0, 2)
M, N = MA + MB, NA + NB
BOUNDS = ((0, MA), (MA, M))
g1 = dsum([(g1A, 0), (g1B, MA)], M, 1)
lam = {2: dsum([(l2A, 0), (l2B, MA)], M, 2), 3: dsum([(l3A, 0), (l3B, MA)], M, 3)}
lamA = {2: l2A, 3: l3A}
lamB = {2: l2B, 3: l3B}

print("composite H3 bent (M=6,N=3) + H2/0.9A (M=4,N=2)  ->  M=%d N=%d, non-interacting" % (M, N))
print()
print("Q0  is the partial trace map Tr_k even surjective?  if not, no exact trace inverse exists")
print("    %-4s %5s %14s %8s %10s" % ("rank", "M", "target dim", "rank", "surjective"))
for k in (2, 3):
    for Mx in (4, 6, 8, 10, 12):
        dim, rk = az.trace_rank_deficiency(Mx, k)
        print("    %-4d %5d %14d %8d %10s" % (k, Mx, dim, rk, "yes" if rk == dim else "NO"))
print()
print("Q1  is the trace projector P extensive?  P lambda = min_norm_trace(Tr_k lambda)")
print("    %-4s %-22s %11s %11s %11s" % ("rank", "object", "|inside|", "leak", "captured"))
for k in (2, 3):
    ex = lam[k]
    ins, lk = leak(ex, M, k, BOUNDS)
    print("    %-4d %-22s %11.6f %11.3e %11s" % (k, "lambda_k exact", ins, lk, "--"))
    Pex = az.project_trace(ex, M, k)
    ins, lk = leak(Pex, M, k, BOUNDS)
    print("    %-4d %-22s %11.6f %11.3e %10.1f%%"
          % (k, "P lambda_k  (global)", ins, lk / max(ins, 1e-30),
             100.0 * (fro(Pex) / fro(ex)) ** 2))
    PA = az.project_trace(lamA[k], MA, k)
    PB = az.project_trace(lamB[k], MB, k)
    Pblk = dsum([(PA, 0), (PB, MA)], M, k)
    ins, lk = leak(Pblk, M, k, BOUNDS)
    print("    %-4d %-22s %11.6f %11.3e %10.1f%%"
          % (k, "P lambda_k  (blockwise)", ins, lk / max(ins, 1e-30),
             100.0 * (fro(Pblk) / fro(ex)) ** 2))
print()

print("Q2  two solutions of the SAME constraint Tr_k lambda_k = S, one extensive one not")
print("    %-4s %-22s %11s %11s %11s" % ("rank", "solution", "rel err", "leak", "trace res"))
for k in (2, 3):
    ex = lam[k]
    S = az.trace_target(k, g1, None if k == 2 else lam[2])
    SA = az.trace_target(k, g1A, None if k == 2 else l2A)
    SB = az.trace_target(k, g1B, None if k == 2 else l2B)
    glob = az.min_norm_trace(S, M, k)
    blk = dsum([(az.min_norm_trace(SA, MA, k), 0),
                (az.min_norm_trace(SB, MB, k), MA)], M, k)
    for nm, L in (("global min-norm", glob), ("blockwise min-norm", blk),
                  ("null (lambda_k = 0)", np.zeros_like(ex))):
        ins, lk = leak(L, M, k, BOUNDS)
        print("    %-4d %-22s %11.4f %11.3e %11.2e"
              % (k, nm, fro(L - ex) / fro(ex), lk / max(ins, 1e-30),
                 az.trace_residual(L, S, k)))
print()

print("Q3  n non-interacting copies of H3 bent: within-block magnitude vs n")
print("    exact is independent of n by clause 2; global min-norm spreads over n times the slots")
print("    %-4s %3s %5s | %13s %13s %13s"
      % ("rank", "n", "M", "exact |blk|", "global |blk|", "ratio glob/ex"))
for k, LA, NC in ((2, l2A, (1, 2, 3, 4)), (3, l3A, (1, 2))):
    for n in NC:
        Mn = MA * n
        g1n = dsum([(g1A, i * MA) for i in range(n)], Mn, 1)
        l2n = dsum([(l2A, i * MA) for i in range(n)], Mn, 2)
        Sn = az.trace_target(k, g1n, None if k == 2 else l2n)
        G = az.min_norm_trace(Sn, Mn, k)
        gb = float(np.abs(G[tuple(slice(0, MA) for _ in range(2 * k))]).max())
        eb = float(np.abs(LA).max())
        print("    %-4d %3d %5d | %13.6f %13.6f %13.4f" % (k, n, Mn, eb, gb, gb / eb))
    print()
