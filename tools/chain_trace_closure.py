#!/usr/bin/env python3
r"""Can the trace identity be satisfied EXACTLY inside the connected span? Common coefficients.

THE UNFINISHED PIECE. The connected-contraction rule says every term of lambda_3 must be a single
connected contraction, and leg counting forces exactly two lambda_2 factors joined through one index
(a tree). The trace identity demands

    Tr_3 lambda_3 = anti(E, lambda_2),      E = gamma_1 - I/2,
    Tr_3 T[pv,qw] = sum_y T[p v y, q w y],  anti(E,X) = E contracted into each of X's 4 slots

Tracing a bilinear-in-lambda_2 structure leaves something bilinear, while the target is linear in
lambda_2. Some terms become effectively linear once the rank-2 identity sum_y lambda_2[ry,wy] =
-(gamma_1 eta_1)[r,w] is used to close a loop; the rest is a genuine BILINEAR RESIDUE. I asserted
that residue could be absorbed and never checked it. tools/connected_chains.py measured a residual
of 7.1e-05 against a target norm of 3.5e-02 at H3 and 8.4e-04 against 5.0e-02 at H4 -- 0.2% to 1.7%,
NOT zero -- but that was with only four structures.

THIS FILE enlarges the connected basis and asks whether the residual goes to zero with coefficients
COMMON ACROSS SYSTEMS AND STATES. A per-state fit proves nothing: this project has already been
burned by fitting a bare scalar multiple of one structure and reading the answer off a single system.

THE BASIS, all connected, all bilinear in lambda_2, all dimensionless:
    chain_a  = A3[ sum_{r,s} lambda_2[p v, q r] (gamma_1^a)[r,s] lambda_2[s x, w y] ]
    chainT_a = the adjoint contraction pattern
    then dressed on the free legs by anti(gamma_1^b, .), b = 0, 1, 2
for a = 0..3, giving 2 x 4 x 3 = 24 structures, each pair-Hermitised. A3 is the antisymmetriser over
the three upper and three lower indices (36 signed terms over 36).

THE FIT. For each state, Tr_3 of every basis structure is stacked as a row block and the target
anti(E, lambda_2) as the right-hand side; the blocks from all states are stacked vertically and one
common coefficient vector is solved for by least squares. Reported: the residual with common
coefficients, the residual allowed to vary per state (a lower bound on what any common fit can do),
and the singular spectrum of the design matrix so rank deficiency is visible rather than hidden.

Ground truth qdiss.eom.fullci.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from ansatz_hierarchy import fro, rdms
from molecular_2rdm import build_H_sector, molecule, sector
from qdiss import ansatz as az
from qdiss.eom import fullci, hamiltonian

H3 = "H 0.00 0.00 0.0; H 1.00 0.00 0.0; H 0.35 0.90 0.0"
H4 = "H 0.00 0.00 0.0; H 0.80 0.00 0.0; H 1.05 0.72 0.0; H -0.15 0.72 0.0"
AMAX, BMAX = 4, 3


def basis_of(g1, l2, M):
    """Connected rank-3 structures: two lambda_2 joined through gamma_1^a, legs dressed by
    anti(gamma_1^b, .). Names carried so the fitted coefficients are interpretable."""
    pw = [np.eye(M, dtype=complex)]
    for _ in range(AMAX - 1):
        pw.append(pw[-1] @ g1)
    out = []
    for a in range(AMAX):
        for tag, sub in (("chain", "pvqr,rs,sxwy->pvxqwy"), ("chainT", "prqw,rs,vxsy->pvxqwy")):
            X = az.antisym3(np.einsum(sub, l2, pw[a], l2, optimize=True))
            X = 0.5 * (X + np.conj(X.transpose(3, 4, 5, 0, 1, 2)))
            for b in range(BMAX):
                Y = X
                for _ in range(b):
                    Y = az.anti3(g1, Y)
                Y = 0.5 * (Y + np.conj(Y.transpose(3, 4, 5, 0, 1, 2)))
                out.append(("%s_g%d.anti%d" % (tag, a, b), Y))
    return out


def state_data(atom, basis, N, spin, label, kt=None, seed=None):
    if atom is None:                                  # random Hamiltonian control
        rng = np.random.default_rng(seed)
        M = 6
        h, V = hamiltonian.random_h(M, rng), hamiltonian.random_V(M, rng)
        psi = fullci.random_state(M, N, rng)
        g1 = fullci.rdm(psi, M, 1)
        l2 = az.cumulant(2, fullci.rdm(psi, M, 2), {("G", 1): g1})
        l3 = az.cumulant(3, fullci.rdm(psi, M, 3), {("G", 1): g1, ("L", 2): l2})
        return label, M, g1, l2, l3
    mol, mf, M, h, V, enuc = molecule(atom, basis, 0, spin)
    dets = sector(M, N)
    es, vec = np.linalg.eigh(build_H_sector(h, V, M, dets))
    if kt is None:
        dm = np.zeros((len(dets),) * 2, dtype=complex)
        dm[0, 0] = 1.0
    else:
        w = np.exp(-(es - es[0]).real / kt)
        w /= w.sum()
        dm = np.diag(w.astype(complex))
    g = rdms(vec @ dm @ vec.conj().T, M, N, dets)
    l2 = az.cumulant(2, g[1], {("G", 1): g[0]})
    l3 = az.cumulant(3, g[2], {("G", 1): g[0], ("L", 2): l2})
    return label, M, g[0], l2, l3


CASES = [state_data(H3, "sto-3g", 3, 1, "H3 pure"),
         state_data(H3, "sto-3g", 3, 1, "H3 kT=0.5", kt=0.5),
         state_data(H4, "sto-3g", 4, 0, "H4 pure"),
         state_data(H4, "sto-3g", 4, 0, "H4 kT=0.5", kt=0.5),
         state_data(None, None, 3, None, "random h,V (6,3)", seed=23)]

rows, rhs, per_state, names = [], [], [], None
print("Can Tr_3(connected span) = anti(E, lambda_2) with COMMON coefficients?")
print("basis: 2 topologies x gamma_1^a (a<%d) x anti(gamma_1^b) (b<%d) = %d structures"
      % (AMAX, BMAX, 2 * AMAX * BMAX))
print()
print("%-18s %3s | %11s %11s | %11s"
      % ("state", "M", "|target|", "per-state res", "rel"))
for label, M, g1, l2, l3 in CASES:
    bas = basis_of(g1, l2, M)
    names = [n for n, _ in bas]
    S = az.trace_target(3, g1, l2)
    A = np.array([np.einsum("pvyqwy->pvqw", X).reshape(-1) for _, X in bas]).T
    rows.append(A)
    rhs.append(S.reshape(-1))
    c, *_ = np.linalg.lstsq(A, S.reshape(-1), rcond=None)
    r = np.linalg.norm(A @ c - S.reshape(-1))
    per_state.append(r)
    print("%-18s %3d | %11.4e %11.4e | %11.4e"
          % (label, M, np.linalg.norm(S), r, r / np.linalg.norm(S)))

Ab = np.vstack(rows)
bb = np.concatenate(rhs)
c, *_ = np.linalg.lstsq(Ab, bb, rcond=None)
res = np.linalg.norm(Ab @ c - bb)
sv = np.linalg.svd(Ab, compute_uv=False)
print()
print("COMMON fit over %d states: residual %.4e against target norm %.4e  (rel %.4e)"
      % (len(CASES), res, np.linalg.norm(bb), res / np.linalg.norm(bb)))
print("  best possible if coefficients were allowed to vary per state: %.4e"
      % np.linalg.norm(per_state))
print("  design matrix: %d equations x %d unknowns, singular values %.2e .. %.2e, "
      "numerical rank %d" % (Ab.shape[0], Ab.shape[1], sv[0], sv[-1],
                             int((sv > 1e-10 * sv[0]).sum())))
print()
print("  coefficients above 1e-3 of the largest:")
big = np.abs(c) > 1e-3 * np.abs(c).max()
for nm, ci in zip(np.array(names)[big], c[big]):
    print("    %-22s %+12.5f" % (nm, ci.real))
