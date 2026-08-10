#!/usr/bin/env python3
"""Can an EXTENSIVE lambda_3 satisfy the trace identity? The connected-contraction space.

THE ARGUMENT BEING TESTED. A reconstruction term leaks across non-interacting fragments exactly
when it is a DISCONNECTED contraction. If two factors share no summed index, each is confined to
its own fragment by block-diagonality and their outer product puts weight on an index tuple that
straddles the two, which the exact cumulant forbids. If they share an index, block-diagonality of
each factor drags all their indices into one fragment. So:

    lambda_k is extensive  <=>  every term is a connected contraction, i.e. all 2k external
    indices are joined through a path of shared internal indices.

Counting free indices then gives a rank-dependent no-go. Writing (u,l) for upper and lower index
counts: gamma_1 is (1,1), lambda_2 is (2,2), V is (2,2).
    rank 2 needs (2,2). A connected diagram of gamma_1 vertices alone is a chain
        gamma_1 gamma_1 ... gamma_1, always (1,1). So NO gamma_1-only connected lambda_2 exists,
        and every gamma_1-only reconstruction must wedge, hence must leak. V is (2,2) and
        connected, at order V^1, which is the order of lambda_2. V is not optional at rank 2.
    rank 3 needs (3,3). Two lambda_2 sharing one index give (2+2-1, 2+2-1) = (3,3), connected,
        at order V^2, which is the order of lambda_3. So rank 3 DOES admit a gamma-only extensive
        reconstruction and rank 2 does not.
Dimensional analysis narrows it further: lambda_k is dimensionless while V is an energy, so a bare
V term cannot appear without an energy denominator. lambda_2 * lambda_2 is the only dimensionless
gamma-only connected rank-3 structure at leading order.

WHAT IS MEASURED HERE.
  1. leak of every basis structure on a composite of non-interacting fragments. Must be 0.
  2. best least-squares fit of lambda_3^exact inside the connected span: the fraction of the
     exact 3-cumulant that an extensive reconstruction can reach at all.
  3. whether the trace identity Tr_3 lambda_3 = anti(E, lambda_2) can be met inside that span.
     This is the crux. The target is formally O(V), while every connected structure is O(V^2);
     at an idempotent reference the four E eigenvalues cancel on lambda_2's amplitude block and
     the target drops to O(V^2), but away from idempotency it need not, and then no extensive
     trace-consistent lambda_3 exists in this space.

SYMBOLS. lambda_2[p v, q w] with p,v upper and q,w lower; E = gamma_1 - I/2; anti(E, X) is E
contracted into each of X's slots and summed; Tr_3 T[pv,qw] = sum_y T[p v y, q w y]. A3 is the
full antisymmetriser over the three upper and the three lower indices, sum of 36 signed
permutations divided by 36. Ground truth qdiss.eom.fullci.
"""
import itertools
import sys
import numpy as np
import scipy.linalg as sla
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from ansatz_hierarchy import fro
from molecular_2rdm import build_H_sector, molecule, sector
from qdiss import ansatz as az
from qdiss.eom import fullci

H3_BENT = "H 0.00 0.00 0.0; H 1.00 0.00 0.0; H 0.35 0.90 0.0"
H2 = "H 0 0 0; H 0 0 0.9"
H4_TRAP = "H 0.00 0.00 0.0; H 0.80 0.00 0.0; H 1.05 0.72 0.0; H -0.15 0.72 0.0"

_P3 = [(p, (-1) ** sum(1 for i in range(3) for j in range(i + 1, 3) if p[i] > p[j]))
       for p in itertools.permutations(range(3))]


def A3(T):
    out = np.zeros_like(T)
    for pu, su in _P3:
        for pl, sl in _P3:
            out += (su * sl) * T.transpose(tuple(pu) + tuple(3 + i for i in pl))
    return out / 36.0


def chains(g1, l2, M):
    """A basis of connected rank-3 structures built from gamma_1 and lambda_2.

    Every entry is a single contraction between two lambda_2 factors, optionally with gamma_1
    inserted on the contraction line or contracted into a free slot. Generated redundantly and
    then orthogonalised, so the reported dimension is the true dimension of the connected span.
    """
    I = np.eye(M, dtype=complex)
    E = g1 - 0.5 * I
    eta = I - g1
    raw = []
    for nm, D in (("bare", I), ("g1", g1), ("eta", eta)):
        # lower index of the first factor joined to an upper index of the second, through D
        raw.append(("chain_%s" % nm,
                    np.einsum("pvqr,rs,sxwy->pvxqwy", l2, D, l2, optimize=True)))
        # the adjoint pattern: upper of the first joined to a lower of the second
        raw.append(("chainT_%s" % nm,
                    np.einsum("prqw,rs,vxsy->pvxqwy", l2, D, l2, optimize=True)))
    base = []
    for nm, T in raw:
        X = A3(T)
        X = 0.5 * (X + np.conj(X.transpose(3, 4, 5, 0, 1, 2)))   # pair-hermitian part
        base.append((nm, X))
        Y = _anti3(E, X)
        base.append(("anti(E)." + nm, 0.5 * (Y + np.conj(Y.transpose(3, 4, 5, 0, 1, 2)))))
    return base


def _anti3(E, X):
    """E contracted into each of the six slots of a rank-3 tensor, summed."""
    out = np.einsum("pr,rvxqwy->pvxqwy", E, X)
    out += np.einsum("vr,prxqwy->pvxqwy", E, X)
    out += np.einsum("xr,pvrqwy->pvxqwy", E, X)
    out += np.einsum("pvxrwy,rq->pvxqwy", X, E)
    out += np.einsum("pvxqry,rw->pvxqwy", X, E)
    out += np.einsum("pvxqwr,ry->pvxqwy", X, E)
    return out


def orthonormalise(base, tol=1e-9):
    V = np.array([T.reshape(-1) for _, T in base])
    G = V.conj() @ V.T
    ev, U = np.linalg.eigh(G)
    keep = ev > tol * max(ev.max().real, 1.0)
    W = (U[:, keep] / np.sqrt(ev[keep])).T @ V
    return W, int(keep.sum())


def frag(atom, basis, chg, spin, Nx):
    mol, mf, Mx, h, V, enuc = molecule(atom, basis, chg, spin)
    dets = sector(Mx, Nx)
    es, vec = np.linalg.eigh(build_H_sector(h, V, Mx, dets))
    psi = np.zeros(1 << Mx, dtype=complex)
    psi[dets] = vec[:, 0]
    g1 = fullci.rdm(psi, Mx, 1)
    l2 = az.cumulant(2, fullci.rdm(psi, Mx, 2), {("G", 1): g1})
    l3 = az.cumulant(3, fullci.rdm(psi, Mx, 3), {("G", 1): g1, ("L", 2): l2})
    return Mx, Nx, g1, l2, l3, es[0].real + enuc


def leak(T, M, k, bounds):
    A = np.abs(T)
    inside, keep = 0.0, np.zeros((M,) * (2 * k), dtype=bool)
    for lo, hi in bounds:
        sl = tuple(slice(lo, hi) for _ in range(2 * k))
        inside = max(inside, float(A[sl].max()))
        keep[sl] = True
    return inside, float(A[~keep].max())


# ---------------------------------------------------------------- 1. extensivity of the basis
MA, NA, g1A, l2A, l3A, eA = frag(H3_BENT, "sto-3g", 0, 1, 3)
MB, NB, g1B, l2B, l3B, eB = frag(H2, "sto-3g", 0, 0, 2)
M, N = MA + MB, NA + NB


def dsum(pairs, Mx, k):
    out = np.zeros((Mx,) * (2 * k), dtype=complex)
    for T, off in pairs:
        m = T.shape[0]
        out[tuple(slice(off, off + m) for _ in range(2 * k))] += T
    return out


g1c = dsum([(g1A, 0), (g1B, MA)], M, 1)
l2c = dsum([(l2A, 0), (l2B, MA)], M, 2)
l3c = dsum([(l3A, 0), (l3B, MA)], M, 3)
BOUNDS = ((0, MA), (MA, M))

print("1. EXTENSIVITY of every connected structure, on H3 + H2 non-interacting, M=%d N=%d" % (M, N))
print("   the exact lambda_3 has leak %.2e; every connected structure must match that"
      % (leak(l3c, M, 3, BOUNDS)[1] / leak(l3c, M, 3, BOUNDS)[0]))
print("   %-22s %13s %13s" % ("structure", "|inside|", "leak/inside"))
basec = chains(g1c, l2c, M)
for nm, T in basec:
    ins, lk = leak(T, M, 3, BOUNDS)
    print("   %-22s %13.3e %13.3e" % (nm, ins, lk / max(ins, 1e-30)))
# the wedge-based trace inverse, for contrast
S = az.trace_target(3, g1c, l2c)
ins, lk = leak(az.min_norm_trace(S, M, 3), M, 3, BOUNDS)
print("   %-22s %13.3e %13.3e" % ("(min-norm trace inv)", ins, lk / max(ins, 1e-30)))
print()

# ---------------------------------------------------------------- 2 & 3. reach and trace identity
print("2/3. REACH of the connected span, coefficient universality, and the parameter-free closure")
print()
print("   Two fits are reported and they are NOT the same thing:")
print("     ORACLE   coefficients least-squares fitted to lambda_3^exact. This is an UPPER BOUND")
print("              on what the span can reach, not a closure -- it has seen the answer.")
print("     CLOSURE  coefficients determined ONLY by requiring Tr_3 lambda_3 = anti(E, lambda_2),")
print("              which never looks at lambda_3. This IS a closure, with no free parameters.")
print("   If the ORACLE coefficients agree across systems and states, the span is universal; if")
print("   they do not, the small oracle error is a fitting artifact and means nothing.")
print()

# A fixed, state-independent subset of the generated structures, chosen once by QR pivoting so
# that coefficients are comparable across systems.
_probeM, _probeN, _pg1, _pl2, _pl3, _ = frag(H4_TRAP, "sto-3g", 0, 0, 4)
_pb = chains(_pg1, _pl2, _probeM)
_Vp = np.array([T.reshape(-1) for _, T in _pb])
_q, _r, _p = sla.qr(_Vp.T, pivoting=True, mode="economic")
_rank = int((np.abs(np.diag(_r)) > 1e-9 * abs(_r[0, 0])).sum())
PICK = sorted(_p[:_rank])
print("   independent structure subset (chosen once, by QR pivoting on H4): %s"
      % ", ".join(_pb[i][0] for i in PICK))
print()

CASES = []
for atom, basis, chg, spin, Nx, label in ((H3_BENT, "sto-3g", 0, 1, 3, "H3 bent"),
                                          (H4_TRAP, "sto-3g", 0, 0, 4, "H4 trapezoid")):
    Mx, Nx2, g1, l2, l3, e0 = frag(atom, basis, chg, spin, Nx)
    CASES.append((label + " (pure)", Mx, g1, l2, l3))

print("   %-22s | %s | %9s %9s | %9s %9s"
      % ("system", "  ".join("%9s" % _pb[i][0][:9] for i in PICK),
         "oracle", "err", "closure", "err"))
for label, Mx, g1, l2, l3 in CASES:
    base = chains(g1, l2, Mx)
    B = np.array([base[i][1].reshape(-1) for i in PICK])
    Sx = az.trace_target(3, g1, l2)
    # ORACLE: fit to the exact 3-cumulant
    co, *_ = np.linalg.lstsq(B.T, l3.reshape(-1), rcond=None)
    fo = (co @ B).reshape(l3.shape)
    # CLOSURE: fit only the partial trace, never seeing lambda_3
    Tr = np.array([np.einsum("pvyqwy->pvqw", b.reshape(l3.shape)).reshape(-1) for b in B])
    cc, *_ = np.linalg.lstsq(Tr.T, Sx.reshape(-1), rcond=None)
    fc = (cc @ B).reshape(l3.shape)
    print("   %-22s | %s | %9s %9.4f | %9s %9.4f"
          % (label, "  ".join("%+9.4f" % c.real for c in co),
             "", fro(fo - l3) / fro(l3), "", fro(fc - l3) / fro(l3)))
    print("   %-22s | %s | trace res oracle %.2e  closure %.2e   target |S| %.2e"
          % ("  closure coeffs", "  ".join("%+9.4f" % c.real for c in cc),
             az.trace_residual(fo, Sx, 3), az.trace_residual(fc, Sx, 3), np.abs(Sx).max()))
print()
print("   reference errors on the same states, from tools/ansatz_hierarchy.py:")
print("     null (lambda_3 = 0) 1.0000    min-norm trace 0.79-0.82    Mueller+Lackner 2.25-2.73")
