#!/usr/bin/env python3
"""Reference-free adiabatic closure for lambda_2: partition the MKNO rank-2 EOM by cumulant
degree, then solve the stationarity / Hylleraas condition.

DERIVATION. qdiss.eom.gwt.terms_cumulant(2) gives the MKNO rank-2 EOM in gamma_1 and cumulants,
    lambda_2-dot = L[lambda_2] + D[gamma_1] + (lambda_3 terms) + (lambda_2^2 terms)
where the split is purely by how many ("L",2) factors a term carries: degree 1 -> L, degree 0 -> D.
lambda_3 is TRUNCATED to zero and the degree-2 terms are dropped, which is what makes the
remaining problem linear. Adiabatic elimination (lambda_2-dot = 0) then gives

    L[lambda_2] = -D        <=>     lambda_2 stationary point of  J = 1/2 <l2, L l2> + <l2, D>

the Hylleraas functional. NO REFERENCE DETERMINANT ENTERS: L and D are built from h, V and
gamma_1 alone, the occupied/virtual projectors of ordinary MP2 being replaced by gamma_1 and
eta_1 = 1 - gamma_1 inside the MKNO term list, and the energy denominators by the spectrum of the
superoperator L built from the generalised Fock operator f = h + V.gamma_1.

ACCEPTANCE GATE. As gamma_1 becomes idempotent this must reduce to ordinary MP2: the correlation
energy (1/4) V.lambda_2 must equal PySCF's MP2 correlation energy. That is an independent check of
the whole construction and is reported first.

L is singular on part of the pair space (zero denominators / the gap condition), so the solve is a
least-squares one with an explicit rank report -- no silent pseudo-inverse.
"""
import sys
import numpy as np
sys.path.insert(0, "."); sys.path.insert(0, "tools")
import stationary_lagrangian as S
from pyscf import mp
from qdiss import closures
from qdiss.eom import fullci, gwt

M, N, h, V, enuc = S.M, S.N, S.h, S.V, S.enuc
g_hf, g_fci, E_HF, E_FCI = S.g_hf, S.g_fci, S.E_HF, S.E_FCI
EYE = np.eye(M, dtype=complex)
L2KEY, L3KEY = ("L", 2), ("L", 3)


def split_terms():
    """Partition the rank-2 MKNO EOM term list by degree in lambda_2."""
    one, two = gwt.terms_cumulant(2)
    deg = {0: ([], []), 1: ([], [])}
    n_l3 = n_hi = 0
    for which, lst in ((0, one), (1, two)):
        for s, f in lst:
            keys = [(n, r) for n, r, u, l in f]
            if L3KEY in keys:
                n_l3 += 1; continue
            d = keys.count(L2KEY)
            if d > 1:
                n_hi += 1; continue
            subs = ["".join(u) + "".join(l) for n, r, u, l in f]
            deg[d][which].append((s, subs, keys))
    return deg, n_l3, n_hi


DEG, N_L3, N_HI = split_terms()

# Second partition of the degree-1 terms, which is what makes the expansion perturbative rather
# than a self-consistent resummation. Count how many of the interaction operator's contracted
# indices the lambda_2 factor shares:
#   ONE  -> the operator reaches lambda_2 through a single line, exactly as h does in the one-body
#           terms (h[r,s] d[p,s] l2[rvqw]). With V this is the mean-field/Fock channel, its other
#           two indices closing on a gamma_1. These belong in L0.
#   TWO  -> V acts inside the pair, the ladder/ring channel. First order in V beyond L0, so it is
#           NOT part of the zeroth-order operator; keeping it turns MP1 into an RPA-like
#           resummation, which is why the first gate attempt missed MP2 by a factor of -1.32.
OPIDX = {0: set("rs"), 1: set("rtsu")}


def _shared(subs, keys, which):
    for sub, k in zip(subs, keys):
        if k == L2KEY:
            return len(OPIDX[which] & set(sub))
    return 0


L0 = ([], [])
LAD = ([], [])
for which in (0, 1):
    for s, subs, keys in DEG[1][which]:
        (L0 if _shared(subs, keys, which) <= 1 else LAD)[which].append((s, subs, keys))
print("second partition of the degree-1 terms by how many operator indices lambda_2 touches:")
print("   L0  (mean-field/Fock channel, 1 shared index): %3d one-body + %3d two-body"
      % (len(L0[0]), len(L0[1])))
print("   ladder (2 shared indices, first order in V)  : %3d one-body + %3d two-body"
      % (len(LAD[0]), len(LAD[1])))
print()


def apply_L0(g1, X):
    return build(L0, {("D", 1): EYE, ("G", 1): g1, L2KEY: X})
print("term partition of the MKNO rank-2 EOM by degree in lambda_2:")
print("   degree 0 (driving D) : %3d one-body + %3d two-body"
      % (len(DEG[0][0]), len(DEG[0][1])))
print("   degree 1 (operator L): %3d one-body + %3d two-body"
      % (len(DEG[1][0]), len(DEG[1][1])))
print("   dropped: %d terms containing lambda_3 (truncated), %d nonlinear in lambda_2"
      % (N_L3, N_HI))
print()


def build(terms, tab):
    one, two = terms
    tot = np.zeros((M,) * 4, dtype=complex)
    for s, subs, keys in one:
        tot += s * np.einsum(",".join(["rs"] + subs) + "->pvqw", h,
                             *[tab[k] for k in keys], optimize=True)
    for s, subs, keys in two:
        tot += s * 0.25 * np.einsum(",".join(["rtsu"] + subs) + "->pvqw", V,
                                    *[tab[k] for k in keys], optimize=True)
    # The EOM is lambda_2-dot = i * tot. The i is NOT optional here: tot is anti-Hermitian in the
    # pair index, so without it L maps into the anti-Hermitian sector, exactly orthogonal to the
    # Hermitian pair basis, and the projection returns zero for every basis element.
    return 1j * tot


def driving(g1):
    return build(DEG[0], {("D", 1): EYE, ("G", 1): g1})


def apply_L(g1, X):
    return build(DEG[1], {("D", 1): EYE, ("G", 1): g1, L2KEY: X})


PAIRS = [(p, v) for p in range(M) for v in range(p + 1, M)]
NP = len(PAIRS)


def pair_to_tensor(Z):
    T = np.zeros((M,) * 4, dtype=complex)
    for i, (p, v) in enumerate(PAIRS):
        for j, (q, w) in enumerate(PAIRS):
            x = 0.5 * Z[i, j]
            T[p, v, q, w] = x; T[v, p, q, w] = -x
            T[p, v, w, q] = -x; T[v, p, w, q] = x
    return T


BASIS = []
for a in range(NP):
    for b in range(a, NP):
        for ph in (1.0, 1j):
            Z = np.zeros((NP, NP), dtype=complex)
            Z[a, b] += ph; Z[b, a] += np.conj(ph)
            if np.abs(Z).max() < 1e-12:
                continue
            BASIS.append(pair_to_tensor(Z))
NB = len(BASIS)
GRAM = np.array([np.real(np.tensordot(b.conj(), b, axes=4)) for b in BASIS])


def coeff(T):
    return np.array([np.real(np.tensordot(b.conj(), T, axes=4)) / g
                     for b, g in zip(BASIS, GRAM)])


def solve_lambda2(g1, tag, op=None):
    op = op or apply_L
    Lm = np.zeros((NB, NB))
    for k, bk in enumerate(BASIS):
        Lm[:, k] = coeff(op(g1, bk))
    d = coeff(driving(g1))
    x, res, rank, sv = np.linalg.lstsq(Lm, -d, rcond=1e-8)
    l2 = sum(xi * b for xi, b in zip(x, BASIS))
    resid = np.abs(op(g1, l2) + driving(g1)).max()
    nz = sv[sv > 1e-12]
    print("   %-22s NB %d  rank(L) %d  |L|_max %.2e  sigma_max %.2e  sigma_min>0 %s"
          "  |D| %.2e  |L l2 + D| %.2e"
          % (tag, NB, rank, np.abs(Lm).max(), sv.max() if len(sv) else 0.0,
             ("%.2e" % nz.min()) if len(nz) else "none", np.abs(driving(g1)).max(), resid))
    return l2


def ecorr(l2):
    return float(0.25 * np.real(np.einsum("pvqw,pvqw->", V, l2)))


# ---- trace_plus_pt: the pt amplitude plus a ker(L0) component carrying the trace identity ----
# Both halves are verified orthogonal: the amplitude lives in range(L0) and contributes exactly
# nothing to the contraction (index disjointness), while the correction lives in ker(L0) and so
# cannot disturb the stationarity condition L0[lambda_2] = -D. Adding them is therefore not a
# compromise -- each metric is fixed by the piece that can see it.
HBn_np = [b / np.sqrt(np.real((b.conj() * b).sum())) for b in S.HB]


def herm_coeff(X):
    return np.array([np.real(np.tensordot(b.conj(), X, axes=([0, 1], [0, 1]))) for b in HBn_np])


def kernel_basis(g1, tol=1e-8):
    """Orthonormal basis of ker(L0) in the antisymmetric Hermitian pair basis."""
    Lm = np.zeros((NB, NB))
    for k, bk in enumerate(BASIS):
        Lm[:, k] = coeff(apply_L0(g1, bk))
    u, s, vt = np.linalg.svd(Lm)
    keep = s <= tol * s.max()
    return vt[keep], s


def trace_plus_pt(g1, l2_pt, verbose=True):
    """Add W in ker(L0) so that sum_x lambda_2[px,qx] = -(gamma_1 - gamma_1^2) exactly."""
    K, s = kernel_basis(g1)
    nk = K.shape[0]
    Kt = [sum(K[a, k] * BASIS[k] for k in range(NB)) for a in range(nk)]
    Amat = g1 - g1 @ g1
    target = herm_coeff(-Amat - np.einsum("pxqx->pq", l2_pt))
    Pm = np.array([herm_coeff(np.einsum("pxqx->pq", w)) for w in Kt]).T
    c, _, rk, _ = np.linalg.lstsq(Pm, target, rcond=None)
    W = sum(ci * w for ci, w in zip(c, Kt))
    if verbose:
        gap = s[s > 1e-8 * s.max()].min() / s.max()
        print("   ker(L0) dim %d of %d (relative gap %.1e); partial-trace map rank %d of %d"
              % (nk, NB, gap, rk, Pm.shape[0]))
        print("   |residual of the trace target within ker(L0)| = %.2e"
              % np.abs(Pm @ c - target).max())
        print("   |W| = %.5f   |L0[W]| = %.2e  (must be zero: W is in the kernel)"
              % (np.abs(W).max(), np.abs(apply_L0(g1, W)).max()))
    return l2_pt + W


if __name__ == "__main__":
    e_mp2 = mp.MP2(S.mf).kernel()[0]
    print("ACCEPTANCE GATE: at the idempotent (HF) density this must reproduce MP2")
    l2_hf = solve_lambda2(g_hf, "L0 only (MP1)", apply_L0)
    l2_hf_full = solve_lambda2(g_hf, "L0 + ladder (RPA-like)", apply_L)
    print("   PySCF MP2 E_corr                      = %+.9f" % e_mp2)
    print("   (1/4) V.lambda_2  from L0 only         = %+.9f   ratio %+.6f"
          % (ecorr(l2_hf), ecorr(l2_hf) / e_mp2))
    print("   (1/4) V.lambda_2  from L0 + ladder     = %+.9f   ratio %+.6f"
          % (ecorr(l2_hf_full), ecorr(l2_hf_full) / e_mp2))
    print()

    g2ex = fullci.rdm(S.psi0, M, 2)
    l2ex = g2ex - closures.wedge(g_fci, g_fci)
    print("at the correlated FCI density, against the exact cumulant (|l2_ex| = %.5f):"
          % np.abs(l2ex).max())
    l2_c = solve_lambda2(g_fci, "at g_FCI, L0 only", apply_L0)
    A = g_fci - g_fci @ g_fci
    print()
    print("building trace_plus_pt at the FCI density:")
    l2_tpt = trace_plus_pt(g_fci, l2_c)
    print("   |L0[l2] + D| after adding W = %.2e  (stationarity preserved)"
          % np.abs(apply_L0(g_fci, l2_tpt) + driving(g_fci)).max())
    print()
    PAIRS_ = PAIRS

    def mineig2(T):
        Aq = np.array([[2.0 * T[p, v, q, w] for (q, w) in PAIRS_] for (p, v) in PAIRS_])
        return np.linalg.eigvalsh(0.5 * (Aq + Aq.conj().T)).min()

    print("   %-14s %10s %11s %11s %10s %10s"
          % ("closure", "|l2-l2ex|", "trace-id", "E-E_FCI", "tr g2 err", "min eig D"))
    for tag, l2 in (("trace_plus_pt", l2_tpt), ("pt", l2_c),
                    ("blocked", closures.lambda2(g_fci, "blocked")),
                    ("blocked+pt", closures.lambda2(g_fci, "blocked") + l2_c),
                    ("mp2", closures.lambda2(g_fci, "mp2", h, V)),
                    ("zero", np.zeros((M,) * 4, dtype=complex))):
        g2 = closures.wedge(g_fci, g_fci) + l2
        print("   %-14s %10.5f %11.5f %+11.6f %10.2e %10.6f"
              % (tag, np.abs(l2 - l2ex).max(),
                 np.abs(np.einsum("pxqx->pq", l2) + A).max(),
                 float(np.real(np.einsum("pq,pq->", h, g_fci)
                               + 0.25 * np.einsum("pvqw,pvqw->", V, g2))) + enuc - E_FCI,
                 abs(np.einsum("pvpv->", g2).real - N * (N - 1)), mineig2(g2)))
