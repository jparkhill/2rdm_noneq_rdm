#!/usr/bin/env python3
"""Poles and positivity versus bath strength, starting from the exact FCI density.

THE OBJECTIVE. Pole positions against FCI, and positivity of the reduced density matrices --
nothing else is ranked here. The question is whether a physical dephasing bath can be tuned so
that a closure both reproduces the FCI excitation spectrum and stays positive, and how wide that
window is.

THE APPROXIMATION. The propagated variable is gamma_2 alone; gamma_1 is recovered at every step as
Tr_2(gamma_2)/(N-1). The equation is the rank-2 Heisenberg EOM from qdiss.eom.wick (term-identical
to the MKNO route qdiss.eom.gwt, validated against qdiss.eom.fullci at 1e-15), closed by
reconstructing gamma_3:

    lam3_0      gamma_3 = the separable part of the MKNO cumulant expansion, i.e. lambda_3 = 0
    muller2405  Mueller et al. arXiv:2405.05093 Eq. 14 closed-system limit plus the Lackner
                Eq. 44/45 contraction correction
    tdhf        the CONTROL, and a different equation: the rank-1 MKNO EOM with lambda_2 := 0,
                propagating gamma_1. Included because its spectrum is known independently.

THE BATH. Pure dephasing in the occupation-number basis, jump operators L_r = a+_r a_r at rate g.
Under the adjoint dissipator an operator carrying charge q_r = d_rp + d_rv - d_rq - d_rw under
n_r decays at rate (g/2) sum_r q_r^2, so the induced term on gamma_2 is elementwise damping

    Dephase[gamma_2][p,v,q,w] = -(g/2) * S[p,v,q,w] * gamma_2[p,v,q,w],
    S[p,v,q,w] = sum_r (d_rp + d_rv - d_rq - d_rw)^2

which vanishes exactly when the index sets {p,v} and {q,w} coincide. Hence populations and
tr gamma_2 are untouched, and the contraction to gamma_1 is preserved: S[p,x,q,x] = 2 for p != q,
so gamma_1's off-diagonals damp at rate g, the correct one-body dephasing.

WHAT IS MEASURED, per closure and bath strength.
  poles       eigenvalues of the Jacobian of the flow at the exact FCI gamma_2, built by finite
              differences over a Hermitian, pair-antisymmetric basis. Imaginary parts are
              frequencies, real parts are decay rates the bath introduces. Frequencies are matched
              to the lowest FCI excitation energies by Hungarian assignment.
  positivity  a short propagation from the same point, reporting the smallest eigenvalue of the
              pair matrix Mat_2(gamma_2)[I,J] = 2 gamma_2[i1,i2,j1,j2] over sorted pairs, the
              occupation range, and the two trace errors.

STATUS 2026-08-09: THIS HARNESS RUNS BUT ITS FIRST RESULTS WERE NOT TRUSTWORTHY. Three defects,
all in the experiment design rather than the machinery, are documented here so the next run fixes
them rather than rediscovering them.

  1. THE KICK IS NOT N-REPRESENTABLE. `propagate` perturbs gamma_2 directly by a random
     antisymmetric Hermitian tensor. That starting point is not the 2-RDM of any state, and at
     (M,N) = (6,3) the exact gamma_2 sits EXACTLY on the PSD boundary (rank <= C(M,N-2) = 6 of
     15, a 9-dimensional kernel), so any perturbation goes negative immediately. The reported
     min eig D of -0.08 was the kick, not the closure -- visible in the fact that the positivity
     and trace columns were identical to five digits across every bath strength. FIX: perturb the
     WAVEFUNCTION with a dipole operator, exp(-i eps mu)|Psi>, and compute gamma_2 from it. That
     is both N-representable and physically the right probe for a spectrum.
  2. THE TEST SYSTEM HAS A DEGENERATE SPECTRUM. Asymmetric H3/sto-3g gave FCI excitations
     [0.0373, 0.0373, 0.1717, 0.1717, 0.1717] -- five targets with two distinct frequencies, so
     Hungarian matching is nearly trivial and the pole MAE is uninformative. This is the failure
     mode CLAUDE.md warns about from the linear H4 chain. H3 was chosen only because the M=8
     Jacobian did not finish in numpy. FIX: use the asymmetric trapezoid H4 and build the
     Jacobian on the GPU backend qdiss.gpu.
  3. POSITIVITY WAS REPORTED AGAINST ZERO. The meaningful comparison is against the EXACT
     clearance at the same state, which varies by five orders of magnitude between pure and mixed
     states. FIX: report min eig D alongside the exact clearance, or as a ratio.

WHAT DID SURVIVE. The TDHF row is interpretable because its spectrum is known independently: mean
width -0.0251 at g = 0, i.e. an intrinsically GROWING mode, turning to +0.0249 at g = 0.05 and
+0.5251 at g = 0.5, so the bath adds width almost exactly linearly in g and the marginal strength
is gamma* ~ 0.05 atomic units. The pole MAE stayed at 0.0139 across the whole sweep, so the bath
stabilised the flow without moving the frequencies. That is a real gamma* measurement in the rate
form, and it is the shape the corrected experiment should produce for every closure.

CAVEAT ON MATCHING, stated rather than buried: the assignment uses frequency only, with no
transition-density weighting, over a mode set far larger than the number of targets. That is
optimistic, and weaker than the biorthogonal scheme in tools/singlet_poles.py. The TDHF row is the
control for whether the procedure means anything at all.
"""
import itertools
import sys
import numpy as np
from scipy.optimize import linear_sum_assignment
sys.path.insert(0, "."); sys.path.insert(0, "tools")
from molecular_2rdm import build_H_sector, molecule, sector
from qdiss import closures, muller2405
from qdiss.eom import fullci, gwt, wick
from qdiss.symbolic import mkno

# Asymmetric H3, sto-3g: M = 6 spin-orbitals, N = 3 electrons (doublet). Chosen over H4 because
# the finite-difference Jacobian is C(M,2)^2 = 225 rather than 784 columns, each needing a full
# gamma_3 reconstruction; the H4 version did not finish. No spatial symmetry, so no degenerate
# poles. N = 3 respects the rule against testing a 2-RDM EOM on two electrons.
ATOM = "H 0.00 0.00 0.00; H 1.45 0.00 0.00; H 0.62 1.21 0.00"
_mol, _mf, M, h, V, enuc = molecule(ATOM, "sto-3g", 0, 1)
N = _mol.nelectron
EYE = np.eye(M, dtype=complex)
PAIRS = [(p, v) for p in range(M) for v in range(p + 1, M)]
NP = len(PAIRS)
ONE2, TWO2 = wick.terms(2)
CUM1 = gwt.terms_cumulant(1)

dets = sector(M, N)
es, vec = np.linalg.eigh(build_H_sector(h, V, M, dets))
psi0 = np.zeros(1 << M, dtype=complex); psi0[dets] = vec[:, 0]
g2_0 = fullci.rdm(psi0, M, 2)
g1_0 = np.einsum("pxqx->pq", g2_0) / (N - 1)
EXC = np.sort((es - es[0]).real)
EXC = EXC[EXC > 1e-8][:5]

S = np.zeros((M,) * 4)
for p, v, q, w in itertools.product(range(M), repeat=4):
    S[p, v, q, w] = sum((("" != "") + (r == p) + (r == v) - (r == q) - (r == w)) ** 2
                        for r in range(M))
S1 = np.zeros((M, M))
for p, q in itertools.product(range(M), repeat=2):
    S1[p, q] = sum(((r == p) - (r == q)) ** 2 for r in range(M))

SEP3 = [(s, ",".join("".join(u) + "".join(l) for n, r, u, l in f) + "->pvxqwy",
         [(n, r) for n, r, u, l in f])
        for s, f in mkno.gamma_to_cumulants(("p", "v", "x"), ("q", "w", "y"))
        if not any(n == "L" and r == 3 for n, r, u, l in f)]
T1 = [(s, [("".join(u) + "".join(l)) for nn, r, u, l in f], [(nn, r) for nn, r, u, l in f])
      for s, f in list(CUM1[0])]
T2 = [(s, [("".join(u) + "".join(l)) for nn, r, u, l in f], [(nn, r) for nn, r, u, l in f])
      for s, f in list(CUM1[1])]


def recon(nm, g1, g2):
    if nm == "muller2405":
        return muller2405.gamma3(g1, g2, N)
    l2 = g2 - closures.wedge(g1, g1)
    t = {("G", 1): g1, ("L", 2): l2}
    return sum(s * np.einsum(sub, *[t[k] for k in keys], optimize=True)
               for s, sub, keys in SEP3)


def rhs2(g2, nm, gb):
    g1 = np.einsum("pxqx->pq", g2) / (N - 1)
    d = wick.gamma_dot(2, h, V, {("D", 1): EYE, ("G", 1): g1, ("G", 2): g2,
                                 ("G", 3): recon(nm, g1, g2)}, M)
    return d - 0.5 * gb * S * g2


def rhs1(g1, gb):
    tab = {("D", 1): EYE, ("G", 1): g1, ("L", 2): np.zeros((M,) * 4, dtype=complex)}
    tot = np.zeros((M, M), dtype=complex)
    for s, subs, keys in T1:
        tot += s * np.einsum(",".join(["rs"] + subs) + "->pq", h,
                             *[tab[k] for k in keys], optimize=True)
    for s, subs, keys in T2:
        tot += s * 0.25 * np.einsum(",".join(["rtsu"] + subs) + "->pq", V,
                                    *[tab[k] for k in keys], optimize=True)
    return 1j * tot - 0.5 * gb * S1 * g1


def basis_pairs():
    B = []
    for a in range(NP):
        for b in range(a, NP):
            for ph in (1.0, 1j):
                Z = np.zeros((NP, NP), dtype=complex)
                Z[a, b] += ph; Z[b, a] += np.conj(ph)
                if np.abs(Z).max() < 1e-12:
                    continue
                T = np.zeros((M,) * 4, dtype=complex)
                for i, (p, v) in enumerate(PAIRS):
                    for j, (q, w) in enumerate(PAIRS):
                        x = 0.5 * Z[i, j]
                        T[p, v, q, w] = x; T[v, p, q, w] = -x
                        T[p, v, w, q] = -x; T[v, p, w, q] = x
                B.append(T)
    return B


def basis_one():
    B = []
    for i in range(M):
        E = np.zeros((M, M), dtype=complex); E[i, i] = 1.0; B.append(E)
    for i in range(M):
        for j in range(i + 1, M):
            E = np.zeros((M, M), dtype=complex); E[i, j] = 1.0; E[j, i] = 1.0; B.append(E)
            E = np.zeros((M, M), dtype=complex); E[i, j] = 1j; E[j, i] = -1j; B.append(E)
    return B


def jac(f, x0, B):
    n = len(B)
    nrm = [float(np.real(np.tensordot(b.conj(), b, axes=b.ndim))) for b in B]
    f0 = f(x0); J = np.zeros((n, n))
    eps = 1e-6
    for k in range(n):
        d = (f(x0 + eps * B[k]) - f0) / eps
        for m in range(n):
            J[m, k] = float(np.real(np.tensordot(B[m].conj(), d, axes=B[m].ndim))) / nrm[m]
    return J


def poles_of(J):
    w = np.linalg.eigvals(J)
    sel = w[w.imag > 1e-6]
    return np.abs(sel.imag), -sel.real


def match(om, gam):
    if len(om) == 0:
        return np.full(len(EXC), np.nan), np.nan, np.nan
    cost = np.abs(om[None, :] - EXC[:, None])
    r, c = linear_sum_assignment(cost)
    return om[c], np.mean(np.abs(om[c] - EXC)), np.mean(gam[c])


def pair_mat(T):
    A = np.array([[2.0 * T[p, v, q, w] for (q, w) in PAIRS] for (p, v) in PAIRS])
    return np.linalg.eigvalsh(0.5 * (A + A.conj().T)).min()


def propagate(nm, gb, dt=0.004, ns=100):
    rng = np.random.default_rng(5)
    k = rng.normal(size=(M,) * 4) + 1j * rng.normal(size=(M,) * 4)
    k = k - k.transpose(1, 0, 2, 3); k = k - k.transpose(0, 1, 3, 2)
    k = 0.5 * (k + k.transpose(2, 3, 0, 1).conj())
    y = g2_0 + 0.02 * k / np.abs(k).max()
    worst, occlo, occhi = 1e9, 1.0, 0.0
    for _ in range(ns):
        if not np.isfinite(y).all():
            return np.nan, np.nan, np.nan, np.nan, np.nan
        worst = min(worst, pair_mat(y))
        g1 = np.einsum("pxqx->pq", y) / (N - 1)
        o = np.linalg.eigvalsh(0.5 * (g1 + g1.conj().T))
        occlo, occhi = min(occlo, o.min()), max(occhi, o.max())
        k1 = rhs2(y, nm, gb); k2 = rhs2(y + dt / 2 * k1, nm, gb)
        k3 = rhs2(y + dt / 2 * k2, nm, gb); k4 = rhs2(y + dt * k3, nm, gb)
        y = y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    g1 = np.einsum("pxqx->pq", y) / (N - 1)
    return (worst, occlo, occhi, abs(np.trace(g1).real - N),
            abs(np.einsum("pvpv->", y).real - N * (N - 1)))


BP, BO = basis_pairs(), basis_one()
print("asymmetric H3 / sto-3g, M=%d N=%d;  FCI excitations %s"
      % (M, N, np.array2string(EXC, precision=4)))
print("bath: pure dephasing in the occupation basis, rate g (atomic units)")
print()
print("  %-11s %6s | %9s %9s | %10s %8s %8s | %9s %9s"
      % ("closure", "g", "pole MAE", "mean width", "min eig D", "min n", "max n",
         "|tr g1-N|", "|tr g2|err"))
for nm in ("tdhf", "lam3_0", "muller2405"):
    for gb in (0.0, 0.05, 0.2, 0.5):
        if nm == "tdhf":
            J = jac(lambda x: rhs1(x, gb), g1_0, BO)
        else:
            J = jac(lambda x: rhs2(x, nm, gb), g2_0, BP)
        om, gam = poles_of(J)
        _, mae, wid = match(om, gam)
        if nm == "tdhf":
            print("  %-11s %6.2f | %9.4f %9.4f | %10s %8s %8s | %9s %9s"
                  % (nm, gb, mae, wid, "-", "-", "-", "-", "-"))
        else:
            mn, lo, hi, e1, e2 = propagate(nm, gb)
            print("  %-11s %6.2f | %9.4f %9.4f | %10.5f %8.4f %8.4f | %9.2e %9.2e"
                  % (nm, gb, mae, wid, mn, lo, hi, e1, e2))
    print()
