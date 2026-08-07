#!/usr/bin/env python3
"""Stationary densities for cumulant closures: energy minimised subject to EOM stationarity and
the Pauli bound, the latter imposed through power-iteration estimates of the extremal occupations.

WHY THIS LAGRANGIAN AND NOT THE HF ONE. The HF variational problem constrains gamma_1 to be
idempotent, and a real HF solver's multiplier enforces exactly that (as orbital orthonormality).
That is the wrong constraint for a cumulant closure. The rank-2 partial-trace identity
    sum_x lambda_2[px,qx] = -(gamma_1 - gamma_1^2)
forces the cumulant's contraction to vanish whenever gamma_1 is idempotent, so every closure
obeying the identity collapses to Hartree-Fock at a determinant; for `blocked`, whose whole tensor
is built from A = gamma_1 - gamma_1^2, the cumulant itself vanishes, verified numerically
(|lambda_2| falls linearly with |gamma_1^2 - gamma_1|). The content of a cumulant reconstruction
lives entirely in the FRACTIONAL occupations, so the constraint must be the INEQUALITY
0 <= gamma_1 <= 1, not the equality gamma_1^2 = gamma_1.

THE LAGRANGIAN. Minimise E subject to one matrix equality, one scalar equality and two scalar
inequalities:

    L = E[g] - Tr(Lambda R[g]) - mu (tr g - N)
          + Psi(n_max - 1, nu_hi, rho) + Psi(-n_min, nu_lo, rho)

    Psi(c, nu, rho) = (1/(2 rho)) [ max(0, nu + rho c)^2 - nu^2 ]

Psi is the Powell-Hestenes-Rockafellar augmented-Lagrangian term for an inequality c <= 0. It is
C^1, it carries a genuine multiplier nu >= 0, and stationarity in nu reproduces the KKT
complementarity condition nu c = 0 exactly: dPsi/dnu = (1/rho)[max(0, nu + rho c) - nu], which is
c when the constraint is active and -nu/rho when it is slack, so ascent in nu with nu clamped at
zero IS the standard multiplier update nu <- max(0, nu + rho c). Equalities keep plain multipliers.

THE PAULI BOUND AS ONE MATRIX-POWER TRACE, not eigh and not power iteration. Since
0 <= n <= 1 is exactly |2n - 1| <= 1, the bound is the spectral-radius condition
||2g - I|| <= 1, for which
    c(g) = ( tr[(2g - I)^(2k)] / M )^(1/(2k)) - 1 <= 0
a normalised POWER MEAN of the |2n - 1|: it needs only k matrix squarings, autograd differentiates
it exactly, and it increases monotonically to the exact spectral radius as k -> infinity. It sits
at equality on any determinant, where every |2n - 1| = 1 and the trace is exactly M, so the
boundary it defines is the right one, and it is strictly slack for fractional occupations. The
2k-th root is essential, not cosmetic: the bare trace tr[(2g-I)^(2k)] - M has the same zero set but
gradients scaling as 2k c^(2k-1), which at fixed step made |dL/dg| reach 1e4 and the iteration
useless. The root restores O(1) gradients. Being a mean it is necessary rather than sufficient at
finite k -- one large violation could in principle be masked by fractional occupations elsewhere --
so the true spectrum is checked by eigh afterwards, as a diagnostic outside the graph.

Shift-and-power-iterate was tried first and abandoned: with the shift needed to make the operator
positive definite the spectrum compresses to [2.01, 2.99], consecutive ratios reach 0.999, and the
iteration does not converge in any affordable number of steps -- it reported n_max = 0.9783 where
eigh gives 0.9877, so the constraint never activated.

APPROXIMATION per row. Solved variable is gamma_1 alone. The equation whose stationary points are
sought is the MKNO rank-1 EOM, qdiss.eom.gwt.terms_cumulant(1), a function of gamma_1 and
lambda_2 alone, validated against exact Fock-space algebra qdiss.eom.fullci at 1e-15. Energy is
the matching functional E = h.g + (1/4) V.(g^g + lambda_2).
    lam2_0   two-particle cumulant TRUNCATED, lambda_2 := 0. This is TDHF, and the CONTROL: with
             the Pauli bound active the answer must be the HF density, known from PySCF SCF.
    blocked  ours: lambda_2 RECONSTRUCTED as A^f + f^A with A = gamma_1 - gamma_1^2 and f fixed
             by demanding the rank-2 partial-trace identity exactly.
Solver: primal-dual residual iteration, descent in gamma_1 and ascent in every multiplier because
L is a saddle, alpha = 1e-3, DIIS on the true residual. Started from the FCI 1-RDM.
"""
import sys
import numpy as np
import torch
sys.path.insert(0, "."); sys.path.insert(0, "tools")
import stationary_lagrangian as S

M, N, enuc, dt64 = S.M, S.N, S.enuc, S.dt64
nB, herm, coeffs = S.nB, S.herm, S.coeffs
g_hf, g_fci, E_HF, E_FCI = S.g_hf, S.g_fci, S.E_HF, S.E_FCI
EYEt = S.EYEt
KPOW = 6                       # exponent p = 2^KPOW = 64
RHO = 20.0
ALPHA_NU = 0.5                 # step for nu, larger than the primal step


def box_violation(G, kpow=KPOW):
    """(tr[(2G-I)^p]/M)^(1/p) - 1 with p = 2^kpow, the power mean surrogate for 0 <= G <= 1.

    kpow squarings give exponent p = 2^kpow, so the root MUST be 1/p, not 1/(2 kpow): with the
    root mismatched the quantity is a power of the true power mean, which has the same zero set
    but is far blunter, and at p = 16 it scored -3.6e-05 (feasible) on a density whose true
    occupations were 1.0124 and -0.0225. p = 64 separates that case decisively.
    """
    Y = 2.0 * G - EYEt
    # Rescale before the power: tr[Y^p] = s^p tr[(Y/s)^p], so the surrogate is exactly
    # s (tr[(Y/s)^p]/M)^(1/p) - 1 with Y/s of spectral radius <= 1. Without this, a runaway
    # iterate overflows float64 at p = 256 and the whole solve goes to NaN.
    s = torch.clamp(torch.linalg.matrix_norm(Y), min=1.0)
    Z = (Y / s) @ (Y / s)
    for _ in range(kpow - 1):
        Z = Z @ Z
    return s * (torch.real(torch.einsum("pp->", Z)) / M) ** (2.0 ** -kpow) - 1.0


def psi(c, nu, rho=RHO):
    """PHR augmented-Lagrangian term for the inequality c <= 0."""
    return (torch.clamp(nu + rho * c, min=0.0) ** 2 - nu ** 2) / (2.0 * rho)


def solve(closure, x0, use_R=True, use_box=True, niter=4000, alpha=1.0e-3, label=""):
    xg = x0.clone().requires_grad_(True)
    xL = torch.zeros(nB, dtype=torch.float64, requires_grad=True)
    mu = torch.zeros(1, dtype=torch.float64, requires_grad=True)
    nu = torch.zeros(1, dtype=torch.float64, requires_grad=True)
    var = [xg] + ([xL] if use_R else []) + [mu] + ([nu] if use_box else [])
    # ascent in every multiplier block, descent in gamma_1: L is a min over gamma, max over duals
    sgn = torch.cat([torch.ones(xg.numel(), dtype=torch.float64)]
                    + [-torch.ones(v.numel(), dtype=torch.float64) for v in var[1:]])
    # per-block step: primal alpha for gamma and the equality duals, a much larger one for nu,
    # because the PHR multiplier update nu <- max(0, nu + rho c) needs O(1) rate to activate.
    stp = torch.cat([torch.full((v.numel(),), alpha, dtype=torch.float64) for v in var])
    if use_box:
        stp[-1] = ALPHA_NU
    inu = sum(v.numel() for v in var[:-1]) if use_box else None
    hist_x, hist_r = [], []
    print("=== %-8s %-18s ===" % (closure, label))
    print("   %5s | %10s %10s %9s | %8s %8s %9s %8s | %12s"
          % ("iter", "|dL/dg|", "|R|", "|tr-N|", "n_max", "n_min", "box c", "nu", "E - E_HF"))
    for it in range(niter):
        G = herm(xg)
        Ls = S.energy_t(G, closure) - mu[0] * (torch.real(torch.einsum("pp->", G)) - N)
        if use_R:
            Ls = Ls - torch.real((herm(xL).conj() * S.rhs_t(G, closure)).sum())
        if use_box:
            Ls = Ls + psi(box_violation(G), nu[0])
        grads = torch.autograd.grad(Ls, var)
        r = torch.cat([g.detach().reshape(-1) for g in grads])
        x = torch.cat([v.detach().reshape(-1) for v in var])
        hist_x.append(x.clone()); hist_r.append(r.clone())
        if len(hist_x) > 10:
            hist_x.pop(0); hist_r.pop(0)
        nh = len(hist_x)
        cand = x - stp * (sgn * r)
        if nh >= 2:
            B = torch.zeros(nh + 1, nh + 1, dtype=torch.float64)
            for a in range(nh):
                for b in range(nh):
                    B[a, b] = torch.dot(hist_r[a], hist_r[b])
            B[nh, :nh] = 1.0; B[:nh, nh] = 1.0
            rv = torch.zeros(nh + 1, dtype=torch.float64); rv[nh] = 1.0
            try:
                c = torch.linalg.solve(B + 1e-10 * torch.eye(nh + 1, dtype=torch.float64), rv)
                cand = (sum(c[a] * hist_x[a] for a in range(nh))
                        - stp * (sgn * sum(c[a] * hist_r[a] for a in range(nh))))
            except Exception:
                pass
        if use_box:                                  # nu >= 0 is part of the KKT conditions
            cand[inu:inu + 1] = torch.clamp(cand[inu:inu + 1], min=0.0)
        with torch.no_grad():
            o = 0
            for v in var:
                v.copy_(cand[o:o + v.numel()].reshape(v.shape)); o += v.numel()
        if it % 1000 == 0 or it == niter - 1:
            with torch.no_grad():
                Gn = herm(xg)
                if not torch.isfinite(Gn).all():
                    print("   %5d | diverged to non-finite" % it); break
                ex = np.linalg.eigvalsh(0.5 * (Gn.numpy() + Gn.numpy().conj().T))
                print("   %5d | %10.3e %10.3e %9.2e | %8.4f %8.4f %9.2e %8.2e | %12.3e"
                      % (it, grads[0].norm(), coeffs(S.rhs_t(Gn, closure)).norm(),
                         abs(float(torch.real(torch.einsum("pp->", Gn))) - N),
                         ex.max(), ex.min(), float(box_violation(Gn)),
                         float(nu[0]) if use_box else 0.0,
                         float(S.energy_t(Gn, closure)) + enuc - E_HF))
    with torch.no_grad():
        Gn = herm(xg).numpy()
        occ = np.sort(np.linalg.eigvalsh(0.5 * (Gn + Gn.conj().T)))[::-1]
        l2 = (S.lam2_blocked_t(herm(xg)).numpy() if closure != "lam2_0"
              else np.zeros((M,) * 4, dtype=complex))
        tr_res = np.abs(np.einsum("pxqx->pq", l2) + (Gn - Gn @ Gn)).max()
        print("   occupations: %s" % np.array2string(occ, precision=5, suppress_small=True))
        print("   |g-g_HF| %.3e  |g-g_FCI| %.3e  |lambda_2| %.3e  trace-identity res %.2e"
              % (np.abs(Gn - g_hf).max(), np.abs(Gn - g_fci).max(),
                 np.abs(l2).max(), tr_res))
    print()
    return Gn


print("references: E_HF = %.9f   E_FCI = %.9f   E_corr = %.9f"
      % (E_HF, E_FCI, E_FCI - E_HF))
print("            box surrogate (tr[(2g-I)^%d]/M)^(1/%d)-1: at g_HF %+.3e (must be 0, determinant)"
      % (2 ** KPOW, 2 ** KPOW, float(box_violation(torch.tensor(g_hf, dtype=dt64)))))
print("                                     at g_FCI %+.3e (must be < 0, fractional)"
      % float(box_violation(torch.tensor(g_fci, dtype=dt64))))
_prev = np.diag(np.array([1.01245, 1.01016, .97124, .96817, .04176, .03838, -.01956, -.02255]))
print("                                     at the p=16 run's own answer %+.3e (must now be > 0)"
      % float(box_violation(torch.tensor(_prev, dtype=dt64))))
_bad = g_fci + 0.10 * (np.diag([1.0] + [0.0] * (M - 2) + [-1.0]).astype(complex))
print("                                     at g_FCI kicked out of range %+.3e (>0)"
      % float(box_violation(torch.tensor(_bad, dtype=dt64))))
print()

x_fci = coeffs(torch.tensor(g_fci, dtype=dt64))
solve("lam2_0", x_fci, True, True, label="EOM + Pauli box")     # gate: must return HF
solve("blocked", x_fci, True, True, label="EOM + Pauli box")
solve("blocked", x_fci, False, True, label="Pauli box only")
