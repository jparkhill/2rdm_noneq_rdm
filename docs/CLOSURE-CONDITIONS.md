# Conditions on cumulant reconstructions

Two tables. The first asks which of the conditions we have used at rank 2 are order-general;
the second inverts the question and asks, for each known physical property of the closed RDM
equation, what it actually demands of the reconstruction.

Ground truth throughout is `qdiss.eom.fullci`. Derivations are `qdiss.eom.wick` (plain vacuum
Wick) and `qdiss.eom.gwt` (MKNO), which produce identical term lists and are each validated
against FullCI at 1e-15 by `qdiss.eom.validate`.

## Table 1: order-generality of the conditions used at rank 2

| condition | rank-2 form | order-general? |
|---|---|---|
| antisymmetry, pair-hermiticity | structural | yes, trivially |
| partial-trace identity | `sum_x lam2[px,qx] = -(g1 - g1^2) = -g eta` | yes as a *relation*, no as a *form*: rank 3 gives `E lam2 + lam2 E` with `E = g1 - 1/2` |
| blocking | a `g*eta` PRODUCT, vanishing at both occupation boundaries | **no**. The rank-3 trace is an anticommutator with the IMBALANCE `(g-eta)/2`, which is maximal at the boundaries |
| exactness where `gamma_{k+1} = 0` | `lam2` vanishes at N=1 (2.8e-17) | in principle; but `lam3` does NOT vanish at N=2, so the naive version fails |
| reduction to a known limit | `lam2 = 0` gives TDHF exactly (3.9e-16) | yes: `lam_k = 0` gives the rank-(k-1) mean field |
| positivity of the reconstructed RDM | D, Q, G on `gamma_2` | yes, but the mechanism changes: at N <= 3 the conic condition degenerates to the linear `X_KK = 0` |
| phase relative to `V` | in-phase component cannot move occupations | presumably, untested above rank 2 |
| smoothness at the reference | `blocked` is 0/0 at Hartree-Fock | yes, but see the caveat below |

Rank 4 was computed (`tools/rank4_contraction.py`) to test whether the trace rule continues.
It does not: the naive analogue `E lam3 + lam3 E` fits with best scalar multiple c = 0.00,
leaving a residual larger than the target. Rank 4 is the first order carrying `lam2*lam2`
partitions with no rank-3 counterpart, so the general rule must include them.

Three of these conditions are weaker than they look, and should not be treated as requirements:

- **smoothness at the reference** presupposes that a stationary state of the dynamics exists,
  as Fermi-Dirac does for the `lam2 = 0` closure. No closure here has one yet, so the condition
  has nothing to be smooth at.
- **blocking** is not known to be a physical necessity. It appears on no row of Table 2.
- **phase relative to V** governs whether occupations can move and whether new poles appear.
  That is physics we WANT the hierarchy to deliver, but at higher order; the rank at which it
  first appears is what decides whether a given level of approximation is practical.

## Table 2: what each physical condition demands of the reconstruction

Measured by `tools/necessary_conditions.py` on H3 linear / sto-3g, M=6 N=3, at the exact FCI
ground state, where every quantity must vanish. The state is `(gamma_1, lambda_2)` and the
reconstructed object is `lambda_3`.

    lambda_3                |  herm g2d   antisym |  tr g1dot  tr g2dot |     dE/dt  contract  min X_KK
    exact                   |  8.47e-16  2.78e-16 |  1.18e-16  9.09e-17 |  0.00e+00  2.67e-16 -7.78e-17
    zero                    |  5.90e-16  3.34e-16 |  1.18e-16  2.57e-16 |  0.00e+00  2.21e-03 -3.28e-02
    -(A+B)  MK_gamma30      |  5.03e-16  1.87e-16 |  1.18e-16  6.66e-18 |  0.00e+00  5.25e-02 -8.72e-17
    closed-form trace-exact |  9.58e-16  2.98e-16 |  1.18e-16  4.86e-16 |  0.00e+00  1.79e-16 -3.39e-02

| physical condition | automatic or conditional | what it demands of `lambda_3` |
|---|---|---|
| hermiticity of `gamma_k-dot` | automatic | only that `lambda_3` be pair-hermitian |
| antisymmetry of `gamma_k-dot` | automatic | only that `lambda_3` be antisymmetric |
| particle number at rank 1, `tr gamma_1-dot = 0` | automatic | nothing |
| particle number at rank 2, `tr gamma_2-dot = 0` | automatic | nothing |
| energy conservation, `dE/dt = 0` | **automatic** | nothing |
| BBGKY contraction consistency | **conditional** | the partial-trace sum rule |
| positivity, Nagumo == CP at a point | **conditional** | the kernel-block condition on `X_KK` |
| exact ground state is a fixed point | **conditional, unmet** | nothing yet known; `\|lambda_2-dot\| = 1.37e-02` |

Energy conservation being automatic is provable, not accidental. `gamma_1-dot` contains no
`gamma_3`, so the whole `gamma_3` dependence of `E-dot = sum h gamma_1-dot + (1/4) sum V
gamma_2-dot` sits in the second term. Since `d<H>/dt = i<[H,H]> = 0` holds for arbitrary states
and therefore for arbitrary `gamma_3` at fixed `gamma_1, gamma_2`, the coefficient of `gamma_3`
must vanish identically. The measurement returns exactly 0.00e+00 rather than 1e-16,
consistent with a symbolic cancellation.

The consequence for design: the space of admissible reconstructions is far less constrained
than the list of physical desiderata suggests. Only two known conditions bite, and no single
reconstruction so far satisfies both -- `MK_gamma30` has positivity and fails the trace by
5.3e-2, the closed-form trace-exact reconstruction is the reverse. A `lambda_3` satisfying both
simultaneously was constructed by linear solve (`tools/lambda3_linear.py`) but is dynamically
unstable, which points at the third row as the real target.

## The stationary-state and dissipation thread

The condition that matters most and is least satisfied is the third: the exact ground state
should be a fixed point of the closed EOM. For `lam2 = 0` the fixed points are the Hartree-Fock
determinants, and Fermi-Dirac in the thermal case, which is why TDHF is usable. No cumulant
reconstruction here has an analogous stationary state.

A related conjecture worth keeping in view: if dissipators obeying detailed balance are
introduced in the higher-order EOMs, they may have a stabilising effect, plausibly by driving
the density toward dominance over the correlations -- the eigenstate-thermalization intuition
that dissipation should reduce many-body correlation. As more stable cumulant approximations
are found, whether they achieve stability by that mechanism is a question to test rather than
assume.
