# Closure specification, rank 2: propagate gamma_2

**Status: PROPOSAL, not run.** Nothing here has been propagated. `tools/propagate_dynamical.py`
exists but is DEFECTIVE and was stopped before producing numbers -- see section 9.

Symbols as in `CLOSURE-SPEC-rank1.md`, plus:

| symbol | definition |
|---|---|
| `gamma_3[p,v,x,q,w,y]` | `<a+_p a+_v a+_x a_y a_w a_q>` |
| `lambda_3` | `gamma_3 - separable_3(gamma_1, lambda_2)`, the 3-cumulant |
| `separable_k` | every term of the MKNO expansion of `gamma_k` with `lambda_k` deleted; `qdiss.ansatz.separable` |
| `Mat_2[(p,v),(q,w)]` | `2 gamma_2[p,v,q,w]` over sorted pairs `p<v`, `q<w`; PSD for representable `gamma_2` |
| `Tr_3 T[pv,qw]` | `sum_y T[p v y, q w y]` |
| `anti(E,X)` | `E` contracted into each of `X`'s four slots and summed |
| `Delta` | `Tr_3 lambda_3 - anti(E, lambda_2)`, zero for the exact state |

## 1. Propagated state

**(a)** `R`, of shape `r x C(M,2)`, with

    X = Mat_2 = R^dag R          so gamma_2 >= 0 (the D condition) holds by construction

`gamma_1` is NOT propagated; it is recovered at every stage by contraction,
`gamma_1 = Tr_2(gamma_2)/(N-1)`, and `lambda_2 = gamma_2 - gamma_1 wedge gamma_1`. This is the
Lackner/Mueller scheme for the state variables.

**(b)** Coefficients `c_a`, `a = 1..A`.

## 2. Exact equation of motion being approximated

    d/dt gamma_2[pv,qw] = +i <[H, a+_p a+_v a_w a_q]>

generated in `qdiss.eom.wick` -- 4 one-body and 12 two-body terms, of which 8 carry `gamma_3` --
term-identical to `qdiss.eom.gwt` and validated against `qdiss.eom.fullci` at 1e-15. The term list is
used as generated; **no prefactor in it is chosen or corrected by hand.**

Structurally, in `Mat_2` variables: `Xdot = A X + X A^dag + Phi[X]`, where the first part is the
`gamma_2`-carrying drift (a Lyapunov map, verified to compress to zero on any kernel) and `Phi` is
the `gamma_3`-carrying collision term. All positivity risk is in `Phi`.

## 3. Integration of the state: the XR construction

`Xdot` from section 2 is lifted to `R` so that positivity cannot be lost. In the eigenbasis of `X`
with eigenvalues `x_i` and a FIXED `sigma > 0`:

    Gamma_ij = -2 Xdot_ij / (x_i + x_j + sigma) ,        Rdot = -(1/2) R Gamma

`Gamma` is Hermitian because `Xdot` is and the denominator is real symmetric. Justification: the exact
flow satisfies `Xdot = i(S^dag X - X S)` at every rank (derived in `XR-POSITIVITY.md`), and
`Rdot = -i R S` reproduces exactly that. Where `x_i + x_j` is healthy the closure's flow is
reproduced unchanged; where the eigenvalues vanish `Gamma` would diverge unless `Xdot_ij -> 0`, which
IS the tangency condition, so precisely the non-tangent components are suppressed by `sigma`.

## 4. The closure basis

`lambda_3 = sum_a c_a S_a`, with `S_a` from `qdiss.ansatz.connected_basis`: 12 generators, each an
antisymmetrised single contraction of two `lambda_2` factors joined through `I`, `gamma_1` or
`eta_1`, optionally dressed on the free legs by `anti(E, .)`, each pair-Hermitised. Then
orthonormalised by SVD to a frame `F` with orthonormal rows, and Procrustes-aligned to the previous
step's frame. `c_a = <F_a, lambda_3>`, `lambda_3 = F^T c`.

**Why this basis and not another.** Every term is a CONNECTED contraction, which is necessary and
sufficient for extensivity (`CONSTRAINT-LEDGER.md`); measured leak `1.8e-13` against `0.98-1.61` for
wedge alternatives. Every term is dimensionless. Leg counting forces exactly two `(2,2)` vertices
sharing one index, so no smaller connected structure exists.

**Numerical rank.** 12 generators span 4 to 6 directions depending on threshold (4 at `3e-5` on
singular values, 6 at `1e-9`). The frame rank must be FIXED across a step or the two frames have
mismatched shapes.

## 5. The coefficient flow

    c-dot  =  <F , lambda_3-dot>  -  <F , sum_b c_b S_b-dot>

with, at every stage:

    gamma_1     = Tr_2(gamma_2)/(N-1)
    lambda_2    = gamma_2 - gamma_1 wedge gamma_1
    gamma_3     = separable_3(gamma_1, lambda_2) + F^T c
    gamma_2-dot = rank-2 MKNO term list applied to (gamma_1, gamma_2, gamma_3)
    gamma_1-dot = Tr_2(gamma_2-dot)/(N-1)
    lambda_2-dot= gamma_2-dot - ( gamma_1-dot wedge gamma_1 + gamma_1 wedge gamma_1-dot )
    gamma_3-dot = rank-3 MKNO term list applied to (gamma_1, gamma_2, gamma_3, gamma_4 = TERMINATION)
    lambda_3-dot= gamma_3-dot - separable_dot(3, {gamma_1, lambda_2}, {gamma_1-dot, lambda_2-dot})
    S_b-dot     = qdiss.ansatz.connected_basis_dot(gamma_1, lambda_2, gamma_1-dot, lambda_2-dot)

`separable_dot` and `connected_basis_dot` are the product rule over the MKNO term lists; both
VALIDATED against finite differences at the truncation floor -- `9.1e-08` against `|dl3| = 0.0065`,
and worst relative structure-rate error `1.4e-05` (`tools/galerkin.py`).

**Termination.** `gamma_4 = 0`, which is EXACT when `N <= 3` and an approximation otherwise. Measured:
at `N = 3` it holds the trace identity at the noise floor; at `N = 4` it breaks it at `3.8e-2`. The
alternative `lambda_4 = 0` is WRONG at `N = 3` (breaks the identity at `2.7e-2`) and must not be used.

## 6. Every approximation, named

1. **Span truncation.** `lambda_3` confined to the connected span. Measured span error: 0.8% (`N=2`),
   2.9% (H3), 8.0% (H4) statically; 16% at `t=1` under a compressed->stretched quench; **83% from a
   rank-2 mixture.**
2. **Galerkin leakage.** The out-of-span residual rate is discarded. Measured 1-3% relative under
   quench (absolute `1.9e-2` against `|c-dot| = 0.56`), `3e-4` absolute adiabatically, 28-35% at low
   `rank(rho)`.
3. **Tower termination** `gamma_4 = 0` as above.
4. **`sigma` in the XR lift**, fixed; introduces an `O(sigma)` flow error concentrated in the
   small-eigenvalue directions, and an `O(sigma)` extensivity error since the eigenbasis of `X` is a
   global object.
5. **Frame rank fixed** across a step; if the true numerical rank changes, one direction is carried
   spuriously or dropped.
6. **`Delta` is not enforced.** See section 8.

## 7. Ledger status

| row | met? | evidence |
|---|---|---|
| `gamma_2 >= 0` (D) | YES, structural | XR construction |
| extensivity | YES, structural | connected basis, leak `1.8e-13` |
| dimensional analysis | YES | construction |
| Lipschitz | YES, fixed `sigma` | `tools/level_shift.py` |
| antisymmetry, hermiticity | YES | sorted-tuple representation |
| exactness at `k >= N` | YES for `N <= 3` via `gamma_4 = 0` | `tools/top_of_hierarchy.py` |
| trace identity | **NO** -- see section 8 | -- |
| `Q` and `G` cones | **NO**, uncontrolled | D+Q+G bound measured at `-1.945` vs `E_FCI = -1.943` |
| convergence | untested | -- |

## 8. The trace identity is NOT enforceable, correcting an earlier claim

I claimed that because `Delta = 0` is an invariant of the exact flow, it becomes "one linear condition
on `c-dot`, elementary to impose". **That is wrong.** With `D_a = Tr_3 S_a`,

    Delta = sum_a c_a D_a - anti(E, lambda_2)

is `C(M,2)^2` scalar conditions -- 784 at `M = 8` -- on `A ~ 6` unknowns. Differentiating gives the
same count on `c-dot`. So `Delta = 0` cannot be imposed exactly on a small coefficient space. Options:
least-squares imposition with a reported residual; enlarge `A` until `{D_a}` spans the needed
subspace (cost not estimated); or leave it violated and MONITOR `|Delta|` along the trajectory.

Separately, and this is what makes it a real defect rather than a formality: `d/dt Delta` depends on
`lambda_4`, so the `gamma_4 = 0` termination breaks the invariance at `3.8e-2` for `N = 4`. It will
drift.

## 9. What exists in code, and why it must not be run as is

`tools/propagate_dynamical.py` is DEFECTIVE and produced no numbers:

- it uses **Euler** for the coefficient leg while RK4 for the state;
- it computes a quantity `Fd` and then **never uses it**, so the basis-motion term
  `sum_b c_b S_b-dot` is **omitted entirely** -- and that term was measured at 60-96% of
  `|lambda_3-dot|`, so the `memory` row would have been meaningless;
- its docstring says `lambda_4 = 0` while the code uses `gamma_4 = 0`; the code is right and the
  docstring is wrong.

It should be rewritten against this specification, not patched.

## 10. Honest expectation

Given the measurements, the scheme should work in the near-adiabatic pure-state regime (span 3-8%,
leakage `3e-4` absolute) and should fail from a low-`rank(rho)` start (span 83%, leakage 30%). The
`Delta` drift is unbounded as specified. **The first run must reproduce a known limit** -- `c = 0`
against the already-measured `null` trajectory, and exact `gamma_3` fed in at every stage against the
exact trajectory -- before any new number is believed. This session produced four vacuous or
mis-normalised tests; a propagation has more places for that to hide.
