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

## 4. The closure, with the trace relation structural

**REVISED.** The earlier version of this file put `lambda_3` in a fitted connected span and left the
trace relation unenforced. That was wrong: fitting it left a 36% residual with common coefficients,
and imposing it on `c-dot` is overdetermined (784 conditions on ~6 unknowns). The relation is now
built in. Full derivation in `STRUCTURAL-TRACE.md`; the construction is

    T           =  anti(E, lambda_2)                                the trace target, 4 indices
    K_sigma     =  - lambda_2  o  ( gamma_1 eta_1 + sigma I )^{-1}  the (2,2) vertex, Tr_2 K = I + O(sigma)
    lambda_3    =  A_3[ T o K_sigma ]  +  sum_a c_a S_a^perp

with `A_3` the antisymmetriser over three upper and three lower indices, `o` a single shared
contraction, and `S_a^perp` connected structures satisfying `Tr_3 S_a^perp = 0`.

Leg counting forces this shape: `T` has 4 free indices, `lambda_3` needs 6, `gamma_1` vertices cannot
change the count, and a `(2,2)` vertex attached by one contraction gives exactly `+2`. So exactly one
`(2,2)` vertex is required -- not a choice.

**Consequences.**

- `Tr_3 lambda_3 = T` holds by construction to `O(sigma)`, with `sigma = 0` exact **subject to the
  VERIFY step of `STRUCTURAL-TRACE.md` section 5, which has NOT been run.**
- The coefficients `c_a` cannot disturb the trace relation, whatever they do, because
  `Tr_3 S_a^perp = 0`. The ancilla is confined to the traceless sector.
- Every term is one connected component, so the construction is extensive.
- The `a_r^{-1}` in `K` is bounded on the N-representable manifold by `D`/`Q` Cauchy-Schwarz
  (`T` and `lambda_2` each contribute `O(sqrt(a_r))`), so it is `O(1)`, not divergent. Whether that
  bound is UNIFORM is the open question -- not the pole itself.

**The traceless connected basis `S_a^perp` does not exist yet.** It is the projection of
`qdiss.ansatz.connected_basis` onto the null space of `Tr_3`, a linear operation, but it must be
built and its dimension measured. Until then `A = 0` and the scheme is the structural trace channel
alone.

## 5. The coefficient flow

Unchanged in form from the earlier version, but now acting only on the traceless sector:

    c-dot  =  <F^perp , lambda_3-dot>  -  <F^perp , sum_b c_b S_b^perp-dot>

with `F^perp` the orthonormalised, Procrustes-aligned frame of `S_a^perp`. At every stage:

    gamma_1     = Tr_2(gamma_2)/(N-1)
    lambda_2    = gamma_2 - gamma_1 wedge gamma_1
    T           = anti(E, lambda_2),   E = gamma_1 - I/2
    lambda_3    = A_3[T o K_sigma] + F^perp^T c
    gamma_3     = separable_3(gamma_1, lambda_2) + lambda_3
    gamma_2-dot = rank-2 MKNO term list on (gamma_1, gamma_2, gamma_3)
    gamma_1-dot = Tr_2(gamma_2-dot)/(N-1)
    lambda_2-dot= gamma_2-dot - ( gamma_1-dot wedge gamma_1 + gamma_1 wedge gamma_1-dot )
    gamma_3-dot = rank-3 MKNO term list on (gamma_1, gamma_2, gamma_3, gamma_4 = 0)
    lambda_3-dot= gamma_3-dot - separable_dot(3, {gamma_1,lambda_2}, {gamma_1-dot,lambda_2-dot})

`separable_dot` and `connected_basis_dot` are the product rule over the MKNO term lists, VALIDATED
against finite differences at the truncation floor (`9.1e-08` against `|dl3| = 0.0065`; worst
structure-rate error `1.4e-05`, `tools/galerkin.py`). The rate of the structural channel
`d/dt A_3[T o K_sigma]` needs the same treatment and **is not written yet.**

**Termination.** `gamma_4 = 0`: EXACT when `N <= 3`, an approximation otherwise. Measured -- at `N = 3`
it holds the trace identity at the noise floor, at `N = 4` it breaks it at `3.8e-2`. The alternative
`lambda_4 = 0` is WRONG at `N = 3` (breaks it at `2.7e-2`) and must not be used.

## 6. Every approximation, named

1. **`sigma` in `K_sigma`.** Makes `Tr_2 K = I + O(sigma)`, so the trace relation holds to `O(sigma)`
   rather than exactly. This is the ONLY place the trace relation is approximated, and the error is a
   dial, not a residual.
2. **Truncation of the traceless sector** to `A` structures. With `A = 0` the whole traceless content
   of `lambda_3` is dropped; the memory measured at 60-96% of `|lambda_3-dot|` under quench lives
   here, so `A = 0` should be expected to fail far from equilibrium.
3. **Galerkin leakage** in the traceless sector. Measured on the full connected span: 1-3% relative
   under quench, `3e-4` absolute adiabatically, 28-35% at low `rank(rho)`.
4. **Tower termination** `gamma_4 = 0` as above.
5. **`sigma` in the XR lift** (section 3), distinct from the `sigma` in `K`. Introduces `O(sigma)`
   flow error in small-eigenvalue directions and `O(sigma)` extensivity error, since the eigenbasis of
   `X` is a global object.
6. **Retraction / integrator error**, unmeasured.

## 7. Ledger status

| row | met? | how |
|---|---|---|
| `gamma_2 >= 0` (D condition) | YES, structural | `X = R^dag R`, XR construction |
| antisymmetry, hermiticity | YES, structural | sorted-tuple representation |
| extensivity | YES, structural | every term one connected component |
| dimensional analysis | YES | `K` dimensionless; no bare `V`, no amplitude |
| **trace relation** | **YES to `O(sigma)`, structural** -- pending VERIFY | `STRUCTURAL-TRACE.md` s.3-5 |
| Lipschitz | bounded by `D`/`Q` Cauchy-Schwarz; uniformity OPEN | `STRUCTURAL-TRACE.md` s.4 |
| exactness at `k >= N` | YES for `N <= 3` via `gamma_4 = 0` | `tools/top_of_hierarchy.py` |
| `Q` and `G` cones | **NO**, uncontrolled | D+Q+G bound `-1.945` vs `E_FCI = -1.943` |
| convergence | in `A` and tower depth; UNTESTED | -- |

Six of nine rows are now structural. The two that are not are `Q`/`G` and convergence.

## 8. What is still owed, in order

1. **VERIFY** the combinatorics of `Tr_3 A_3[T o K]` (`STRUCTURAL-TRACE.md` s.5). Cheap, and nothing
   downstream is meaningful without it.
2. **Build `S_a^perp`**, the traceless connected basis, and measure its dimension and how much of
   `lambda_3`'s traceless part it spans.
3. **Write `d/dt A_3[T o K_sigma]`** by the product rule, and validate against finite differences the
   way `separable_dot` was.
4. Only then propagate, and only after reproducing two known limits: `c = 0` with `K` switched off
   against the measured `null` trajectory, and exact `gamma_3` fed in at every stage against the exact
   trajectory.

## 9. Code status

`tools/propagate_dynamical.py` is DEFECTIVE and produced no numbers: Euler on the coefficient leg
against RK4 on the state; it computes `Fd` and never uses it, so the basis-motion term -- 60-96% of
`|lambda_3-dot|` -- is omitted entirely; and its docstring contradicts its own code on the termination.
It must be rewritten against this specification, not patched.

## 10. Honest expectation

With `A = 0` this is a one-shot closure with a structurally exact trace channel and nothing else. It
should be better than every closure measured so far on the trace relation (which none satisfied) and
no better than them on the traceless content (which it discards). Far from equilibrium it should fail,
because that is where the traceless memory lives. **It is a first rung, not a finished method.**
