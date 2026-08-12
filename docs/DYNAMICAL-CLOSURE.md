# The dynamical closure: cumulants as ancillae, not functionals

John's reframing, 2026-08-12. It supersedes the search for a one-shot reconstruction and gives the
programme the thing the ledger never had: **a hierarchy of approximations that converges to the exact
equation of motion in polynomial compute.**

## The error in every closure attempted so far, including all of ours

The exact hierarchy is

    lambda_2-dot = f(gamma_1, lambda_2, lambda_3)
    lambda_3-dot = g(gamma_1, lambda_2, lambda_3, lambda_4)

so `lambda_3` is an INDEPENDENT DYNAMICAL OBJECT. Writing `lambda_3 = F[gamma_1, lambda_2]` asserts
that it instantaneously follows the lower ranks. Every closure in this programme -- `null`, `trace`,
Mueller, the connected chains, the ladder -- is therefore the ZEROTH TERM OF A MEMORY EXPANSION.
That is why `pt` reproduced MP1/MP2 exactly at an idempotent density and was no better than
`lambda_2 = 0` away from it: the adiabatic limit is only valid near a stationary point. And it is why
the closures failed worst at low `rank(rho)` (`Neff ~ 2-3`), which is precisely where no timescale
separation exists.

## What functional is being expanded

The Mori-Zwanzig memory kernel of the projection of the exact hierarchy onto the propagated
variables. With `P` the projection onto `(gamma_1, lambda_2)`,

    lambda_3(t) = [homogeneous term from the initial condition]
                + integral_0^t K(t-s)[lambda_2(s)] ds

Zeroth order, `K -> K_0 delta(t-s)`, is the Markovian limit: the one-shot reconstruction. Higher
orders are the temporal structure of the kernel, carried by auxiliary variables.

**The expansion parameter is NOT V.** It is the timescale separation between the propagated and the
closed-off ranks. This matters because it explains the observed failures without leaning on
perturbation theory, which John has repeatedly and correctly refused as a validity criterion.

## The construction

The connected span at rank 3 is small -- 4-dimensional -- and was measured to CONTAIN 99.4% of
`lambda_3` (tools/connected_chains.py). So do not reconstruct `lambda_3`; propagate coefficients in
that span:

    lambda_3(t) = sum_a c_a(t) S_a[gamma_1(t), lambda_2(t), V]

with `S_a` the connected, dimensionally admissible structures: `lambda_2 o lambda_2`,
`lambda_2 o T`, `T o T` with `gamma_1` dressings, where `T = V / DeltaE` is an amplitude carrying a
FIXED level shift `sigma ~ 1e-2` (measured to cost nothing in accuracy while passing the Lipschitz
screen). The coefficients obey a Galerkin projection of the MKNO rank-3 equation of motion:

    lambda_3-dot = sum_a c_a-dot S_a + sum_a c_a S_a-dot        (S_a moves through gamma_1, lambda_2)
    c-dot = Gram^{-1} <S, g(gamma_1, lambda_2, lambda_3, lambda_4) - sum_a c_a S_a-dot>

a small linear solve per step. This is HEOM in structure: `S_a` is the auxiliary basis, `c_a` the
auxiliary amplitudes, and the tower closes at finite depth.

**Positivity is not the ancilla's problem.** `gamma_2 = R^dag R` in the XR representation is positive
by construction whatever `lambda_3` does, and `lambda_3` is not a density matrix so it need not be
positive. This is exactly why the scheme escapes the objection that time-dependent perturbation
theory is non-positive.

## Measured: the trace identity is conserved, but NOT by truncation

With `Delta = Tr_3 lambda_3 - anti(E, lambda_2)`, tools/trace_conservation.py at exact
non-stationary FCI states:

    === H3 bent  M=6 N=3 ===                      (gamma_4 == 0 exactly here)
        |Delta| at t=0                               2.776e-16
        d/dt Delta, EXACT rate                       1.550e-08     <- finite-difference noise floor
        d/dt Delta, gamma_4 exact                    1.550e-08     conserved
        d/dt Delta, gamma_4 = 0                      1.550e-08     (identical: exact at N=3)
        d/dt Delta, lambda_4 = 0                     2.662e-02     BREAKS

    === H4 trapezoid  M=8 N=4 ===
        |Delta| at t=0                               7.017e-16
        d/dt Delta, EXACT rate                       4.328e-08     <- noise floor
        d/dt Delta, gamma_4 exact                    4.328e-08     conserved
        d/dt Delta, gamma_4 = 0                      3.830e-02     BREAKS
        d/dt Delta, lambda_4 = 0                     1.016e-02     BREAKS

So the identity IS an invariant of the exact flow -- `d/dt Delta` at the noise floor against
`|d gamma_3/dt| = 0.18` -- and `d/dt Delta` DOES depend on `lambda_4`, so every approximate
termination breaks it at the 1-4e-2 level. It must therefore be imposed on the coefficient flow.

**That is the good outcome.** Compare where the identity can be imposed:

| where | what it demands | verdict |
|---|---|---|
| on a reconstruction `lambda_3 = F[gamma_1,lambda_2]` | a functional satisfying it POINTWISE | IMPOSSIBLE: 36% residual with common coefficients over five states, `tools/chain_trace_closure.py` |
| on the auxiliary flow `c-dot` | one LINEAR condition on a handful of numbers | elementary: project `c-dot` onto the tangent space of `Delta = 0` |

The identity moves from an unsatisfiable constraint on a functional form to a trivially satisfiable
constraint on an ODE. The obstruction that consumed several turns of work was an artifact of
insisting on a one-shot formula.

## Where every ledger row lands

| row | home in the dynamical picture | cost |
|---|---|---|
| trace identity | linear projection on `c-dot` | free |
| positivity of `gamma_2` | XR representation, `X = R^dag R` | free |
| extensivity | connectedness of the `S_a` | free |
| dimensional analysis | `S_a` construction; forces `T = V/DeltaE` not bare `V` | free |
| Lipschitz | fixed level shift in `T` | measured to cost nothing |
| exactness at `k >= N` | the tower terminates there exactly | free |
| CONVERGENCE | number of auxiliaries `A`, and tower depth | the hierarchy that was missing |

## Next

1. The Galerkin machinery: `<S_a, S_b>`, `<S_a, g(...)>`, and `S_a-dot`, validated term-by-term
   against `qdiss.eom.wick` / `fullci` before use.
2. Integrate `(R, c)` jointly and compute POLES AGAINST FCI -- the programme's primary metric, never
   yet computed for any new construction.
3. Convergence study in `A` and tower depth, including a mixed-state row and the low-`rank(rho)`
   regime (`Neff ~ 2-3`) where the one-shot closures failed worst and where the memory expansion
   should help most.
