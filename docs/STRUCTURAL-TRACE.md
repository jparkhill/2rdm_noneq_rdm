# Making the trace relation structural

**Status: DERIVATION with a named verification step. Nothing here is run.** The one bookkeeping step
that needs symbolic checking before use is marked VERIFY in section 5.

## 1. What went wrong before

The trace relation was previously handled three ways, all wrong:

| attempt | what happened |
|---|---|
| minimum-norm partial-trace inverse (`closures.Corr`, Lackner Eq. 43) | satisfies it exactly, but by WEDGING with `I` -- disconnected, so leak 0.73-0.88 |
| fit coefficients of connected structures to it | 36% residual with coefficients common across five states; not representable |
| penalise / project `c-dot` | `C(M,2)^2 = 784` conditions at `M=8` on `A ~ 6` unknowns; overdetermined, cannot be imposed |

The mistake in the second and third is the same: I treated the relation as a CONSTRAINT to be
satisfied by an object built for other reasons, instead of BUILDING the object out of the relation.

## 2. The leg-counting fact that forces the construction

Write `(u,l)` for free upper and lower index counts. Contractions remove one of each.

- `gamma_1` is `(1,1)`. A connected graph of `(1,1)` vertices is a path or a cycle, so it is `(1,1)`
  or `(0,0)`: **`gamma_1` vertices cannot increase the total leg count.**
- A `(2,2)` vertex attached by ONE contraction adds 4 legs and consumes 2: **net `+2`.**

The trace target at rank `k+1` is

    T_{k+1}  =  anti(E, lambda_k) + chains          E = gamma_1 - I/2

which has `2k` free indices, while `lambda_{k+1}` needs `2(k+1)`. So exactly ONE `(2,2)` vertex,
attached by one contraction, is required -- no more, no fewer. This is forced, not chosen.

## 3. The construction

Let `K` be the `(2,2)` vertex. Build

    lambda_{k+1}  =  A_{k+1} [  sum_r  T_{k+1}[... r ...]  K[r ... , ...]  ]

with `A_{k+1}` the full antisymmetriser over the `k+1` upper and `k+1` lower indices. Now take the
partial trace. `Tr_{k+1}` closes one upper against one lower and sums. The antisymmetriser
distributes the traced pair over the slots, giving three kinds of term:

  (i)  **both traced indices on `K`** -> `T_{k+1}` contracted with `Tr_2 K`, a `2k`-index object;
  (ii) **both on `T_{k+1}`** -> `(Tr_k T_{k+1})` contracted with `K`;
  (iii)**one on each** -> `T_{k+1}` and `K` contracted twice.

Term (i) returns `T_{k+1}` itself, up to whatever `Tr_2 K` is. So:

> **If `Tr_2 K = I` in the contracted slot, and if the coefficients of (ii) and (iii) cancel against
> the combinatorial factor of (i), then `Tr_{k+1} lambda_{k+1} = T_{k+1}` BY CONSTRUCTION.**

That is the structural enforcement. The trace relation is not imposed on the answer; the answer is
assembled from the target.

## 4. Fixing `K`

`K` must be `(2,2)`, connected as a vertex, dimensionless, extensive, and normalised so `Tr_2 K = I`.

The only `(2,2)` primitives available are `lambda_2` and `V`. `lambda_2` has a known partial trace,

    Tr_2 lambda_2 = -(gamma_1 - gamma_1^2) = -gamma_1 eta_1

so the normalised vertex is

    K  =  - lambda_2  o  (gamma_1 eta_1)^{-1}        (inverse on the contracted slot only)

which has `Tr_2 K = I` by construction. **This reintroduces `a_r^{-1}` with `a_r = n_r(1-n_r)`, and
that is now expected rather than fatal**: on the N-representable manifold the `D` and `Q` Gram
conditions bound every `lambda_2` element carrying index `r` by `O(sqrt(a_r))`, and
`T_{k+1} = anti(E,lambda_k)` carries a further `O(sqrt(a_r))` through its own `lambda_k`, so

    T  o  a_r^{-1}  o  lambda_2   =   O(sqrt(a_r)) . a_r^{-1} . O(sqrt(a_r))   =   O(1)

bounded. The regularised form with a FIXED shift, used in practice,

    K_sigma  =  - lambda_2  o  ( gamma_1 eta_1 + sigma I )^{-1}

satisfies `Tr_2 K_sigma = I + O(sigma)`, so the trace relation holds to `O(sigma)` rather than
exactly. **The exact-trace version is `sigma = 0` and is admissible only if the boundedness above is
uniform; that is the open question, not the pole itself.**

## 5. What must be verified before this is used

**VERIFY (symbolic, finite, cheap).** The combinatorial coefficients of terms (i), (ii) and (iii) in
section 3, for the specific antisymmetriser `A_{k+1}` and the specific contraction pattern. Concretely:
build `lambda_3 = A_3[T o K]` from exact `gamma_1, lambda_2` at several states, compute `Tr_3` of it,
and compare against `T = anti(E, lambda_2)`. If the residual is at round-off the construction is
correct as written; if it is a clean rational multiple, the normalisation of `K` absorbs it; if it is
state-dependent, terms (ii) and (iii) do not cancel and `K` must be augmented.

Until that check is run, the claim "trace-structural" is a DERIVATION, not a fact. I am not going to
assert it works, because the last three times I asserted a trace result it was wrong.

## 6. Extensivity of the construction

Every term is a single connected contraction: `T` and `K` share exactly one summed index, and `K`'s
own normaliser `(gamma_1 eta_1 + sigma)^{-1}` is a `(1,1)` factor on that same line, so the whole
diagram is one component. By the theorem in `CONSTRAINT-LEDGER.md` (additive under direct sum iff the
contraction graph is connected), the construction is extensive.

Note the contrast with the wedge inverse: that one attaches `I` as a SEPARATE component, which is why
it leaks. Here the added `(2,2)` vertex is attached, not appended.

## 7. What this removes from the earlier specifications

- No coefficient fitting for the trace relation.
- No penalty term, no `Pi_Delta` projection, no least-squares residual.
- No MP2 amplitude. The `(2,2)` vertex is `lambda_2`, a propagated state variable, not a perturbative
  object. See `CLOSURE-SPEC-rank1.md` section 4 for why the rank-1 case is different and harder.

The coefficient/ancilla machinery is NOT removed -- it remains the mechanism for the memory content
measured at 60-96% of `|lambda_3-dot|` under a quench. But it now applies to the part of
`lambda_{k+1}` ORTHOGONAL to the trace channel, with the trace channel fixed structurally:

    lambda_{k+1}  =  A_{k+1}[ T_{k+1} o K_sigma ]   +   sum_a c_a S_a^perp

where `S_a^perp` are connected structures with `Tr_{k+1} S_a^perp = 0`, so the coefficients cannot
disturb the trace relation no matter what they do. **Constructing that traceless connected basis is
the second piece of algebra still owed** and is a projection of the existing `connected_basis` onto
the null space of `Tr_{k+1}`, which is a linear operation.
