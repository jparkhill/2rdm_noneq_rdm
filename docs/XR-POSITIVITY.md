# The XR positivity construction

> **Attribution note (2026-08-13).** Internal name. Ingredients (RDMs as Gram matrices, GKSL-form tangency) are
> standard; whether the rank-general refactoring below appears in Mazziotti (ACSE), Lackner et al. or Joost et al.
> is UNCHECKED. See docs/PRIOR-ART-TODO.md item 2 before citing this as new.

Named for future reference: **the XR construction**, `X = R^dag R` with `Rdot = -i R S`. It makes
positivity of the rank-k reduced density matrix structural in the REPRESENTATION rather than a
constraint on the closure.

## Symbols

For a pair index `I = (p,v)` with `p < v`, `A_I = a_v a_p`, and `zeta_I = A_I |Psi>`. Then
`X[I,J] = <zeta_I, zeta_J>` is `Mat_2` up to convention, and is positive semidefinite because it
is a Gram matrix. `chi_J = [A_J, H]|Psi>`. Inner product `<u,v> = u^dag v`. At rank k, `A_I` is
the k-fold annihilator over a sorted k-tuple and `X = Mat_k`.

## The structure theorem

Differentiating the Gram matrix with `|Psidot> = -i H |Psi>` gives
`zetadot_I = -i H zeta_I - i chi_I`, and the two `H`-sandwich terms cancel identically:

    Xdot = -i ( B - B^dag ),        B[I,J] = <zeta_I, chi_J>

Expand `chi_J = [A_J, H]|Psi>`. The one-body part of `H` returns a linear combination of pair
annihilators, contributing `B_1 = X F^T`. The two-body part, after anticommuting, gives terms with
ONE creator and THREE annihilators, and those refactor as

    a+_t a_w a_u a_s |Psi>  =  (a+_t a_w) · (a_u a_s |Psi>)  =  Omega · zeta_K,    K = (s,u)

with `Omega = a+_t a_w` a ONE-BODY operator. Since `zeta_I` lies in `span{zeta}`, the inner product
sees only the projection, so writing `P_zeta Omega zeta_K = sum_K' m_KK' zeta_K'`,

    B_3[I,J] = sum V(...) <zeta_I, Omega zeta_K> = sum V(...) m_KK' X[I,K'] = (X R)[I,J]

Both parts factor through `X`. Hence with `S = R + F^T`:

    Xdot  =  i ( S^dag X - X S )     EXACTLY, for the true dynamics.

ORDER-GENERAL. At rank k the two-body commutator gives `a+_t` times `k+1` annihilators, and `k+1`
annihilators on `|Psi>` is `a_x` times a k-annihilated state. Same factorisation at every rank.

## Three consequences

**Tangency is derived.** If `X v = 0` then `v^dag Xdot v = i(v^dag S^dag X v - v^dag X S v) = 0`
for ANY `S`. No appeal to reversibility, no kernel projection, no rank argument. This supersedes
the earlier derivation of tangency.

**Exact Cholesky integrability.** With `X = R^dag R`,

    Rdot = -i R S    =>    Xdot = Rdot^dag R + R^dag Rdot = i(S^dag X - X S)

so propagating the factor `R` reproduces the flow exactly while `X = R^dag R >= 0` by construction.
No projection, no semidefinite program, no eigenvalue constraint. Trace conservation is the single
scalar condition `Im Tr(S^dag X) = 0`, which the closure's own tracelessness supplies.

**Free in the interior, binding only at the boundary.** Writing `S = K - i Gamma / 2` gives
`Xdot = i[K,X] - ½{Gamma, X}`. The map `Gamma -> {Gamma, X}` is invertible whenever `X > 0`, so at
any full-rank `X` EVERY Hermitian `Xdot` is representable. The form therefore costs nothing in
accuracy in the interior and becomes automatically tangent as the boundary is approached. This is
the opposite of GKSL, which constrains the flow at every point, and is why GKSL was over-strong.

## The prescription, and what gets tossed

Given any closure's `Xdot`, set `K = 0` and solve `½{Gamma, X} = -Xdot` for `Gamma`. In the
eigenbasis of `X` with eigenvalues `x_i`:

    Gamma_ij = -2 Xdot_ij / (x_i + x_j + sigma)          then      Rdot = -½ R Gamma

`Gamma` is Hermitian because `Xdot` is and the denominator is real symmetric. Then:

- where `x_i + x_j` is healthy, `Gamma` is bounded and the closure's flow is reproduced UNCHANGED;
- where `x_i, x_j -> 0`, `Gamma` would diverge UNLESS `Xdot_ij -> 0`, which is exactly the tangency
  condition. A tangent closure is reproduced exactly; a non-tangent one has precisely its
  non-tangent components suppressed by `sigma`.

So the terms with no tangent representation are discarded, only in the near-kernel directions, and
the discarded amount measures how badly the closure violates positivity. This is NOT a moment
negation: no `gamma_k` is set to zero. It removes the component of the FLOW that has no tangent
representation.

Cost: one eigendecomposition of `X` per step, `O(M^{3k})`, against `O(M^{6k})` plus `n_k^2` flow
evaluations for the Kossakowski test. The scaling obstruction dissolves in this formulation.

## Scope: what XR does and does not cover

| case | covered? |
|---|---|
| closed dynamics, pure state | YES. `X` is rank-deficient, tangency is active and structural |
| closed dynamics, mixed state | YES, and more comfortably: `X` is full rank so the Sylvester solve is well conditioned and `sigma` can be tiny |
| open dynamics (a real bath) | NOT ALONE. `Rdot = -½ R Gamma` preserves `rank(R)`, correct for closed dynamics, but a bath must be able to INCREASE rank. That needs a CP gain term added on top, with rectangular `R` that grows or a stochastic unravelling |

Note that MIXED is not the same as OPEN. A mixed initial state under closed dynamics is fully
covered; only genuine bath coupling needs the CP addition, and there inwardness is physical so
GKSL is the right form.

## What it does not fix

- **Accuracy.** This is a reformulation of the same closure, so the measured cumulant errors
  (0.08 to 1.13 relative) are untouched.
- **Q and G.** XR guarantees `D >= 0`, the pair matrix. The two-hole matrix `Q` and the
  particle-hole matrix `G` are different affine functions of `gamma_2` and `gamma_1` and are NOT
  controlled by this construction. They can still go negative.
- **Extensivity at finite sigma.** Extensivity is inherited from the closure's `Xdot` and is exact
  at `sigma = 0`. The regularisation is a global operation in the eigenbasis of `X`, so it
  introduces an `O(sigma)` extensivity error. For a tangent closure `sigma = 0` is admissible and
  nothing is lost.
- **Convergence of the hierarchy.** Untouched.
- **Nothing here is tested numerically.** The decisive test is the pure-state H4 case where all
  four closures reached `min eig D` of `-0.0059` to `-0.0178` against an exact clearance of
  `4.2e-7`; XR should hold the clearance at the exact value.
