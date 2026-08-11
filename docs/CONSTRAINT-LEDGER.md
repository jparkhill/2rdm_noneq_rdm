# The constraint ledger for a cumulant reconstruction

The objective is not a good closure but a SET OF CONSTRAINTS that fixes one which is stable and
convergent toward exact non-equilibrium RDM dynamics. This file is that ledger. It supersedes the
older per-condition tables in `CLOSURE-CONDITIONS.md`, which remain valid as measurements.

Symbols: `gamma_1[p,q] = <a+_p a_q>`, `eta_1 = I - gamma_1`, `E = gamma_1 - I/2`,
`a_p = n_p (1 - n_p)` with `n_p` the natural occupations; `lambda_k` the MKNO k-cumulant;
`Tr_k T[p_1..p_{k-1}, q_1..q_{k-1}] = sum_y T[p_1..p_{k-1} y, q_1..q_{k-1} y]`;
`anti(E, X)` is `E` contracted into each of `X`'s slots and summed; `A3` is the antisymmetriser
over three upper and three lower indices.

| constraint | type | what it fixes |
|---|---|---|
| antisymmetry, pair-hermiticity | structural | free |
| **connectedness (extensivity)** | physical | the FORM: connected contraction tree; forces `n_lam >= 2` at rank 3 |
| **dimensional analysis** | physical | no bare `V`; inverses only as energy denominators |
| primitive set `{gamma_1, lambda_2, V}` | modelling choice | reopens if a new primitive is admitted |
| polynomiality in `lambda_2` | modelling choice | permits order-by-order matching |
| **trace identity** | physical | the dressing, forcing `a_r^{-1}` |
| **Lipschitz** | physical | NOT violated; `D`/`Q` Cauchy-Schwarz cancels the pole to `O(1)`. Continuity open |
| perturbative limit | **desideratum** | recovery of known limits. Perturbation theories are approximations, not physical conditions |
| N-representability of the manifold | **global qualifier** | every screen is evaluated here, never on free tensors |
| positive generator | physical | not yet analysed; its job is to fix the regularisation |
| convergence `\|lam_{k+1}\|/\|lam_k\| -> 0` | physical | never measured; a bath is the only mechanism on the table |

## Connectedness, derived rather than screened

Cumulants are the logarithm of the moment generating functional. With `: :` vacuum normal
ordering and a source matrix `J`,

    Z[J] = < : exp( sum_pq J_qp a+_p a_q ) : >  =  sum_k (1/k!) J...J gamma_k
    W[J] = log Z[J]                             =  sum_k (1/k!) J...J lambda_k

and the set-partition sum `gamma_k = sum_partitions prod lambda_{|block|}` is the exponential
formula written out. For `rho = rho_A (x) rho_B` on an orbital space `A + B`, every correlation
function with an odd number of A-operators vanishes and `a+_p a_q` is even, so
`Z_AB = Z_A Z_B`, hence `W_AB = W_A + W_B`. Differentiating gives both extensivity clauses at
once: no cross-fragment elements, and no dependence on the spectator.

This is the coupled-cluster argument with the same exponential. CC is extensive because the state
composes multiplicatively, `Omega_AB = Omega_A Omega_B`, so the parameters must compose
additively, and `exp` is the homomorphism between those structures. Here `gamma` is the group
object, `lambda` the algebra object, and the wedge is the product on the group side. A wedge
inside a reconstruction of `lambda_k` is a product inside a logarithm.

**Theorem.** Evaluating a polynomial term on block-diagonal inputs, each factor is confined to one
fragment and a shared contracted index forces two factors into the same fragment. The term's value
is a sum over assignments of the connected components of its contraction graph to fragments, and
only a single-component graph admits no mixed assignment. So

    additive under direct sum   <=>   the contraction diagram is connected.

## What the constraints leave at rank 3

Leg counting with connectivity: with `n_lam` factors of `lambda_2`, `n_g` of `gamma_1` and `c`
contractions, `2 n_lam + n_g - c = 3` and `c >= n_lam + n_g - 1`. For `n_lam = 1` this needs
`c = n_g - 1 < n_g`, impossible. So **`n_lam >= 2`**, with no perturbative input. At `n_lam = 2`
the two inequalities are equalities, so the diagram is a TREE: exactly one internal line joining
the two `lambda_2` vertices, plus dressings. In the natural-orbital basis every dressing is a
scalar function of one occupation, so

    lambda_3[pvx,qwy] = A3[ sum_r  d_out(n_p)... lambda_2[pv,qr] d_0(n_r) lambda_2[rx,wy] ... ]

The unconstrained `lambda_3` has `C(M,3)^2` components (3136 at M=8); connectedness plus
dimension replaces that by a handful of scalar functions of `n`.

## The dressing, and why the pole is not fatal

Any linear-in-`lambda_2` term in `Tr_3` of a connected structure arises only by closing one
`lambda_2` into a loop, and the rank-2 trace identity fixes that loop exactly,
`sum_y lambda_2[ry,wy] = -(gamma_1 eta_1)[r,w]`, weight `a_r`. The target `anti(E, lambda_2)`
carries weight `E`, which is `+-1/2` where `a_r = 0`. So trace-exactness forces
`d_0(n) ~ (n - 1/2)/a_n`, singular at both occupation boundaries.

That does NOT break Lipschitz, because the screen must be applied on the N-representable
manifold. `D[(p,v),(q,w)] = 2 gamma_2[p,v,q,w]` is a Gram matrix of `zeta_ab = a_a a_b |Psi>`, so
`|gamma_2[pv,qw]| <= sqrt(gamma_2[pv,pv] gamma_2[qw,qw])` with `gamma_2[pv,pv] = <n_p n_v> <=
min(n_p,n_v)`; the two-hole matrix `Q` gives the mirror bound as `n_r -> 1`. Uniformly, every
`lambda_2` element carrying index `r` is `O(sqrt(a_r))`, so

    lambda_2 . a_r^{-1} . lambda_2  =  O(sqrt(a_r)) a_r^{-1} O(sqrt(a_r))  =  O(1)

The pole cancels exactly, the two square roots supplied by the `D` and `Q` conditions. An earlier
claim in this programme that connectedness, trace-consistency and Lipschitz were mutually
incompatible was wrong, and wrong specifically because it screened the dressing in isolation.

Three things remain open, and they are the honest residue:

1. **Bounded is not continuous.** The limit has the form `xy/(x^2+y^2)`: finite but possibly
   direction-dependent. What decides it is whether `lambda_2`'s `r`-components scale
   proportionally to `sqrt(a_r)` or are merely bounded by it.
2. **Magnitude.** The cancellation leaves `O(1)` in the weak-correlation limit. This is a
   limit-matching concern, not a stability one, and sits on the desideratum row.
3. **The residual bilinear part** of `Tr_3(chain)` must vanish or be absorbed; this has not been
   verified symbolically, and is a finite MKNO calculation.

## Retrodictions, none of them fitted

Every extensivity measurement follows from the connectedness rule: `null` and the connected
chains leak exactly `0.000e+00`; Mueller Eq. 14, Lackner Eq. 43 and `muller_rdmft` leak
`0.73` to `1.15`. And in the rank-4 trace-rule fit (`tools/trace_rule_pattern.py`, run earlier for
an unrelated purpose over eight candidate structures with coefficients common to four systems) the
four CONNECTED structures took nonzero coefficients and the four WEDGE structures came out
exactly zero, at relative residual 8.6e-15.
