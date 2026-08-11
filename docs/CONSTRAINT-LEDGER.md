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
| **exactness at k >= N** | physical | `gamma_{k+1} == 0` when `k >= N`, so `lambda_{k+1} = -separable`. `null` fails this maximally |
| perturbative limit | **desideratum** | recovery of known limits. Perturbation theories are approximations, not physical conditions |
| N-representability of the manifold | **global qualifier** | every screen is evaluated here, never on free tensors |
| complete positivity (GKSL) | **rejected as a design condition** | sufficient for positivity but far too strong; see below |
| Nagumo tangency | physical | exact, but only bites on the PSD boundary, which is measure zero for mixed states |
| clearance erosion rate | physical | the operative stability condition for mixed states |
| convergence `\|lam_{k+1}\|/\|lam_k\| -> 0` | physical | never measured; a bath is the only mechanism on the table |

## Exactness at k >= N, and why it is not in tension with extensivity

For an N-electron state `gamma_{N+1} == 0` identically, since `N+1` annihilation operators kill it.
So at `k = N` the `4k` terms vanish, the rank-N equation of motion closes with no approximation,
and it must reproduce the von Neumann equation. Verified in `tools/exactness_limit.py`: H3 bent,
M = 6, N = 3, at a NON-STATIONARY state (the ground state gives `d/dt gamma_3 = 0` and the test
would be vacuous), the rank-3 wick term list with `("G",4)` set to zero reproduces
`fullci.rdm_dot` at 2.2e-16 against `|d/dt gamma_3| = 0.188`.

In cumulant language the condition reads `lambda_{k+1} = -separable(gamma_1, ..., lambda_k)` for
`k >= N`, which is NOT zero. So `null` violates it maximally, and Mueller's `(N-2)` prefactor
exists precisely to satisfy it at rank 3, `N = 2`.

There is no conflict with extensivity, though there appears to be one. A global `N`-dependent
prefactor is non-extensive; but the condition itself is about each fragment's own `N_f`, and
`lambda^{AB} = lambda^A (+) lambda^B` handles that automatically. Two 2-electron fragments have
`gamma_3^f = 0` each, while the composite has `gamma_3^{AB} != 0` entirely through the separable
part. So the condition is extensive; only its implementation as a global moment-level prefactor
is not. That is exactly why Mueller fails screen 2 while the condition it encodes is sound.

## Why complete positivity is the wrong design condition

The positive-generator test used is conditional complete positivity: a generator `L` is GKSL iff
it is Hermiticity-preserving, trace-annihilating, and `Q Choi(L) Q >= 0` with
`Choi(L) = sum_ij |i><j| (x) L[|i><j|]` and `Q = 1 - |Omega><Omega|`,
`|Omega> = n^{-1/2} sum_i |i>|i>`. The eigenvalues of `Q Choi(L) Q` are the Kossakowski spectrum.

Machinery validation, H3 bent, `Mat_3` of dimension `C(6,3) = 20`: the exact rank-N generator is
the von Neumann commutator, and its Kossakowski spectrum comes out `+-1.1e-15` against a Choi norm
of 1.347. The rank-2 Lyapunov drift alone gives `+-1.5e-15`, confirming in CP language that all
positivity risk is in the `4k` terms.

Kossakowski spectra of the linearised rank-2 flow under each closure, same system and state:

    generator          min eig      max eig    |Choi|_max
    drift only      -1.494e-15    1.494e-15     1.386
    null            -7.698e-01    7.420e-01     1.354
    trace           -7.029e-01    7.016e-01     1.346
    scaled          -1.784e-01    1.742e-01     1.378
    scaled+trace    -5.857e-01    5.851e-01     1.346

Every spectrum is essentially symmetric, `|min| ~ |max|`, so the Kossakowski matrix is nearly
traceless: these generators are oscillatory non-Markovian corrections, not net-dissipative ones.
That is what reversibility requires. No closure is remotely CP, and this is NOT a defect, because:

1. At `k < N` the exact reduced flow is not autonomous -- `gamma_{k+1}` is independent data, not a
   function of `gamma_k` -- so there is no exact generator, and asking whether it is CP-divisible
   is a category error. CP cannot be validated against the exact answer at `k < N`.
2. Closed-system reduced dynamics is tangent at the PSD boundary, never inward, so a nonzero CP
   gain part is forbidden for a closed system. Imposing GKSL forces either a zero dissipator --
   back to mean field -- or spurious damping.

The quantitative consequence. Reading `gamma*` as the depolarising bath strength needed to make
the Kossakowski matrix PSD gives `gamma*_CP = -lambda_min = 0.18 to 0.77 hartree`, which is 5 to
20 eV as a relaxation rate: unphysical for a molecule. The rate-competition figure measured by
propagation, clearance erosion against the exact clearance, gives `gamma*_rate ~ 1.6e-4 to 7.6e-4`
hartree on a thermal state -- THREE ORDERS OF MAGNITUDE smaller, and the thermal propagation
indeed stayed positive to `t = 1` with no bath at all.

So CP overestimates the required dissipation by ~10^3 and must not be used as the design
condition. The operative condition is the erosion rate against the exact clearance.

One suggestive but uncontrolled observation: the CP-defect ordering (`scaled` 0.178 best, `null`
0.770 worst) matches the mixed-state erosion ordering at both endpoints (`scaled` 1.6e-4 best,
`null` 7.6e-4 worst) but swaps the middle two. The two measurements are on different systems and
states, so this is not a controlled comparison and should not be relied on.

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
