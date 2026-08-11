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
| **complete positivity (GKSL)** | physical, REQUIRED | `Kossakowski >= 0`. The only structural guarantee of positivity for an autonomous generator |
| size of the dissipator | the dial, not a constraint | 0 gives a Hamiltonian reduced flow, tangent and reversible; matched to a bath gives the open case |
| Nagumo tangency | physical | implied by CP: `P_K L[X] P_K = sum_i (P_K L_i) X (P_K L_i)^dag >= 0` |
| clearance erosion rate | diagnostic | the REALISED erosion on a given trajectory; bounded above by `-lambda_min(Kossakowski)` |
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

## Complete positivity: required, and satisfied by the exact solution

The positive-generator test used is conditional complete positivity: a generator `L` is GKSL iff
it is Hermiticity-preserving, trace-annihilating, and `Q Choi(L) Q >= 0` with
`Choi(L) = sum_ij |i><j| (x) L[|i><j|]` and `Q = 1 - |Omega><Omega|`,
`|Omega> = n^{-1/2} sum_i |i>|i>`. The eigenvalues of `Q Choi(L) Q` are the Kossakowski spectrum.

THE EXACT SOLUTION SATISFIES IT. At `k = N` the reduced dynamics genuinely is an autonomous
generator -- the von Neumann commutator -- and `tools/exactness_limit.py` finds its Kossakowski
spectrum at `+-1.1e-15` against a Choi norm of 1.347, i.e. positive semidefinite to round-off.
Machinery checks alongside: Hermiticity preserving 1.6e-15, trace annihilating 6.3e-16, and the
rank-2 Lyapunov drift alone at `+-1.5e-15`, confirming in CP language that all positivity risk
lives in the `4k` terms.

CP CONTAINS REVERSIBLE DYNAMICS; IT DOES NOT EXCLUDE IT. At a boundary point with kernel projector
`P_K`, using `X P_K = P_K X = 0`, both the commutator and the anticommutator parts compress to
zero and

    P_K L[X] P_K  =  sum_i (P_K L_i) X (P_K L_i)^dag  >=  0

with equality iff every `L_i` maps `range(X)` into itself -- exactly what rank-preserving unitary
dynamics does. Kossakowski `= 0` IS a purely Hamiltonian generator. So CP implies Nagumo, contains
tangency as its boundary case, and forbids only outward flow.

Kossakowski spectra of the linearised rank-2 flow under each closure, H3 bent, same state:

    generator          min eig      max eig    |Choi|_max
    drift only      -1.494e-15    1.494e-15     1.386
    null            -7.698e-01    7.420e-01     1.354
    trace           -7.029e-01    7.016e-01     1.346
    scaled          -1.784e-01    1.742e-01     1.378
    scaled+trace    -5.857e-01    5.851e-01     1.346

Every existing closure is badly non-CP. Read this as a figure of merit, not as an argument against
the criterion: the task is to build a closure with `Kossakowski >= 0`, which is a semidefinite
feasibility problem over the 4-dimensional connected span.

THE TWO gamma* ARE NOT IN CONFLICT. `gamma*_CP = -lambda_min(Kossakowski) = 0.18 to 0.77` hartree
is the WORST-CASE clearance erosion rate, attained when the state sits on the PSD boundary along
the offending direction. `gamma*_rate ~ 1.6e-4 to 7.6e-4` hartree is the REALISED rate on one
benign interior thermal trajectory. Necessarily `gamma*_rate <= gamma*_CP`; the gap measures how
far that trajectory stayed from the bad direction, not that the criterion is too strong. The
confirmation is in the propagation data: on the PURE state, near the boundary, all four closures
went negative, ordered as their CP defects predict.

(An earlier version of this file recorded the opposite conclusion, that CP should be rejected as
too strong. That was wrong on all three counts and is retracted.)
