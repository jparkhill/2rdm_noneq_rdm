# The G-transfer form of the hierarchy

The resolution of the gauge obstruction in `XR-POSITIVITY.md`. The exponentially tall isometry `P`
enters `gamma_{k+1}` only through small blocks, so the whole construction stays polynomial.

## Derivation

With `u_I = a_{i_k}..a_{i_1}|Psi>` the k-annihilated states (`I` a sorted k-tuple), these live in the
`(N-k)`-particle sector and `Mat_k = k! Gram(u)`. Write `R = sqrt(k!) * (coefficients of u)` and
factor `R = P R_c` with `P` an isometry onto `range(R)`, `r = rank R`, so `R_c` is `r x C(M,k)` and
`X = R^dag R = R_c^dag R_c`.

The `(k+1)`-annihilated states are ONE more annihilator applied to those,
`w_L = a_{l_{k+1}} u_{(l_1..l_k)}` with `L` sorted so `l_{k+1}` is the largest element. Hence

    Mat_{k+1}[(I,x),(J,y)]  =  (k+1) * R_c(:,I)^dag  G_xy  R_c(:,J)
    G_xy  =  P^dag a_x^dag a_y P                                       an r x r block

Assembled as an `(M r) x (M r)` matrix, `G` is the Gram matrix of `{a_x e_alpha}` for `{e_alpha}` a
basis of `range(R)`. Three properties follow, and each replaces a constraint we previously had to
impose by hand:

| property of `G` | what it delivers |
|---|---|
| `G >= 0` (it is a Gram matrix) | positivity of `gamma_{k+1}`, structural |
| `sum_x G_xx = (N-k) I_r`, since `sum_x a_x^dag a_x = Nhat` and the `e_alpha` lie in the `(N-k)`-sector | `Tr_{k+1} gamma_{k+1} = (N-k) gamma_k` EXACTLY -- the moment form of the single-trace cumulant identity |
| at `N = k` the prefactor is 0, so `G >= 0` with zero trace forces `G = 0` | `gamma_{k+1} = 0`: exactness at the top, automatic, with no global prefactor and hence nothing non-extensive |

The single-trace identity cost a great deal of effort to impose directly, forced the `a_r^{-1}` pole,
and broke extensivity in every wedge implementation. Here it is one linear condition on `G`.

## Sizes

`r = rank(Mat_k) <= min(C(M,k), C(M,N-k))`, so `G` is at most `M*C(M,k)` square: polynomial in `M` at
fixed `k`. The exponentially large `(N-k)`-particle sector never appears in the closure data.

## Verification

`tools/G_transfer.py`, exact FCI states at a NON-STATIONARY point (an eigenstate would make the
dynamical comparisons vacuous):

    case                    |X - Rc^dag Rc|  |Mat_{k+1} - G form|  min eig G   |sum_x G_xx - (N-k)I|
    H3 bent k=2  M=6 N=3          2.220e-15             2.286e-15  -1.620e-15             6.661e-16
    H4 trap  k=2  M=8 N=4         9.998e-16             1.578e-15  -8.852e-15             2.442e-15
    H3 bent k=1  M=6 N=3          6.661e-16             8.882e-16  -1.424e-15             1.776e-15

Sizes found: `r = 6` of 15, `28` of 28, `6` of 6; `G` is `36 x 36`, `224 x 224`, `36 x 36`. Note the
rank bound bites unevenly -- at H3 with `k=2`, `rank Mat_2 = 6 = C(6,1)` against `C(6,2) = 15`,
while at H4 with `k=2` it is full.

## What the closure problem has become

Model `G`: an `(M r) x (M r)` positive semidefinite block matrix with `sum_x G_xx = (N-k) I_r`. That
is the entire remaining freedom, and it is polynomial.

The simplest ansatz is `G = D (x) I_r` with `D` an `M x M` matrix, `D >= 0`, `tr D = N-k`. It gives

    Mat_{k+1}[(I,x),(J,y)] = (k+1) D_xy X[I,J]

i.e. `gamma_{k+1}` as an antisymmetrised product of `D` with `gamma_k`. That is a WEDGE, hence
disconnected, hence non-extensive -- so this ansatz reproduces the Lackner/Mueller-type
reconstructions and inherits their leak. The G-transfer form therefore CONTAINS the known
reconstructions as its crudest member, and extensivity remains the discriminator among `G` ansatze,
now acting on an `M r` square matrix instead of a rank-`(k+1)` tensor.

Interpretation of `D`: taking `D_xy = Tr_r(G_xy)/r` makes `D` the one-particle density matrix of the
`(N-k)`-particle states spanning `range(R)` -- a HOLE-side 1-RDM, with `tr D = N-k` as it must be.

## Open

- No `G` ansatz beyond `D (x) I_r` has been constructed, and that one leaks.
- Extensivity of a general `G` ansatz has not been characterised; the connectedness rule was derived
  for cumulant contractions, and its translation into a condition on `G` is not yet written down.
- Accuracy is entirely untested.
- The dynamics is not implemented. `G` must either be propagated or reconstructed at each step, and
  which of those is correct has not been settled.
- `Q` and `G`-cone N-representability (the two-hole and particle-hole matrices) remain uncontrolled.
  Note the unfortunate name collision: the `G` of this file is the transfer block matrix, NOT the
  particle-hole N-representability matrix.

## What extensivity demands of G: the ladders must be fragment-local

DERIVATION. For two non-interacting fragments in a product state the k-annihilated states factorise.
Splitting a k-tuple `I` into `I_A`, `I_B` with `j = |I_A|`,

    u_I = +- (A_{I_A}|Psi_A>) (x) (A_{I_B}|Psi_B>)

so `range(R) = (+)_j V_A^{(j)} (x) V_B^{(k-j)}` with `V_f^{(j)}` spanned by the j-fold annihilations
of fragment f. Since annihilators map `V_f^{(j)}` INTO `V_f^{(j+1)}` exactly, the compressed ladders
`frak_a_x = a_x` restricted to `range(R)` come out factorised:

    x in A :  frak_a_x = alpha^A_x (x) I_B            raises j by 1
    y in B :  frak_a_y = Pi_A (x) alpha^B_y           leaves j alone, Pi_A = (-1)^{N_A} the parity

Same-fragment blocks of `G` are then (ladder^dag ladder on that fragment) tensor identity on the
other; cross-fragment blocks are (ladder^dag on one) tensor (ladder on the other). Hence

    EXTENSIVITY  <=>  G_xy = frak_a_x^dag frak_a_y with every frak_a_x FRAGMENT-LOCAL.

This is the same statement as "connected contraction" for cumulants and "fragment-local jump
operator" for a GKSL generator, in a third language. The ladder form also gives `G >= 0` for free,
since `G` is then a Gram matrix of the `frak_a_x`, and `sum_x frak_a_x^dag frak_a_x` is the modelled
number operator, so the contraction identity is one condition on the ladders.

VERIFIED, tools/G_locality.py, H3 bent x H2 as an exact product state, M = 10, N = 5, k = 2.
Every row satisfies `G >= 0` and `sum_x G_xx = (N-k) I_r` to 1e-14, so locality is the only
difference:

    G ansatz                                    |inside|        leak
    exact G (local ladders)                    1.764e-02   1.826e-13
    V^dag G V, V random on range(R)            4.497e-01   1.610e+00
    V^dag G V, V random but j-block diagonal   7.594e-01   9.847e-01
    G = D (x) I_r, D the hole 1-RDM            4.633e-01   1.278e+00

Preserving the particle-number grading is NOT sufficient (leak 0.985); the tensor factorisation
inside each graded block is what is required. The dimension count confirms the structure
independently: `r = 31 = 6 + 24 + 1 = sum_j dim V_A^{(j)} dim V_B^{(k-j)}` for `j = 2,1,0`, with the
A-particle numbers on `range(R)` measured as `[1 x6, 2 x24, 3 x1]`.

## The construction, complete

STATE: `R_c`, an `r x C(M,k)` matrix, with `X = Mat_k = R_c^dag R_c`.
CLOSURE DATA: ladder operators `{frak_a_x}`, one per spin orbital, from the level-k compressed space
to level `k+1`, subject to

    sum_x frak_a_x^dag frak_a_x = (N-k) I_r        and        each frak_a_x fragment-local

RECONSTRUCTION: `G_xy = frak_a_x^dag frak_a_y`, then
`Mat_{k+1}[(I,x),(J,y)] = (k+1) R_c(:,I)^dag G_xy R_c(:,J)`.
PROPAGATION: XR, `Rdot = -i K R - i R S`, with `K` the compressed Hamiltonian supplying the gauge.

Every derived row of the ledger is then structural:

| ledger row | how it is satisfied |
|---|---|
| antisymmetry, pair-hermiticity | the sorted-tuple representation |
| connectedness (extensivity) | ladders fragment-local; leak 1.8e-13 against 0.98-1.61 otherwise |
| dimensional analysis | ladders are dimensionless; no bare `V` |
| trace identity | `sum_x frak_a_x^dag frak_a_x = (N-k) I_r` |
| Lipschitz | polynomial in `R_c` and the ladders; no inverse of anything |
| exactness at `k >= N` | at `N = k` the constraint forces every `frak_a_x = 0`, so `gamma_{k+1} = 0` |
| positivity of `D` | `X = R_c^dag R_c` |
| positivity of `gamma_{k+1}` | `G = frak_a^dag frak_a` |
| Nagumo tangency | derived from `Xdot = i(S^dag X - X S)` |

WHAT REMAINS is entirely quantitative, and it is one object: the ladders themselves. `r` is
polynomial even in the exact case (`r = rank Mat_k <= min(C(M,k), C(M,N-k))`), so nothing exponential
is required; the approximation is the choice of `{frak_a_x}`, whose defect is measurable as
`||{frak_a_x, frak_a_y^dag} - delta_xy||`. Accuracy is untested, the ladder equation of motion is not
derived, and `Q` and the particle-hole N-representability matrix remain uncontrolled.

## CORRECTION: the number condition does NOT give the contraction identity

An earlier section of this file, and commits e256383 / 08462ee / f84d58f, claimed that
`sum_x frak_a_x^dag frak_a_x = (N-k) I_r` delivers `Tr_{k+1} gamma_{k+1} = (N-k) gamma_k` exactly.
THAT IS WRONG for a MODELLED ladder. It holds for the exact ladder only because the exact ladder
also satisfies CAR anticommutation.

Measured, H4 trapezoid M=8 N=4, level 1, ladders drawn at random from the Stiefel manifold so that
`sum_x frak_a_x^dag frak_a_x = (N-1) I` holds by construction:

    ladder                       tr Mat_2 |Tr2 g2-(N-1)g1|    min eig D     (N(N-1) = 12)
    random chi'=4               15.787393      1.563e+00   -1.694e-15
    random chi'=8               15.977133      1.684e+00   -6.826e-16
    random chi'=28              16.148988      1.755e+00    1.605e-05
    EXACT ladder                12.000000      4.441e-15    4.214e-07

The cause is that the G-transfer formula peels the LARGEST index of the (k+1)-tuple. The contraction
identity requires peeling arbitrary indices and getting the same answer, which is exactly
`{frak_a_x, frak_a_y} = 0` across levels. Without it the construction is peel-order dependent and the
trace comes out wrong.

CORRECTED ACCOUNTING of what the ladder ansatz gives for free:

| property | automatic? |
|---|---|
| `gamma_{k+1} >= 0` | YES. `G = frak_a^dag frak_a` is a Gram matrix; `min eig D >= 0` in every sample |
| exactness at `k = N` | YES. `sum_x frak_a_x^dag frak_a_x = 0` with PSD terms forces every `frak_a_x = 0` |
| contraction identity | NO. Must be imposed |
| hole condition `Gtilde <= I` | NO. It followed from CAR, so it goes too |
| extensivity | NO. Needs the locality condition |

The contraction identity is however expressible AT A SINGLE LEVEL, so a one-rung closure can still
impose it without touching level k+1. Writing `sum_x gamma_2[p,x,q,x]` in the antisymmetric extension
gives `sum_{x > max(p,q)} Rc(:,p)^dag G_xx Rc(:,q)` plus `sum_{x < min(p,q)} Rc(:,x)^dag G_pq Rc(:,x)`
plus mixed-sign terms -- a LINEAR condition on `G`. Likewise the hole condition. So the closure
remains a well-posed problem: linear constraints plus the conic/rank structure `G = frak_a^dag frak_a`,
i.e. a constrained Stiefel optimisation or a rank-constrained SDP.

CONSEQUENCE for the level-1 experiment (tools/ladder_optimise.py): the ORACLE objective, fitting the
exact Mat_2, remains valid. The VARIATIONAL objective without the contraction constraint is
meaningless and was observed running away to E = -4.87 against E_FCI = -1.943, because the
unconstrained family contains gamma_2 with tr gamma_2 ~ 16 instead of 12.

## Exact CAR forces an active space, so a CAR-exact ladder is not a closure

Two results that together end the ladder route as a vehicle for the trace identity.

FIRST, the peel convention is not the problem. The G-transfer formula peels the largest index, which
is arbitrary; the peel-symmetric form is

    Mat_{k+1}[L,L'] = c sum_{i,j} (-1)^{i+j} R_c(:,L\l_i)^dag G_{l_i l'_j} R_c(:,L'\l'_j)

and at k = 1 with c = 1/2 it reproduces the exact Mat_2 to 3.8e-15, so the symmetrisation is
consistent. But it does NOT rescue the contraction identity for a modelled ladder -- it makes it
worse (H4, N(N-1) = 12):

    ladder / formula                  tr Mat_2   |Tr2 g2-(N-1)g1|
    random chi'=8   single-peel      16.586770         1.523e+00
    random chi'=8   PEEL-SYMMETRIC    4.890902         2.179e+00
    random chi'=28  single-peel      16.844916         1.558e+00
    random chi'=28  PEEL-SYMMETRIC    5.217381         2.083e+00

So the failure recorded above is a genuine defect, not a convention error.

SECOND, and this is the reason. Collect the levels into W = (+)_k V^{(k)}. Then {frak_a_x, frak_a_y}
= 0 says the ladders generate an EXTERIOR ALGEBRA action on W, whose general form is
frak_a_x = sum_i C_xi b_i with b_i standard anticommuting generators. Imposing the number condition,

    sum_x frak_a_x^dag frak_a_x = sum_ij (C^dag C)_ij b_i^dag b_j = (N-k) I   requires  C^dag C = I_m

and C^dag C = I_m forces m <= M. So the ladder is an isometric embedding into an m-mode Fock space
with m <= M: an ACTIVE SPACE.

    Exact CAR plus the number condition makes the ladder a wavefunction truncation (CASCI),
    not a BBGKY closure. Any genuine closure must break CAR, and once CAR is broken the trace
    identity is no longer structural and must be imposed explicitly.

That is the structural cause of the retreat into penalties and semidefinite programs recorded in
ccc7791, 840606f and 90b66c8: the trace identity was never going to be free in this representation.
Stated as an argument, not a theorem -- reducible and non-isometric representations have not been
ruled out.

CONSEQUENCE FOR THE PROGRAMME. The trace identity's natural home is the cumulant form
Tr_k lambda_k = anti(E, lambda_{k-1}) + chains, where it is a LINEAR condition on the coefficients of
connected structures rather than something to penalise. The connected-chain construction lives there:
extensive by index structure, correct perturbative order, and a 4-dimensional span containing 99.4%
of lambda_3. Its one unfinished piece is whether the bilinear residue of Tr_3(chain) can be absorbed,
which is a finite MKNO calculation and was never done.
