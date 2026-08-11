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
