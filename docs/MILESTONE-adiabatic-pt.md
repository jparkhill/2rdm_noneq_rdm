# Reference-free perturbative cumulant from the MKNO EOM

## What was built

`tools/hylleraas.py` derives `lambda_2` by adiabatic elimination of the MKNO rank-2 EOM, with no
reference determinant anywhere. The term list `gwt.terms_cumulant(2)` is partitioned twice:

1. **By degree in `lambda_2`.** Degree 0 is the driving `D[gamma_1]` (8 one-body + 56 two-body),
   degree 1 is the operator (4 + 76). 8 terms carrying `lambda_3` are truncated; 0 terms are
   nonlinear in `lambda_2`, so the equation is exactly linear once `lambda_3` is dropped.
2. **By how many of the interaction operator's indices the `lambda_2` factor touches.** One shared
   index is the Fock/mean-field channel and belongs in the zeroth-order operator `L0`
   (4 + 16 terms). Two shared indices is the ladder/ring channel, first order in `V` beyond `L0`
   (0 + 60 terms).

Setting `lambda_2-dot = 0` gives `L0[lambda_2] = -D`, which is the stationarity condition of the
Hylleraas functional `J = 1/2 <l2, L0 l2> + <l2, D>`. No determinant enters: the occupied/virtual
projectors of ordinary MP2 are replaced by `gamma_1` and `eta_1 = 1 - gamma_1` inside the MKNO
term list, and the energy denominators by the spectrum of `L0` built from the generalised Fock
operator `f = h + V.gamma_1`.

## Acceptance gate: PASSED

At the idempotent (HF) density the solved cumulant equals the MP1 amplitude exactly.

    vv-oo block elementwise ratio solved/MP1 : mean 1.000000   std 3.77e-08
    block norms: vv,oo 0.237010  oo,vv 0.237010  ovov 0  vvvv 0  oooo 0

Every second-order block is exactly zero, as a first-order cumulant requires. The energy is
`(1/4) V.lambda_2 = -0.066918314` against PySCF `E_MP2 = -0.033459145`, a ratio of `+2.000001`;
the factor 2 is the standard difference between the variational (Hylleraas) second-order energy
and the expectation-value formula evaluated with a first-order cumulant, since
`lambda_2 = t + t^dag` contributes to `(1/4) V.lambda_2` twice. `J_min` equals `E_MP2` exactly.

Keeping the ladder terms in `L0` instead gives ratio `-1.322572`: that solve is a self-consistent
RPA-like resummation, not MP2. The second partition is what makes the scheme perturbative.

## Cross-validation

At the correlated FCI density the reference-free construction reproduces `closures.lambda2(...,
"mp2")`, which was built by an entirely different route (explicit occupation factors in the
natural-orbital basis): cumulant error `0.08790` vs `0.08798`, trace-identity residual `0.03061`
vs `0.03083`. Both routes validate each other.

## The structural result: the two condition families are complementary

Measured on the jittered asymmetric-trapezoid H4/sto-3g, `M=8 N=4`, at the FCI density:

    closure      |l2 - l2_exact|   trace-id residual   E - E_FCI
    hylleraas        0.08790            0.03061         +0.077460
    mp2              0.08798            0.03083         +0.079831
    blocked          0.13998            0.00000         +0.100005
    zero             0.13998            0.03083         +0.141855

`hylleraas`, `mp2` and `lambda_2 = 0` all sit at the same trace-identity residual `~0.0308`: the
amplitude sector contributes nothing to the contraction. That is forced, not accidental --
contracting `lambda_2[ab,ij]` needs an index that is simultaneously occupied and virtual, so the
contraction of the amplitude block vanishes identically (measured `0.00e+00`). Conversely
`blocked` fixes the trace exactly while leaving `|l2 - l2_exact|` equal to `|l2_exact|`, i.e. it
captures none of the dominant block, which carries 68.3% of the exact cumulant's Frobenius norm.

`ker L0` is 100-dimensional (`rank 684` of `784`) and the minimum-norm solve puts nothing in it.

**So stationarity fixes `range(L0)`, where the norm is; the contraction identities fix part of
`ker L0`, where the trace is. They constrain orthogonal subspaces and cannot compete.** This
retires the question of which condition is "right": the single-trace family we derived at ranks
2, 3, 4 is order-general but structurally blind to the leading block of every cumulant, and was
never capable of fixing it alone.

## Naming

- `pt` -- the amplitude alone, `L0[lambda_2] = -D`. Verified against MP2 above.
- `trace_plus_pt` -- `pt` plus a `ker L0` component chosen to satisfy the partial-trace identity.
  Not yet built; this is the next step, and it should inherit `pt`'s amplitude accuracy and
  `blocked`'s exact trace simultaneously.

## Why this replaces `blocked`

`blocked` was shown to have three independent defects, all traceable to enforcing the trace
identity by a *global* rescaling `f_i = -a_i/[(1+Q)(T-2a_i)]`:

- **A pole surface inside the feasible set.** `T - 2a_i = 0` whenever one orbital carries half the
  total blocking, e.g. any density with two equally fractional occupations. `min E` over
  `0 <= gamma_1 <= 1` is therefore ill-posed, and the constrained solve sat at a saddle with
  471 of 1200 random perturbations lowering the energy.
- **Non-extensivity.** `T` is a global sum, so adding a non-interacting fragment changes the
  cumulant here: cross-fragment `|l2| = 0.021` and within-fragment drift `0.029` / `0.056`.
  Mueller fails only in the cross terms (`0.112`) with no within-fragment drift.
- **No `V`-dependence**, hence no perturbative limit and no coherences.

The adiabatic construction has none of these: `L0` and `D` are built from connected contractions
of `h`, `V`, `gamma_1`, so both are fragment-block-diagonal and so is the (pseudo-)inverse, which
makes the solution block-additive; and the amplitude is `V`-linear by construction.

## Higher ranks

The same partition applies to `terms_cumulant(k)`: `lambda_k-dot = L0[lambda_k] + D_k[gamma_1,
lambda_2..lambda_{k-1}] + ...`, truncate `lambda_{k+1}`, solve. Each step remains a *linear*
solve because lower cumulants are already known. The MBPT order counting falls out unforced:
`D_2 ~ V` gives `lambda_2 = O(V)`, `D_3` contains `V.lambda_2` giving `lambda_3 = O(V^2)`, and
generally `lambda_k = O(V^{k-1})`, which is the standard order of a connected `k`-body amplitude.
Truncating `lambda_{k+1}` costs `O(V^k)`, one order beyond the leading term, so the truncation is
order-consistent at every rank.

Expected difficulties at higher rank: `ker L0` grows, so more of `lambda_k` is left unfixed by
stationarity; dense construction becomes infeasible (`[C(M,k)^2]^2`) and needs Krylov solves with
the cheap einsum application; and `L0` is Hermitian only in the `gamma_1 eta_1`-weighted metric,
which degrades near integer occupations. One thing must be *checked* at each rank rather than
assumed: `D_k` must lie in `range(L0)`, or no adiabatic solution exists and that cumulant has to
be propagated instead. At rank 2 it does (`|L l2 + D| = 4.4e-16`).

## Trap worth recording

`lambda_2-dot = i * (term sum)`. The `i` is load-bearing: the term sum is anti-Hermitian in the
pair index, so omitting it makes `L0` map entirely into the anti-Hermitian sector, exactly
orthogonal to a Hermitian pair basis. Every projection then returns `~1e-16` and `L0` looks
identically zero. The same trap exists at every rank and anywhere the EOM is projected onto a
Hermitian basis.
