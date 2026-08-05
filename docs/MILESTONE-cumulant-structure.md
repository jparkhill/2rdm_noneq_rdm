# Milestone: cumulant trace identities and where blocking lives

Recorded state of understanding. Every statement below names its derivation method and is
backed by a tool in `tools/`. Ground truth throughout is `qdiss.eom.fullci`, exact Fock-space
linear algebra; the two symbolic derivations `qdiss.eom.wick` (plain vacuum Wick, bare
moments) and `qdiss.eom.gwt` (MKNO generalized normal ordering) are each validated against it
at 1e-15 by `qdiss.eom.validate` and never against each other.

## 1. The two derivations agree, untruncated

`tools/gno_vs_bare.py`. The GNO route collapses to the identical term list as plain Wick,
4 one-body and 12 two-body terms with factors D1 x16, G2 x4, G3 x8, `|bare - GNO| = 0.0`.
Pushed to cumulant form the rank-2 equation is 12 + 140 terms over gamma_1, lambda_2 and
lambda_3, of which 8 carry lambda_3, matching at 1.8e-15.

## 2. Closures, named by procedure

A gamma-zero closure means: start from the unapproximated MKNO EOM of the relevant rank,
solve the identity `gamma_{k+1} = 0` for the top cumulant, substitute back. Nothing is set to
zero by hand.

- `MK_gamma30` — `lambda_3 := -(A + B)`, state `(gamma_1, lambda_2)`, `gamma_1-dot` exact.
  Exact for two electrons: RK4 convergence ratios 16.1 and 16.0, so the closure error is
  identically zero (`tools/mk_gamma30.py`, `tools/mk_gamma30_convergence.py`).
- `MK_gamma20` — `lambda_2 := -(gamma_1 wedge gamma_1)`. The 4 lambda_2 terms cancel the 8
  gamma_1-product terms, the entire two-body contribution collapsing to 8.9e-16, leaving
  `i[h^T, gamma_1]`: free evolution. NOT TDHF (`tools/mk_gamma20_procedure.py`).
- TDHF — the same unapproximated MKNO 1-RDM EOM with `lambda_2 := 0`. Neglecting the
  CUMULANT leaves `gamma_2 = gamma_1 wedge gamma_1`, the determinant form, which is the mean
  field; the 8 surviving terms assemble into `W[p,q] = V[p,r,q,s] gamma_1[r,s]` and the
  equation becomes `i[(h+W)^T, gamma_1]`, confirmed at 3.9e-16 (`tools/tdhf_identify.py`).

Neglecting the cumulant is not neglecting the 2-RDM. `lambda_2 = 0` gives TDHF, `gamma_2 = 0`
gives free evolution, and they differ by 0.867.

General pattern: `gamma_{k+1} = 0` yields exactly `i[H_k^T, gamma_k]` with `H_k` the bare
k-particle Hamiltonian, verified at 8.9e-16 for k=1 and 1.2e-15 for k=2, and exact at N = k.
`H_k` contains V only for k >= 2, which is why rank 1 degenerates.

## 3. The trace identities

Definitional plus the contraction sum rules, no dynamics.

    sum_x lambda_2[px,qx] = -(gamma_1 - gamma_1^2) = -gamma_1 eta_1

verified at 1e-16 (`tools/lambda2_contraction.py`). The contraction of lambda_2 IS the
particle-hole blocking product, and since `gamma_1 eta_1 >= 0` exactly when
`0 <= gamma_1 <= 1`, the identity is the Pauli bound written in cumulant variables.

At rank 3, with A the gamma_1^3 orbit and B the gamma_1*lambda_2 orbit
(`tools/rank3_trace_identity.py`, verified at 1e-16 for (4,2), (6,3), (8,4), (10,4)):

    sum_y A = N (g1^g1) - (g1^g1^2) - (g1^2^g1)
    sum_y B = N lam2 + [four gamma_1 contractions] + (g1^G) + (G^g1),   G = -(g1 - g1^2)

and the whole thing collapses (`tools/lambda3_factorization.py`, verified at 1e-16 for
(4,2), (6,3), (8,4), (10,4), (10,5)) to

    sum_y lambda_3[pvy,qwy] = E lambda_2 + lambda_2 E ,    E = gamma_1 - 1/2 = (gamma_1 - eta_1)/2

with E acting on each of the four slots; in pair-space form `(E (+) E) = Gam - 1 = (Gam - Eta)/2`.

This is NOT a product of geminal blocking factors, and the difference is structural. `-g1 eta`
vanishes at both boundaries, empty and full, which is what makes it a blocking factor.
`E = gamma_1 - 1/2` is maximal at the boundaries and vanishes at half filling. So the rank-3
trace carries no blocking of its own; blocking reaches it only through lambda_2, which must
itself vanish at saturation by the rank-2 identity.

## 4. Where the blocking is, and what breaks it

`tools/blocking_decomposition.py`. The rank-1 EOM splits into 8 two-body terms carrying no
lambda_2 plus the 2 one-body terms, which are exactly TDHF and move NO occupation at any time
(max |ndot| 7e-16 throughout a propagation), against 4 terms carrying lambda_2 which supply
all occupation motion (1.1 to 2.7). The blocking factors are present and work perfectly; they
simply cannot constrain lambda_2, which is an independent tensor rather than a product of
gamma and eta.

`tools/trace_identities.py`. For the residual
`R = sum_x lambda_2[px,qx] + (g1 - g1^2)`, whose rate vanishes on exact dynamics:

    lambda_3 used         |Rdot|
    exact                 3.0e-15
    zero                  3.7e-01
    MK_gamma30, -(A+B)    1.7e+00

    split with exact lambda_3:  no-lambda_3 terms 3.7e-01,  lambda_3 terms 3.9e+00

The lambda_3 terms are not what breaks the identity, they are what preserves it. The 132
terms without lambda_3 break it on their own, and the 8 lambda_3 terms cancel them to 3e-15.
lambda_3 restores the lower-rank blocking.

`MK_gamma30` at N=3 preserves D, since gamma_2-dot remains the isospectral commutator on the
rank-6 face, but Q reaches -1.96, G reaches -0.69, and gamma_1 breaks Pauli at -0.233
(`tools/mk_gamma30_dqg.py`). The cumulant contraction goes negative at t = 0.10 while gamma_1
is still admissible, and gamma_1 follows at t = 0.20: the contraction is the leading
indicator.

## 5. The open conjecture

lambda_3 enters `Rdot` linearly. So `Rdot = 0` is a linear condition on lambda_3, defining an
affine family of reconstructions that preserve the partial-trace identity exactly -- a
constraint derived from the structure rather than a projection applied after the fact or a
factor inserted anywhere. Whether that affine family also contains a member preserving
positivity is the open question, and it is an intersection problem rather than a search.
