# Closure specification, rank 1: propagate gamma_1

**Status: PROPOSAL, not run.** Nothing in this file has been propagated. Every approximation is
named in section 6 and every unvalidated step is marked UNVALIDATED.

## 1. Symbols

| symbol | definition |
|---|---|
| `M` | number of spin orbitals; indices `p,q,r,s,t,u,v,w,x,y` run `1..M` |
| `N` | number of electrons (sharp, not averaged) |
| `h[p,q]` | one-body integrals in the convention of `qdiss.eom.hamiltonian` |
| `V[p,t,q,u]` | antisymmetrised two-body integrals, same convention; `H = sum h a+ a + (1/4) sum V a+ a+ a a` |
| `gamma_1[p,q]` | `<a+_p a_q>`, `tr gamma_1 = N` |
| `eta_1` | `I - gamma_1` |
| `gamma_2[p,v,q,w]` | `<a+_p a+_v a_w a_q>`, `tr gamma_2 = N(N-1)` |
| `lambda_2` | `gamma_2 - gamma_1 wedge gamma_1`, the 2-cumulant |
| `a wedge b` | `a[p,q]b[v,w] - a[p,w]b[v,q]` |
| `W[p,q]` | `sum_{r,s} V[p,r,q,s] gamma_1[r,s]`, the mean field |
| `E` | `gamma_1 - I/2` |
| `Mat_1` | `gamma_1` itself (the `k=1` case of `Mat_k(T)[I,J] = k! T[I+J]`) |
| `T_sigma` | the amplitude defined in section 4 |
| `<A,B>` | `sum` over all indices of `conj(A) * B` (Hilbert-Schmidt) |

## 2. Propagated state

Two objects.

**(a) An isometry** `Theta` of shape `(r + rbar) x M` with `Theta^dag Theta = I_M`, split as
`Theta = [R ; Rbar]` with `R` of shape `r x M` and `Rbar` of shape `rbar x M`. Then

    gamma_1 = R^dag R ,        eta_1 = I - gamma_1 = Rbar^dag Rbar

so `0 <= gamma_1 <= I` -- BOTH Pauli bounds -- hold by construction. This is the rank-1 analogue of
the XR construction; one Stiefel condition does both bounds. `r + rbar >= M` is required;
`r = rbar = M` is the safe choice.

**(b) Coefficients** `c_a`, `a = 1..A`, complex.

## 3. Exact equation of motion being approximated

The rank-1 MKNO equation of motion, generated symbolically in `qdiss.eom.wick` (2 one-body and 4
two-body terms), term-identical to the MKNO route `qdiss.eom.gwt`, and validated against exact
Fock-space linear algebra `qdiss.eom.fullci` at 1e-15 by `qdiss.eom.validate`. In closed form,
established in `tools/collision_form.py` at 5.4e-15:

    gamma_1-dot  =  i [ (h + W)^T , gamma_1 ]  +  C
    C[p,q]       =  (i/2) ( X - X^dag )[p,q]
    X[p,q]       =  sum_{r,t,u} V[r,t,p,u] lambda_2[r,t,q,u]

The first term is TDHF and is exact given `gamma_1`; verified to give TDHF at 3.9e-16 when
`lambda_2 = 0`. **`C` contains no explicit `gamma_1`**, so the entire closure content sits in
`lambda_2`.

## 4. The closure

`lambda_2` is NOT a functional of `gamma_1`. Two facts force this:

- **No connected `gamma_1`-only structure with 4 free indices exists.** A connected diagram of
  `gamma_1` vertices (each `(1,1)`) is a path or a cycle and always has `(1,1)` or `(0,0)` free
  indices. Proof in `CONSTRAINT-LEDGER.md`. So any `lambda_2[gamma_1]` is a wedge and LEAKS
  (measured: `muller_rdmft` leak 0.99).
- **Dimensional analysis** then forces the `(2,2)` primitive to be `V` divided by an energy.

So the basis is amplitudes. With `f = h + W` the generalised Fock matrix, `f = U e U^dag` its
eigendecomposition, and a FIXED level shift `sigma > 0`:

    T_sigma[p,v,q,w]  =  antisym[ V~[p,v,q,w] / ( e_q + e_w - e_p - e_v + sigma ) ]

where `V~` is `V` in the `U` basis and `antisym` antisymmetrises in `(p,v)` and in `(q,w)`; the
result is rotated back. Then

    lambda_2  =  sum_a c_a S_a ,     S_a in { T_sigma_a , dressings of T_sigma_a by powers of gamma_1 }

`S_a` is orthonormalised by SVD to a frame `F` (rows orthonormal), so `c` is unique; `c_a = <F_a, lambda_2>`
and `lambda_2 = F^T c`. The frame is aligned to the previous step's by the orthogonal Procrustes
solution, so `c` measures the state moving and not the frame being re-chosen.

## 5. The coefficient flow

Writing `lambda_2 = sum_a c_a S_a + rho` with `rho` orthogonal to the span, the exact rank-2
equation of motion gives `lambda_2-dot`, and projecting onto the frame:

    c-dot  =  <F , lambda_2-dot>  -  <F , sum_b c_b S_b-dot>

`S_b-dot` is computed by the product rule over the structure definition -- analytically, via
`qdiss.ansatz.connected_basis_dot` for the chain family; the amplitude family needs the analogous
routine, which **does not exist yet (UNVALIDATED, NOT WRITTEN)**.

`lambda_2-dot` comes from the rank-2 MKNO equation of motion, which requires `gamma_3`. The tower is
terminated by

    gamma_3 = 0      when N <= 2 (EXACT there, since gamma_3 == 0 identically)
    gamma_3 = separable(gamma_1, lambda_2)  otherwise, i.e. lambda_3 = 0   (APPROXIMATION)

The second is the `null` closure at rank 3 and is measured to be wrong by 10-19% of the flow at
`k = N` (`tools/top_of_hierarchy.py`). Note the termination is at the MOMENT level where it can be:
`gamma_{k+2} = 0` is exact when `N <= k+1`, whereas `lambda_{k+2} = 0` is not.

## 6. Every approximation, named

1. `lambda_2` is restricted to the `A`-dimensional amplitude span. Unmeasured at rank 2. (At rank 3
   the analogous connected span held `lambda_3` to 0.8-8% statically, 16% under quench, 83% from a
   rank-2 mixture.)
2. The coefficient flow is a Galerkin projection; the out-of-span residual `rho-dot` feeds back and
   is discarded. At rank 3 this leakage measured 1-3% under quench, 28-35% at low `rank(rho)`.
3. The tower is terminated at `gamma_3` as in section 5.
4. The level shift `sigma` is fixed, not variational. Measured at rank 2 elsewhere to cost nothing in
   accuracy while passing the Lipschitz screen (`tools/level_shift.py`, `sigma = 1e-2`).
5. `Theta` is integrated with a retraction (QR or polar); the retraction error is `O(dt^2)` and
   **unmeasured**.

## 7. Ledger status

| row | how it is met | evidence |
|---|---|---|
| `0 <= gamma_1 <= I` | structural, `Theta` isometry | construction |
| extensivity | `T_sigma` is built from `V`, which vanishes between separated fragments | NOT YET MEASURED at rank 2 |
| dimensional analysis | `V` over an energy | construction |
| Lipschitz | fixed `sigma`, no state-dependent inverse | `tools/level_shift.py` |
| trace identity `Tr_2 lambda_2 = -gamma_1 eta_1` | **NOT MET.** See section 8 | -- |
| exactness at `k >= N` | `gamma_3 = 0` termination is exact for `N <= 2` | `tools/top_of_hierarchy.py` |
| convergence | in `A` and in the tower depth | NOT MEASURED |

## 8. The trace identity is NOT enforceable here, and this corrects an earlier claim

I previously said the trace identity becomes "one linear condition on a handful of numbers, elementary
to impose". **That is wrong.** With `lambda_2 = sum_a c_a S_a`, define `D_a = Tr_2 S_a` (an `M x M`
object) and `Tgt = -gamma_1 eta_1`. The identity is

    sum_a c_a D_a  =  Tgt

which is `M^2` scalar conditions on `A` unknowns -- at `M = 8`, 64 conditions on `A ~ 6`. Massively
overdetermined. Differentiating gives the same count on `c-dot`. So it can only be

  (a) imposed in least squares, leaving a residual that must be reported; or
  (b) satisfied exactly by enlarging `A` until `{D_a}` spans the required subspace, at a cost that
      has not been estimated; or
  (c) violated, with `|sum_a c_a D_a - Tgt|` monitored along the trajectory.

Option (c) with monitoring is the only honest default until (b) is costed. **This is the largest
single defect of the rank-1 scheme.**
