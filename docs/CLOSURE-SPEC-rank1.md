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

## 4. The closure -- and why rank 1 is the HARD case, not the easy one

**REVISED.** The earlier version of this file used an MP2-type amplitude
`T_sigma = V / (orbital energy denominator)`. **That is withdrawn.** Its order relation is anchored to
an idempotent reference; away from one it carries no justification at all, and during a propagation
far from the initial state there is no order relation to the equation of motion. It must not be used
as a closure ingredient.

What remains is a genuine obstruction, and it is worth stating plainly rather than papering over.

**Leg counting.** `Tr_2 lambda_2 = -gamma_1 eta_1` is the trace target: a `(1,1)` object. `lambda_2`
needs `(2,2)`. `gamma_1` vertices cannot change the leg count (a connected graph of `(1,1)` vertices is
a path or a cycle). So the structural-trace construction of `STRUCTURAL-TRACE.md` needs exactly one
`(2,2)` vertex attached by one contraction:

    lambda_2  =  A_2[ (-gamma_1 eta_1)  o  K ]  ,      Tr_2 K = I

**But the only `(2,2)` primitives are `lambda_2` and `V`.** Using `lambda_2` is circular -- it is the
object being constructed. So `K` must be built from `V`, and since `lambda_2` is dimensionless while
`V` is an energy, `K` requires an energy denominator.

**This is where MP2 came from, and the objection is to the JUSTIFICATION, not to the resolvent.** A
denominator built from orbital energies of an idempotent reference is illegitimate. A denominator built
from the INSTANTANEOUS one-body generator of the equation of motion is a different object: it is defined
at any density, it is not a perturbative expansion, and it claims no order relation. Concretely, with
`f = h + W[gamma_1]` the generalised Fock matrix at the current density and `L_0` the superoperator
`L_0[X] = i(f-commutator on all four slots of X)`, a candidate is

    K  =  ( L_0^dag L_0 + sigma^2 )^{-1} L_0^dag [ V ]        (Tikhonov-regularised resolvent)

normalised so that `Tr_2 K = I`. **This is a PROPOSAL and is not derived.** It is dimensionally sound
and Lipschitz for fixed `sigma > 0` (`tools/level_shift.py`: the Lipschitz constant falls from `1.6e3`
to `0.85` at `sigma = 1e-2` at no measurable cost in accuracy), but:

- I have NOT shown that `Tr_2 K = I` is achievable for this `K` without reintroducing a global
  normalisation, which would break extensivity;
- I have NOT shown the construction is extensive; `L_0` is built from `f`, which is local, so it is
  plausible, but plausible is not measured;
- the `sigma -> 0` limit is where `L_0` is singular, which is exactly the rank-jumping pseudo-inverse
  that disqualified the unregularised `pt` closure under screen 1.

## 5. Honest assessment of the rank-1 scheme

**The rank-1 closure is not ready and should not be implemented.** The rank-2 scheme is, because there
`lambda_2` is a propagated state variable and can serve as the `(2,2)` vertex without circularity. At
rank 1 there is no such vertex, and every route to one passes through `V` with a denominator, which is
where the previous generation of closures (`pt`, `mp2`) failed.

This is not a defect of the specification; it is the rank-2 no-go asserting itself. It says the cheapest
rung of the hierarchy is the hardest to close honestly, and that the natural entry point is rank 2.

## 6. Ledger status

| row | met? | how |
|---|---|---|
| `0 <= gamma_1 <= I` | YES, structural | `Theta` isometry, section 2 |
| antisymmetry, hermiticity | YES | representation |
| dimensional analysis | YES for the resolvent form; NO for bare `V` | section 4 |
| **trace relation** | **UNRESOLVED** -- needs a `(2,2)` vertex that is not `lambda_2` | section 4 |
| **extensivity** | **UNRESOLVED** -- depends on the normalisation of `K` | section 4 |
| Lipschitz | plausible for fixed `sigma`; `sigma -> 0` disqualified | `tools/level_shift.py` |
| exactness at `k >= N` | YES for `N <= 2` via `gamma_3 = 0` | `tools/top_of_hierarchy.py` |
| perturbative limit | DELIBERATELY NOT CLAIMED -- see section 4 | -- |
| convergence | untested | -- |

Three rows are structural, two are unresolved, and the unresolved ones are the ones that matter.

## 7. What is owed before this rung is attempted

1. Decide whether a `(2,2)` vertex with `Tr_2 K = I` can be built from `V` and the instantaneous `f`
   WITHOUT a global normalisation. This is the whole question; everything else follows.
2. If yes, measure its extensivity on the H3 x H2 product-state harness
   (`tools/ansatz_extensivity.py`) before any dynamics.
3. Only then specify the coefficient flow. It is deliberately omitted here because specifying a flow
   for a closure whose leading term is unresolved would be exactly the getting-ahead-of-the-algebra
   that this revision is meant to stop.
