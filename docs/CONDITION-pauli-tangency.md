# The Pauli Tangency Condition

## The algebraic result

The MKNO rank-1 equation of motion is exact in `(gamma_1, lambda_2)`:

    gamma_1-dot = i [f, gamma_1] + C[lambda_2]

`f[p,q] = h[p,q] + sum_{r,s} V[p,r,q,s] gamma_1[r,s]` is the generalised Fock operator: the
one-body Hamiltonian plus the mean field built from the one-particle reduced density matrix
`gamma_1[p,q] = <a^dag_p a_q>`. `C` is the collision term, defined as exactly those terms of
`gwt.terms_cumulant(1)` carrying a two-particle cumulant factor `("L",2)`.

There are exactly four such terms:

    +1  V[r,t,s,u] d[p,s] l2[r,t,q,u]        -1  V[r,t,s,u] d[p,u] l2[r,t,q,s]
    -1  V[r,t,s,u] d[r,q] l2[p,t,s,u]        +1  V[r,t,s,u] d[t,q] l2[p,r,s,u]

and by the antisymmetry of `V` they collapse in pairs, giving the closed form

    C[p,q] = (i/2) ( X - X^dagger )[p,q],    X[p,q] = sum_{r,t,u} V[r,t,p,u] lambda_2[r,t,q,u]

verified in `tools/collision_form.py` to `5.4e-15` against the assembled term list, on a generic
antisymmetric Hermitian `lambda_2` (so it is an identity, not a coincidence of some state).

**And `|C(gamma_1) - C(2 gamma_1)| = 0` exactly.** The collision term contains no explicit
`gamma_1` whatsoever; it is a fixed bilinear map from `(V, lambda_2)` to a Hermitian matrix.

## What this rules out

A fermionic Lindblad master equation with loss channels (jump operators `a_p`) and gain channels
(`a^dag_p`) yields a 1-RDM equation of the gain-loss form

    gamma_1-dot = i[f, gamma_1] + (1/2){Gam_plus, eta_1} - (1/2){Gam_minus, gamma_1}

with `eta_1 = 1 - gamma_1` the one-hole density matrix and `Gam_plus`, `Gam_minus` positive
semidefinite. That form makes `0 <= gamma_1 <= 1` automatic: at an empty orbital the rate is
`<phi|Gam_plus|phi> >= 0`, at a full one `-<phi|Gam_minus|phi> <= 0`. It was the candidate
"Lindblad analogue" -- a parameterisation making positivity structural rather than monitored.

It cannot be realised here. `C` has no `gamma_1` in it, so there is no `eta_1` to write. One can
always *define* `Gam_plus`, `Gam_minus` to absorb a `gamma_1` dependence, since the decomposition
is underdetermined, but then the form is not manifest and guarantees nothing. **This route is
closed.**

## The condition that replaces it

Since `C` is a fixed map, all Pauli blocking must come from the closure's dependence of
`lambda_2` on `gamma_1`. The commutator `i[f, gamma_1]` is unitary and cannot move the spectrum of
`gamma_1`, so forward invariance of the Pauli cone is entirely a condition on `C`. With
`<phi|C|phi> = -Im <phi|X|phi>`:

    Pauli Tangency Condition (PTC)

        Im <phi|X|phi>  <= 0   for every natural orbital with gamma_1|phi> = 0     (empty)
        Im <phi|X|phi>  >= 0   for every natural orbital with gamma_1|phi> = |phi> (full)

        <phi|X|phi> = sum_{r,t,u} V[r,t,phi,u] lambda_2[r,t,phi,u]

the boundary orbital being contracted into the third slot of both `V` and `lambda_2`. This is
Nagumo tangency specialised to the Pauli cone and reduced to a sign condition on the imaginary
part of a single `V`-`lambda_2` contraction. It is checkable for any candidate closure.

## Why "tangency"

For the exact dynamics the boundary rate is exactly zero, not merely correctly signed: `n_p = 0`
implies `a_p|psi> = 0`, hence `gamma_1[p,:] = 0` and every `lambda_2` component with `p` in a
lower index vanishes, so `<phi|X|phi> = 0`. The exact flow is tangent to the Pauli boundary, and
the orbital is populated only at second order in time through `lambda_2`'s own evolution.

**Corollary, and a warning about adiabatic closures.** If `lambda_2 = F[gamma_1]` vanishes
identically whenever `n_p = 0`, then `n_p = 0` is a forward-invariant manifold of the closed
dynamics and the orbital stays empty forever -- which the exact dynamics does not do. Exact
boundary behaviour and correct population transfer are in tension for *any* adiabatic closure.
This is a property of adiabatic elimination rather than of the blocking structure, and it predicts
systematic under-population of weakly occupied orbitals. It is the same degeneracy already seen as
`R = 0` at rank 1 and as the 100-dimensional `ker L0` at rank 2, in a third guise.

## The Lipschitz screen, which should be applied before any of this

Nagumo's theorem, existence and uniqueness for the truncated ODE, and the frozen-coefficient
argument that would carry Lindblad-type guarantees to a nonlinear propagator all require the
closure to be locally Lipschitz in its arguments. Both families investigated so far fail:

- `pt = -L0^{-1}[D]`: `L0` depends on `gamma_1` through `f`, and its eigenvalues
  `i(eps_p + eps_v - eps_q - eps_w)` cross zero as `gamma_1` varies. At a crossing the
  pseudo-inverse's norm diverges and its rank jumps, so the closure is discontinuous, not merely
  non-Lipschitz. Warning signs were already measured: `rank(L0) = 684` of `784` with
  `sigma_min>0 = 3.76e-03`.
- `blocked`: pole surface `a_i = (1/2) sum_j a_j` strictly inside the feasible set.

What survives the screen: polynomials in `gamma_1` and `V`; truncated Neumann or Chebyshev series
in place of a full inverse; and level-shifted inverses such as `(L0^2 + sigma^2)^{-1} L0` with a
FIXED `sigma > 0`, Lipschitz for every positive `sigma`, recovering the pseudo-inverse only in the
singular limit. The regularisation parameter is part of the closure's definition, not a numerical
detail.

## Status of the numerical test in tools/collision_form.py

The block test built into that file is **vacuous as written** and is retained only as a record.
It freezes spin-orbitals 0,1 occupied and 6,7 empty, which makes them non-interacting spectators;
by extensivity `lambda_2` then has no components touching them, `C` has zero rows and columns
there, and both blocks come out exactly zero rather than correctly signed. A meaningful test needs
a boundary orbital that is coupled, and exact emptiness fights coupling -- so the right experiment
is the *approach* to the boundary, `n_p -> 0`, measuring the sign of `Im <phi|X|phi>` and the rate
at which it vanishes. Not yet run.
