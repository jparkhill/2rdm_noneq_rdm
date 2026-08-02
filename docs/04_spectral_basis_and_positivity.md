# Working in the spectral basis of the density matrix

JP's proposal: derive the 2-RDM EOM, then apply the spectral decomposition and
look for positivity violation there — e.g. as geminal occupation numbers going
out of range.  This note works out why that is the right move, what it makes
exact, and where it stops.

## 1. The spectral basis turns the cone condition into a diagonal element

For a Hermitian flow Ṁ = F(M) the PSD cone {M ⪰ 0} is forward invariant iff
(Nagumo / subtangentiality):

    for every M ⪰ 0 and every v with Mv = 0:   ⟨v| F(M) |v⟩ ≥ 0.

In the eigenbasis of M this is *literally the diagonal element of the EOM in
the null space*.  Writing M = U Λ U†, first-order eigenvalue perturbation gives

    λ̇_g = ⟨g| Ṁ |g⟩,        U̇ = U X,  X_{gh} = ⟨g|Ṁ|h⟩/(λ_h − λ_g)  (g≠h)

so the eigenvalue equation is *free of energy denominators* — all the denominators
live in the eigenvector rotation, which positivity does not care about.  This is
exactly why JP's intuition is right: the offending terms cannot hide.  You need
one quadratic form per null direction, and nothing else.

**Rank-1 sanity check.** At rank 1 the natural-orbital basis makes this the
statement n_p = 0 ⟹ ṅ_p ≥ 0 and n_p = 1 ⟹ ṅ_p ≤ 0 — i.e. exactly the invariant
region argument of `01_derivation_1rdm.md` §4.2.  The spectral view reproduces
the one result we have, which is the minimum bar for trusting it at rank 2.

### 1.1 Degeneracy is not optional
Near the boundary many eigenvalues are 0, so the null space is massively
degenerate and one must use *degenerate* first-order perturbation theory:
the correct object is the **block** P_ker Ṁ P_ker, and the condition is that
this matrix is PSD, not that its diagonal entries are non-negative.  Tracking
individual λ_g through a degeneracy is meaningless (and numerically explosive
through the 1/(λ_h−λ_g) denominators).  Practical consequence: propagate ²D (or
its square-root factor), and diagonalise only as a diagnostic.

## 2. This makes the HF boundary problem completely explicit

For a single determinant, γ is idempotent and

    ²D_HF = γ ∧ γ  is supported on occ⊗occ only.

So of the C(2r,2) geminals, exactly C(N,2) have nonzero occupation and the rest
are exactly zero.  Measured (see §4): for H4/STO-3G, 6 nonzero out of 28 — a
22-dimensional kernel; for H2O/STO-3G, 45 of 91.  **HF is not near ∂C, it is an
extreme point with a kernel covering most of the space.**

The kernel is the (occ⊗virt) ⊕ (virt⊗virt) sector.  Therefore:

> **Positivity of any HF-based 2-RDM perturbation theory, at the first order at
> which it acts, is exactly the statement that the ov/vv block of ²Ḋ is PSD.**

That is a finite, explicit matrix condition — no SOS program, no cone geometry.

### 2.1 It predicts which terms matter, and the old code agrees
Sandwiching by γ and η = 1 − γ is a projection onto p-h blocks: γ·X·γ is the
oo block, η·X·η the vv block.  So the boundary flux at HF is controlled *only*
by the η-sandwiched structures.  In `TightBindingTCL2/TCL_EE2.h` the assembly is

    Out ± ( Rho*I1*Eta + Eta*I2*Rho + Rho*I3*Rho + Eta*I4*Eta )

with `Eta*I4*Eta` commented **"Ex.State Correlation Terms"** — that is the vv
block, i.e. precisely the term whose sign decides whether the flow enters or
leaves the cone at HF.  And that overall sign is `+=` in `RhoDot` and `-=` in
`RhoDotGood`: the two variants differ exactly on the quantity this analysis says
is decisive.  Worth resolving before anything else is rebuilt.

## 3. Purification / Stinespring: what it buys and the precise reason it is not enough

Any PSD matrix is a partial trace of a pure state, so
"²D ⪰ 0" ⟺ "²D = Tr_anc |Ψ⟩⟨Ψ| for some ancilla".  Propagating a factor X with
²D = XX† (equivalently propagating |Ψ⟩) makes D-positivity structural: *any* Ẋ
gives ²Ḋ preserving D ⪰ 0.  Given a target ²Ḋ one solves a Sylvester equation
Ẋ X† + X Ẋ† = ²Ḋ for Ẋ; the gauge freedom X → XV (V unitary on the ancilla) is
fixed by demanding Ẋ X† Hermitian.  The natural-geminal basis *is* the SVD of X.

**The sharp statement of the limitation.**  N-representability says ²D is the
two-particle contraction of an *N-fermion* state.  Purification says ²D is the
contraction of a two-particle-plus-ancilla state with *no antisymmetry between
system and ancilla* and no constraint that the ancilla live in the same
one-particle space.  The whole gap between D ⪰ 0 and N-representability is
exactly that dropped antisymmetry.  Re-imposing it fully = re-inventing the
wavefunction and losing polynomial scaling.  So purification is not a loophole;
the interesting question is what a *partially* antisymmetrised, small ancilla
buys — that is an open and, as far as the literature sweep found, unclaimed
direction.

### 3.1 Where the particle-hole conjecture cashes in
D and Q are related by an affine map (Q = D + terms linear in γ + constant), so
the natural-geminal basis that diagonalises D does **not** diagonalise Q.  But
the p-h conjugate basis does.  Hence: *if* the generator is p-h covariant, the
boundary-flux argument in the conjugate spectral basis gives Q for free from the
D argument.  This is the operational form of JP's conjecture, and it is a Lean
target.  G is a genuinely third decomposition (mixed p-h) and still needs its
own argument.

## 4. Measured: is the reference we perturb about even representable?

`scripts/check_reference_representability.py`, STO-3G, spin-orbital DQG:

| system | reference | min eig D | min eig Q | min eig G |
|---|---|---|---|---|
| H4 linear 1.0 Å | HF   | 0 (exact) | 0 (exact) | −2e−16 |
| | MP2  | **−3.32e−2** | **−5.12e−2** | **−1.08e−2** |
| | FCI  | +2.7e−6 | +5.3e−6 | −1.5e−16 |
| H4 square 1.23 Å | MP2 | **−3.96e−2** | **−5.77e−2** | **−8.60e−3** |
| | FCI  | −2.4e−16 | +2.2e−16 | −4.0e−16 |
| H2O STO-3G | MP2 | **−7.83e−2** | **−9.53e−2** | **−8.27e−2** |
| | FCI  | +7.4e−4 | +1.6e−14 | −5.1e−16 |

FCI passes to machine precision (validates the code); HF sits exactly on the
boundary (validates §2).  **The MP2 2-RDM is badly non-representable** — for
H2O, 69 of 91 Q eigenvalues are negative, total negative mass 0.67.

This matters because EE2's `CorrelateInitial()` offers `Mp2NOGuess` as an
initial condition.  Any run started there began outside the cone, and no
property of the generator can rescue that.  Screening the reference is a
one-line check that should precede blaming the EOM.  Note also that Q is the
worst-violated in every case: MP2's failure is concentrated in the **hole**
sector, which is suggestive given the p-h theme.
