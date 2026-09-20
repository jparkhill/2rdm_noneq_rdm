# Rigorous derivation of the bath dissipator for the 2-RDM hierarchy

**Status: DERIVATION, not run.** Every approximation is named where it is made. The places where a numerical
check is required before use are marked CHECK.

## 0. Why the rank-1 result does not transfer trivially

The published rank-1 theory (JCP 142, 134113) derives a TCL2 master equation for the *system density operator*
ρ, and then takes its one-body contraction; generalised normal ordering with respect to the correlated γ₁ puts
the blocking factors `n(1−n)` on the rates, which makes `0 ≤ γ₁ ≤ 1` invariant. **The same object — the TCL2
generator on ρ — contracts to every rank.** So the correct procedure is not "add a dissipator to the rank-2
equation"; it is: derive the dissipator ONCE at the level of ρ, then contract it to rank 1 and rank 2 and ask
which higher RDMs it introduces. That is what §2–§4 do. The EE2 failure was that the rank-2 contraction was
never done consistently with the rank-1 one.

## 1. Setup

System: `H_S = Σ h_pq a†ₚa_q + ¼ Σ V_pvqw a†ₚa†ᵥa_w a_q` on M spin-orbitals.
Bath: harmonic, `H_B = Σ_α ω_α b†_α b_α`, in a thermal state `ρ_B ∝ e^{−βH_B}`.
Coupling, linear in bath coordinates and **one-body** in the system:

```
H_SB = Σ_α  A_α ⊗ (b_α + b†_α),        A_α = Σ_pq g^α_pq a†ₚ a_q      (Hermitian one-body system operators)
```

This is the standard electron–phonon / system–boson coupling. (A two-body system coupling would give rank-3
terms at the rank-2 level; it is excluded here, as it was in the rank-1 paper.) Bath correlation functions
`C_αβ(t) = tr_B[ ρ_B B_α(t) B_β(0) ]`, `B_α = b_α + b†_α`, with `C_αβ(t) = δ_αβ C_α(t)`,
`C_α(t) = coth(βω_α/2) cos ω_α t − i sin ω_α t`, and spectral density `J(ω) = Σ_α g²_α δ(ω − ω_α)` once a
continuum is taken.

## 2. TCL2 generator on ρ (exact through second order in H_SB)

Time-convolutionless projection with `P ρ_tot = ρ_S ⊗ ρ_B`, second order:

```
d/dt ρ = −i[H_S, ρ] + D[ρ],
D[ρ]   = −Σ_α ∫₀^∞ ds  { C_α(s) [A_α, Ã_α(−s) ρ] − C_α(−s)... }      → in standard form
D[ρ]   = Σ_α ( [Ã_α ρ, A_α] + [A_α, ρ Ã_α†] ),        Ã_α = ∫₀^∞ ds C_α(s) e^{−iH_S s} A_α e^{iH_S s}
```

**Approximation A1 (Born / second order in g).** **Approximation A2 (Markov: upper limit ∞).** Both are the same
as in the rank-1 paper. The result is a Redfield generator: Hermiticity- and trace-preserving, NOT completely
positive in general.

**The obstacle that is specific to RDMs.** `Ã_α` requires `e^{−iH_S s} A_α e^{iH_S s}`, the Heisenberg evolution
of a one-body operator under the interacting `H_S`. For a one-body `H_S` this stays one-body (u(M) closes); for
the interacting `H_S` it develops two-, three-,… body components. **Every truncation of the bath term for an
interacting system is a truncation of this Heisenberg evolution**, and that is where the closure enters.

## 3. Closing `Ã_α` with the effective one-body generator (the choice that makes contraction consistent)

**Approximation A3.** Replace `e^{−iH_S s} A e^{iH_S s}` by `e^{−i f s} A e^{i f s}` with `f = h + W[γ₁]` the
instantaneous mean-field (Fock) one-body operator. Then `Ã_α` stays one-body:

```
Ã_α = Σ_pq  g̃^α_pq a†ₚ a_q,        g̃^α = Σ_ab  |a⟩⟨a| g^α |b⟩⟨b|  Γ_α(ε_b − ε_a)
Γ_α(ω) = ∫₀^∞ ds C_α(s) e^{iωs} = π J-type half-Fourier transform:  Γ(ω) = γ(ω)/2 + i Δ(ω)
γ(ω) = 2π J(|ω|) [ n_B(|ω|) + θ(ω) ]     (emission for ω>0, absorption for ω<0; n_B Bose function)
```
where `|a⟩` are eigenvectors of `f` with eigenvalues `ε_a`. This is the rank-1 paper's choice, evaluated at the
CURRENT γ₁ (not at a fixed reference), and it is what makes the rank-1 contraction come out with `n(1−n)` factors.
Its justification is the same as the mean-field justification of TDHF: it is the exact Heisenberg evolution of a
one-body operator under the *one-body part* of the generator, and it neglects the two-body scattering of the
coupling operator during the bath memory time. **It is an approximation and it is where the two-body content of
the system enters the bath term only through `f`.** A more accurate alternative — evolve `A` under the effective
PAIR Hamiltonian of episode 5 — is discussed in §6.

With A3, `D[ρ]` is a **quadratic Lindblad-type generator**: `A_α` and `Ã_α` are one-body, so `D` is a sum of
terms `a†a ρ a†a`. This is decisive for §4: a quadratic dissipator, contracted to rank k, closes on ranks
≤ k+1 exactly and — after GNO — on rank k plus products of lower ranks.

## 4. Contraction to the hierarchy: what the dissipator does at every rank

For any operator O, `d/dt ⟨O⟩_D = tr(O D[ρ]) = Σ_α ⟨ [A_α, O] Ã_α + Ã_α† [O, A_α] ⟩` (cyclicity). With
`A = Σ g a†a`, `Ã = Σ g̃ a†a`:

**Rank 1** (`O = a†ₚa_q`): `[A, a†ₚa_q]` is one-body, times the one-body `Ã` gives ⟨one-body · one-body⟩ = a
TWO-body expectation. So `d/dt γ₁|_D` is linear in γ₂:

```
d/dt γ₁[p,q]|_D = Σ_α Σ_rs ( g^α[q,r] ⟨a†ₚ a_r Ã_α⟩ − g^α[r,p] ⟨a†_r a_q Ã_α⟩ + h.c.-type )
                = (terms linear in γ₁) + (terms linear in γ₂, via ⟨a†a a†a⟩ = γ₂-contraction + δ·γ₁)
```
Writing `γ₂ = γ₁∧γ₁ + λ₂`: the `γ₁∧γ₁` part produces exactly the Pauli-blocked rates
`Σ γ(ε_b−ε_a) |g_ab|² [ n_b(1−n_a) − n_a(1−n_b) ]`-type of the rank-1 paper (the hole factors ARE the wedge term),
and the `λ₂` part is a **new term, linear in λ₂, that the rank-1 paper set to zero** (it propagated γ₁ only).
CHECK 1: reproduce the published rank-1 dissipator from this contraction at λ₂ = 0, term by term.

**Rank 2** (`O = a†ₚa†ᵥa_w a_q`): `[A, O]` is two-body, times `Ã` gives ⟨two-body · one-body⟩ = a THREE-body
expectation. So

```
d/dt γ₂|_D  =  (linear in γ₂, with g, g̃)  +  (linear in γ₃, with g g̃)
```

**This is the point that EE2 missed and that must be handled.** The bath term at rank 2 carries γ₃, exactly as
the Hamiltonian term does. Substituting the MKNO expansion `γ₃ = A₃(γ₁γ₁γ₁) + A₃(γ₁λ₂) + λ₃`:
- the `γ₁γ₁γ₁` and `γ₁λ₂` parts are polynomial in `(γ₁, γ₂)` and produce, at rank 2, the pair-space analogue of the
  Pauli blocking: rates dressed by `(1 − n)` factors on the pair indices — the same structure as `K_V(I − N_op)`
  in F0, now on the bath transition amplitudes;
- the `λ₃` part is dropped (**Approximation A4**, the same λ₃ = 0 truncation as in F0, applied to the dissipator).

Because `D` is quadratic (one-body × one-body), its `γ₃` dependence is only through the SINGLE contraction
`Σ_s g̃[s,·] γ₃[..s..]`, so after A4 the rank-2 dissipator is a finite polynomial in `(γ₁, γ₂, g, g̃)` — closed-form
identifiable in pair/ph space exactly as `Phi_ggg`, `Phi_gl` were. CHECK 2: build it from the Wick term list of
`tr(O D[ρ])` (same machinery as `qdiss.eom.wick`, with the two one-body operators in place of H's terms), verify
against the exact dissipative rate computed on the full Fock space for a small system.

**Rank k in general**: `d/dt γ_k|_D` is linear in `γ_k` and `γ_{k+1}`, never higher — because `D` is quadratic.
So the bath enlarges no rank's coupling beyond what the Hamiltonian already couples. **The dissipative hierarchy
has the same BBGKY structure as the reversible one, one rank up at each level.** Whatever closure is used for the
Hamiltonian part must be used for the dissipative part at the same rank, or the two are inconsistent (this is
exactly the inconsistency in EE2, where γ₂ was propagated with a rank-1-derived dissipator).

## 5. Properties of the closed rank-2 dissipator

With A1–A4 the rank-2 equation becomes `d/dt X = F0(X) + δ + D₂(X)` where `D₂` is:
- linear in `g g̃`, polynomial in `(γ₁, γ₂)`;
- **trace-preserving** on pairs (CHECK 3: `tr D₂ = 0` identically, follows from `tr D[ρ] = 0` provided the
  γ₃ → separable substitution preserves the partial trace, which `A₃(γ₁γ₁γ₁)+A₃(γ₁λ₂)` does NOT exactly — the
  trace defect is `Tr₃λ₃`-sized and must be measured; if nonzero, the contraction constraint (b) in the δ
  machinery absorbs it);
- **contraction-consistent**: `Tr₂ D₂/(N−1)` must equal the rank-1 dissipator of §4 evaluated with the same γ₂;
  this is automatic if both come from the same `D[ρ]` and the same A4 (CHECK 4: verify numerically).
- **not CP as written** (Redfield). Positivity of D, Q, G along the dissipative flow is therefore again the job
  of the δ correction, OR of a secular/Davies approximation (**Approximation A5**, optional): keep only
  `g̃_ab` with `|ε_a − ε_b|` resolved, which turns `D` into GKSL form on ρ. GKSL on ρ contracts to
  positivity-preserving flows at every rank (the contraction of a CP map is a well-defined RDM of a legitimate
  state), so with A5 **the dissipative part is inward-pointing at all ranks by construction** and δ is only needed
  for the Hamiltonian part's λ₃ error.

**Fixed point.** With A5 and detailed balance `γ(−ω)/γ(ω) = e^{−βω}`, the fixed point of `D[ρ]` on ρ is the Gibbs
state of the *effective one-body* Hamiltonian `f` (Davies' theorem applied to a quadratic generator) — a
fermionic Gaussian state with Fermi–Dirac occupations of `f`'s eigenvalues, NOT the FCI thermal state. The
Hamiltonian part `F0` then competes: the combined fixed point is where `F0(X) + δ + D₂(X) = 0`, a correlated
state that is neither Gaussian nor exactly FCI-thermal. **How close it is to the FCI thermal γ₂ at temperature
T is the quantitative question for episode 10**, and it is a property of the closure, not of the bath: with the
exact hierarchy the fixed point would be exactly the FCI Gibbs state (the exact `D[ρ]` with a full Heisenberg
`Ã` has that as its fixed point under detailed balance).

## 6. Beyond A3: evolving the coupling operator under the pair Hamiltonian

A3 evolves `A_α` under `f` (one-body). Episode 5 gives the exact one-body-plus-pair-space generator of the λ₃=0
flow; evolving `A_α` under it, `Ã_α` acquires a two-body part `Σ g̃₂ a†a†aa`, and then `D` contains
`a†a†aa ρ a†a†aa`-type terms. Contracted to rank 2 this brings in γ₄; closed with λ₄ = 0 (consistent with the
hierarchy's own truncation at the next rank) it is again polynomial. This would make the dissipator's fixed
point closer to the correlated thermal state at the price of one more rank of algebra. It is the natural
second step, not the first.

## 7. Plan for episode 10 (in order)

1. Build `D[ρ]` on the full Fock space for H4 with a chosen `J(ω)` (Ohmic, `J = η ω e^{−ω/ω_c}`), one coupling
   operator (e.g. the x-dipole), A1–A3 as above, and propagate ρ exactly → the "exact-within-Born-Markov"
   reference RDMs. Verify Hermiticity/trace; measure its fixed point vs the FCI Gibbs state of `H_S` (they differ
   because of A3; that difference is the bath model's own error and is separate from the closure's).
2. Derive `D₁` and `D₂` by Wick contraction of `tr(O D[ρ])`; CHECK 1–4.
3. Identify `D₂` in closed pair/ph form (as for `Phi`), with the `(1−n)` structure made explicit.
4. Propagate `F0 + δ + D₂` from the HF determinant at `kT = 0.5` and `0.2`; report: dS/dt sign, cones, distance of
   the endpoint γ₂ from (a) the reference of step 1 and (b) the FCI Gibbs γ₂; then kick the endpoint and compare
   the dynamics with the kicked reference. **System, temperatures, propagation length in fs, and all four
   approximations A1–A4 (A5 if used) stated in every table.**
