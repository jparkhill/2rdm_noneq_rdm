# Where the rank-2 flow stops being a commutator

First measurement for route F (Hamiltonian / metriplectic closures). Tool:
`tools/commutator_structure.py`; full output `results/commutator_structure.txt`.

## The question, and the diagnostic

At rank 1 the mean-field flow is a commutator, `gamma_1-dot = i[f^T, gamma_1]`, so the spectrum of
`gamma_1` is conserved for any energy functional and `0 <= gamma_1 <= 1` is structural. The rank-2
analogue would be a closure with `Xdot = i[K, X]`, `X = Mat_2`, `K` Hermitian on pair space.

Symbols: `Mat_2[(p,v),(q,w)] = 2 gamma_2[p,v,q,w]` over sorted pairs; `gamma_2[pvqw] = <a+_p a+_v a_w
a_q>`; `x_i` eigenvalues of `Mat_2`. In the eigenbasis of `X`, `Xdot_ij = i K_ij (x_j - x_i) -
(1/2) Gamma_ij (x_i + x_j)`, so the block-diagonal of `Xdot` over the eigenspaces of `X` -- the
**spectral motion**, `x_i-dot` -- cannot come from any Hermitian `K`. The diagnostic is

    spectral fraction  =  |blockdiag_eig(Xdot)|_F / |Xdot|_F

which is 0 iff a Hermitian-`K` representation exists, and measures the distance from all of them.
**The split `Xdot = D + Phi`.** The rank-2 EOM `d/dt gamma_2[pvqw] = i<[H, a+_p a+_v a_w a_q]>` has 4
one-body and 12 two-body Wick terms, each containing EITHER `gamma_2` OR `gamma_3`. Hence exactly

    d/dt gamma_2  =  D[gamma_1, gamma_2]  +  Phi[gamma_3]

    D[pvqw]   =  i ( sum_r h[p,r] gamma_2[rvqw] + h[v,r] gamma_2[prqw] - h[r,q] gamma_2[pvrw] - h[r,w] gamma_2[pvqr] )
               + (i/2) ( sum_rs V[p,v,r,s] gamma_2[rsqw] - V[r,s,q,w] gamma_2[pvrs] )
    Phi[pvqw] =  the eight terms of the form (i/2) sum_rst V[r,s,q,t] gamma_3[p,v,t,r,s,w] and index partners

`D` ("drift") is `h` acting on each index plus `V` acting on the ket pair and on the bra pair; in
`Mat_2` language it is `i[K_2, Mat_2]` with `K_2` the Hermitian pair-space matrix `h (+) h + V`.
`Phi` ("collision part") is LINEAR in `gamma_3`, so every closure differs from the exact flow by
`Phi[gamma_3^recon - gamma_3^exact]` and nothing else. In code `D` is the term list evaluated with the
`gamma_3` slot zeroed (a term-list operation, not a moment negation used as a closure). Ground truth `qdiss.eom.fullci`; the wick term list
agrees with it at 1e-16 on every row.

States: dipole-kicked ground state (pure; and thermal `kT = 0.5` for H4); the HF determinant at
`t = 0`; an interaction quench (ground state of `h + 0.3 V` propagated exactly under `h + V`) at
`t = 0, 1, 2`, pure and thermal. H3 bent (`M=6, N=3`) and H4 trapezoid (`M=8, N=4`).

## Results

**1. The drift is EXACTLY a commutator flow.** Spectral fraction `< 1e-16` at all twelve states, both
molecules, pure and mixed. So `drift = i[K_2, X]` with `K_2` Hermitian -- the one-body Fock part plus
`V` restricted to pair space (the "pair Hamiltonian"). Every bit of spectral motion of `Mat_2` is in
`Phi[gamma_3]`. This is stronger than the earlier statement that the drift "compresses to zero on any
kernel": it conserves the whole spectrum of `Mat_2`, not just its kernel.

**2. Near equilibrium the exact flow is nearly a commutator; far from it, it is not.**

| state | rank-1 spectral frac (exact) | rank-2 spectral frac (exact total) | of `Phi` |
|---|---|---|---|
| H3 kicked ground | 0.024 | 0.016 | 0.028 |
| H4 kicked ground | 0.004 | 0.004 | 0.004 |
| H4 kicked thermal | 0.001 | 0.001 | 0.002 |
| H3 V-quench t=1 | 0.58 | 0.28 | 0.22 |
| H3 V-quench t=2 | 0.45 | 0.23 | 0.20 |
| H4 V-quench t=1 | 0.57 | 0.35 | 0.45 |
| H4 V-quench t=2 | 0.27 | 0.13 | 0.19 |
| H4 V-quench thermal t=1 | 0.48 | 0.28 | 0.27 |

Under the quench 13-35% of the exact rank-2 flow, and 27-58% of the exact rank-1 flow, is spectral.
The natural and geminal occupations genuinely move. **TDHF is not accurate because occupation motion
is small; it is stable because it is a commutator regardless.** A Hermitian-`K` rank-2 closure -- a
"geminal TDHF" -- would be positive and spectrum-conserving and would drop 13-35% of the flow in a
quench, the same kind of error TDHF makes at rank 1 (27-58% here).

**3. Time-reversal makes `t = 0` uninformative.** At the HF determinant and at the quench's `t = 0`
the exact spectral motion vanishes identically (real initial state, `x_i-dot = 0` by time reversal).
The scattering of `V` into the cumulants shows up at second order in `t`, not first. Any test of
"does the closure build correlation from a determinant" must be run at finite `t`.

**4. The existing closures have the RIGHT AMOUNT of spectral motion and the WRONG PLACE for it.**
`lambda_3 = 0` and Mueller-2405 reproduce the spectral fraction of `Phi` (H4 t=1: 0.449 and 0.336
against 0.447 exact). Their failure is at the kernel: at H3 pure, where `Mat_2` has a 9-dimensional
exact kernel, the exact `x_i-dot` on the kernel is `1e-17` and the closures give `2e-3` to `2e-2`
(rank-jumping, the positivity-violation rate). At H4 pure there is no exact kernel, but the closures
push the smallest eigenvalues at up to `max|x_i-dot / x_i| = 29` (Mueller, t=1) against `2.7` exact --
ten times faster toward the boundary than the true dynamics.

**5. The spectral motion is not monotone-entropic in these systems.** `dS/dt` of the occupation
spectrum (fermionic entropy at rank 1, `-sum x ln x` at rank 2) is positive at `t = 1` after both
quenches, pure and thermal, but negative at H4 `t = 2` and slightly negative for the kicked ground
states. Four electrons in eight orbitals recur; an entropy-gradient ("metric") model of the spectral
motion is a coarse-graining assumption, not an exact statement, and will be right on average and wrong
at recurrences -- exactly as at rank 1.

## What this says about route F

- The reversible half is already structural: `drift = i[K_2, X]`. A closure of the form
  `Xdot = i[K_2 + K_3(gamma_1, gamma_2), X]` with Hermitian `K_3` approximating the non-spectral part
  of `Phi` is positive, energy-conserving and real-poled by construction, and is the rank-2 analogue
  of TDHF. Its known defect is the frozen geminal spectrum, which is 13-35% of the flow far from
  equilibrium.
- The spectral half of `Phi` is the rank-2 analogue of the collision integral. Modelling it as a
  gradient flow (metric bracket, entropy of the geminal spectrum, Fermi-type blocking at `x_i = 0`)
  gives the metriplectic form: positivity by boundary repulsion, energy conserved by construction,
  monotone entropy, a thermal fixed point. It is the same trade the published rank-1 theory made, one
  rank up, and result 5 says it is a trade, not a derivation.
- What is NOT yet known: whether the non-spectral part of `Phi` admits a `K_3` of bounded norm
  (near-degenerate `x_i` make `K_ij = Phi_ij / (i(x_j - x_i))` large), whether the contraction to
  `gamma_1` and the rank-2 trace identity are preserved by `i[K_2 + K_3, X]` (they are conserved by
  the exact flow; a commutator on pair space need not respect them), and whether a Jacobi-consistent
  `K_3(gamma_1, gamma_2)` exists at all. Those are the next three measurements, in that order.

---

## Episode 2 (2026-08-12): bounded `K_3`, and what the commutator closure costs the 1-RDM

Tool `tools/commutator_k3.py`; output `results/commutator_k3.txt`. Eight states (episode-1 set minus
the `t = 0` rows, which are uninformative by time reversal).

**`K_2` in closed form.** `D = i[K_2, X]` with

    K_2[(p,v),(q,w)] = h[p,q] d_vw + h[v,w] d_pq - h[p,w] d_vq - h[v,q] d_pw + V[p,v,q,w]

verified at `1e-16` on all eight states (coefficient 1 on `V` in the `qdiss.eom.hamiltonian`
convention). `|K_2|_2 = 2.8` (H3), `4.1` (H4).

**(i) Boundedness of `K_3`.** Split `Phi = Phi_spec + Phi_off` (block-diagonal / off-block in the
eigenbasis of `X`, tolerance `1e-6`). `Phi_off = i[K_3, X]` with the unique off-block Hermitian
`K_3[ij] = Phi_off[ij] / (i(x_j - x_i))`, reproduced to `1e-16`.

| state | `|K_3|_2` | `|K_3|_F` | min gap used | `sigma=1e-2`: `|K_3|_2`, frac of `Phi_off` kept |
|---|---|---|---|---|
| H3 kicked ground | 0.33 | 0.66 | 2.2e-2 | 0.32, 0.997 |
| H3 quench t=1 | 0.19 | 0.34 | 1.5e-2 | 0.18, 0.999 |
| H3 quench t=2 | 0.58 | 1.15 | 3.7e-2 | 0.55, 0.995 |
| H4 kicked ground | 0.54 | 1.58 | 1.1e-5 | 0.53, 0.996 |
| H4 kicked thermal | 0.25 | 0.49 | 2.6e-3 | 0.25, 0.999 |
| H4 quench t=1 | 0.37 | 1.11 | 1.9e-6 | 0.29, 0.995 |
| **H4 quench t=2** | **7.05** | **17.8** | 2.1e-5 | 1.46, **0.858** |
| H4 quench thermal t=1 | 0.28 | 0.56 | 8.3e-4 | 0.28, 0.999 |

In seven of eight states `K_3` is an order of magnitude smaller than `K_2`. In one (H4 quench `t=2`)
a near-degenerate pair of geminal occupations with a non-vanishing `Phi_off` element makes
`|K_3|_2 = 7`, and the Tikhonov regularisation that tames it to `1.5` discards 14% of `Phi_off`.
**The bound is not uniform.** This is the Lipschitz screen reappearing in F: the exact non-spectral
collision flow is representable by a Hermitian generator, but not by a bounded one near degeneracies
of `Mat_2`. A closure that SPECIFIES `K_3(gamma_1, gamma_2)` is not subject to this (it is bounded by
construction); the point is that reproducing the exact off-block flow is not always possible with a
bounded one.

**(ii) Contraction to the 1-RDM flow.** With `c(T)[p,q] = sum_v T[p,v,q,v] / (N-1)`, the induced
`gamma_1-dot` of the exact flow is `c(D) + c(Phi)` and equals `fullci` at `1e-16` on every row. A
commutator closure keeps `c(D) + c(Phi_off)` and loses `c(Phi_spec)`:

| state | `|c(Phi_spec)| / |c(Phi)|` | `|c(Phi_spec)| / |gamma_1-dot|` |
|---|---|---|
| H3 kicked ground | 0.058 | 0.011 |
| H4 kicked ground | 0.003 | 0.002 |
| H4 kicked thermal | 0.003 | 0.001 |
| H3 quench t=1 | 0.47 | 0.29 |
| H3 quench t=2 | 0.38 | 0.22 |
| H4 quench t=1 | 0.43 | 0.38 |
| H4 quench t=2 | 0.14 | 0.08 |
| H4 quench thermal t=1 | 0.42 | 0.38 |

Near equilibrium the loss is 0.1-1% of the 1-RDM flow. Far from it, **22-38% of `gamma_1-dot` is
lost**, and the lost part is essentially pure natural-occupation motion (spectral fraction of
`c(Phi_spec)` w.r.t. `gamma_1` is 0.74-1.00 on every row). Since `gamma_1` is obtained by contraction,
the rank-2 trace identity `Tr_2 lambda_2 = -gamma_1 eta_1` holds identically; the cost appears not as a
broken identity but as a WRONG 1-RDM flow of the same relative size as the rank-2 error found in
episode 1. Note also that the induced flow `c(D) + c(Phi_off)` still moves the natural occupations
(spectral fraction 0.12-0.34 under the quench, against 0.45-0.58 exact): a pair-space commutator does
not freeze `gamma_1`'s spectrum, it just moves it by the wrong amount.

**Standing of route F after two episodes.**

- Reversible half: structural and explicit (`K_2`), with a bounded Hermitian `K_3` recovering the
  non-spectral collision flow at 7 of 8 states and requiring regularisation at the eighth.
- Cost of a pure commutator closure: 13-35% of the rank-2 flow and 22-38% of the induced rank-1 flow
  far from equilibrium; ~1% near it. Those two numbers are the size of the dissipator the metriplectic
  half must supply, and they are dominated by occupation (spectral) motion at both ranks -- which is
  exactly what a metric bracket with an occupation entropy acts on.
- Not yet done: the Jacobi question for a state-dependent `K_3`; and the actual construction of a
  metric bracket for the spectral part with the right fixed point. Both are now well-posed.

---

## Episode 3 (2026-08-12): energy, the Hamiltonian pair flow, and a contraction identity

Tool `tools/hamiltonian_pair_flow.py`; output `results/hamiltonian_pair_flow.txt`.
Symbols as above, plus: `K_h`, `K_V` the one- and two-body parts of `K_2` (`K_2 = K_h + K_V`);
`E = tr(h gamma_1) + (1/4) V.gamma_2` the physical energy, linear in `X`; `grad E` its (constant)
gradient on pair space; `c(T) = Tr_2 T / (N-1)`.

**Energy.** `D = i[K_2, X]` does NOT conserve `E`: `dE/dt|_D = -4e-3 .. -7e-2`, cancelled exactly by
`Phi`. So the drift conserves the spectrum of `Mat_2` and the quantity `tr(K_2 X)` (an energy-like
invariant with one- and two-body parts weighted `2(N-1) : 2` instead of `1 : 1/4`), but not `E`.

**`grad E` in closed form.** Measured: `grad E = K_h / (2(N-1)) + K_V / 2`, linearity check 1e-15.

**Contraction identity (1e-16 on all rows).**

    c( i[K_h, X] ) = i[h, gamma_1]
    c( Phi )       = (N-2) . c( i[K_V, X] )
    =>  gamma_1-dot = i[h, gamma_1] + (N-1) . c( i[K_V, X] )        exactly

The entire two-body part of the 1-RDM flow is the contraction of a pair-space commutator, with weight
`N-1` on `K_V` where the drift has weight 1. Consequences:

- The `u(pairs)` Hamiltonian flow `i[grad E, X]` induces exactly `gamma_1-dot / (2(N-1))`: the exact
  1-RDM flow uniformly slowed. (This is why its induced rank-1 error was exactly `1 - 1/(2(N-1))`.)
- `K* = K_h + (N-1) K_V = 2(N-1) grad E` gives a commutator flow with EXACT induced rank-1 EOM,
  exact energy conservation, and structural positivity, no parameters. The 22-38% 1-RDM loss of
  episode 2 was the wrong weight on `K_V`, not a property of commutator closures.
- The exact flow is `Xdot = i[K*, X] + R`, `R = Phi - i[(N-2) K_V, X]`, with `c(R) = 0` identically.
  Any modelling of `R` leaves the rank-1 EOM exact.

**But `K*` is a bad rank-2 flow.** Relative error of `i[K*, X]` against `Xdot_exact`: 1.2-4.3 (worse
than `D` at 0.56-1.29, and than `lambda_3 = 0` at 0.06-0.19). The remainder `R` is 1.2-4.3 times the
exact flow, nearly all of it non-spectral (spectral fraction of `R` 0.0004-0.14), so `R` is very
nearly a commutator `i[K_3, X]` with a LARGE `K_3` that must cancel most of the `(N-2) K_V`
overweight everywhere in pair space while contributing nothing to the contraction.

**Reading.** Two exact linear constraints on a commutator closure `Xdot = i[K, X]` are now known:
`c(i[K, X]) = gamma_1-dot` (contraction) and `tr(X [grad E, K]) = 0` (energy). `K*` satisfies both and
is the worst flow measured; `D` satisfies neither and is better; the best commutator approximation to
the exact flow (`D + Phi_off`) has error equal to the spectral fraction (0.1-35%). The closure content
is a `K_3(gamma_1, gamma_2)` that interpolates: it must supply the pair-space structure of `Phi` (which
`lambda_3 = 0` already approximates to 6-19%) while respecting the two constraints. That is episode 4.

---

## Episode 4 (2026-08-12): anatomy of `Phi`; the geminal mean-field flow in closed form; the ph channel

Tools `tools/phi_anatomy.py`, `tools/phi_sandwich_fit.py`, `tools/phi_identify.py`,
`tools/geminal_meanfield.py`; outputs in `results/` under the same names.

Symbols as above, plus: `X_d = Mat(γ_1∧γ_1)`, `L = X − X_d = Mat(λ_2)`, `N = γ_1⊕γ_1` (pair-space
occupation matrix, `N[(pv),(qw)] = γ_1[pq]δ_vw + γ_1[vw]δ_pq − γ_1[pw]δ_vq − γ_1[vq]δ_pw`),
`K_{W^T} = W^T⊕W^T` with `W[p,q] = Σ_rs V[p,r,q,s] γ_1[r,s]`. By MKNO,
`γ_3 = A_3(γ_1γ_1γ_1) + A_3(γ_1λ_2) + λ_3`, so `Phi = Phi_ggg + Phi_gl + Phi_l3` (`Phi` is linear in `γ_3`).

**4a. Anatomy (8 states).** `Phi_ggg` is nearly spectral-free on H4 (0.0006–0.07; on H3 0.11–0.17
because of its 9-dim exact kernel block). `Phi_gl` is heavily spectral (0.4–0.84) — **the
geminal-occupation motion lives in the `γ_1λ_2` part of the separable 3-RDM.** `|Phi_gl| ≈ ⅓|Phi_ggg|`,
`|Phi_l3| ≈ 0.07–0.2 |Phi|`. `Phi_sep` reproduces the exact contraction identity
`c(Phi) = (N−2)c(i[K_V,X])` to 1e-3; its two halves do not individually.

**4b. Hermitian polynomial `K` in pp pair space cannot carry `Phi`.** Common-coefficient fits of
`Phi_piece ≈ Σ_a c_a i[K_a, X]` over ten Hermitian structures (`K_h, K_V, W⊕W, C_h⊕C_h, C_a⊕C_a`,
one- and two-leg `γ_1`-dressed `V`, `λ_2`-dressed `V` in pp and ph form): residual 0.38 (`Phi_ggg`),
0.85 (`Phi_gl`), 0.44 (`Phi`). Closures `i[K_2 + K_fit, X]` have rank-2 error 0.12–0.54 and
contraction error 0.10–0.55 — worse than `λ_3 = 0` (0.06–0.19, 0.005–0.06) on every state.
Fitting non-Hermitian `S` in `i(S†X − XS)` over 15 monomials did no better (0.31–0.70).

**4c. Exact identification of `Phi_ggg` (residual 1e-14 on six states, verified directly at 1e-15 on
eight).** The class of pair-space products `A·B·C`, `B ∈ {K_V, K_{W^T}}`, `A, C ∈ {I, N, X_d}` contains
`Phi_ggg` exactly, and the identity is

    Phi[A_3(γ_1γ_1γ_1)]  =  i( S† X_d − X_d S ),        S = K_{W^T} − K_V N

(with `W`, not `W^T`, it fails at 61%: the repo's `γ_1[p,q] = <a+_p a_q>` convention makes the
transpose essential for complex `γ_1`, as in the rank-1 EOM `i[(h+W)^T, γ_1]`). Together with the
drift, the entire `λ_2 = λ_3 = 0` rank-2 flow is one XR sandwich on `X_d`:

    D + Phi_ggg  =  i[K_h + K_V, X]  +  i(S† X_d − X_d S),     S_tot = K_h + K_{W^T} + K_V (I − N)

`K_V(I − N)` is the Bethe–Goldstone / pp-ladder phase-space factor `(1 − n_p − n_v)V`, and the
anti-Hermitian part `Γ = i(S − S†) = i[N, K_V]` — non-commutativity of occupation with interaction —
is exactly the spectral (occupation-moving) part of the mean-field-level flow. At a determinant
(`X = X_d`, `γ_1` real) this reduces to a rational identity with coefficients ±½, ±¼.

**`Phi_gl` is NOT in the pp product class** (53% residual even with `L` as a side and `C, C^T` as
cores). Its Wick patterns contract `V` with `λ_2` on one ket and one bra index — the
**particle-hole channel** — which is a different index reordering (the G-matrix) and cannot be a
product of pp pair-space matrices. This is structural: the pp Lie algebra `u(C(M,2))` sees the ladder
channel and is blind to the ring/exchange channel.

**4d. The geminal mean-field flow as a closure.** Two variants: GMF-d (as derived, sandwich on `X_d`
plus drift on `X`) and GMF-x (sandwich `i(S_tot†X − XS_tot)` on `X` itself, XR-structural).

| state | D | GMF-d | GMF-x | λ_3=0 | GMF-x `tr Xdot` |
|---|---|---|---|---|---|
| H3 kicked ground | 0.56 | 0.15 | 1.19 | 0.055 | −1.4e-2 |
| H3 quench t=1 | 1.29 | 0.38 | 2.01 | 0.085 | +1.4e-1 |
| H4 kicked ground | 0.96 | 0.35 | 2.52 | 0.135 | −1.2e-3 |
| H4 kicked thermal | 0.57 | 0.21 | 0.54 | 0.062 | −1.7e-4 |
| H4 quench t=1 | 0.79 | 0.41 | 1.76 | 0.096 | +2.2e-1 |
| H4 quench thermal t=1 | 1.05 | 0.40 | 1.03 | 0.079 | +4.0e-2 |

GMF-d halves the drift's error but is 2–3× worse than `λ_3 = 0`: the `γ_1λ_2` (ph-channel) part is
essential, not a correction. GMF-x is wrong: moving the sandwich from `X_d` to `X` costs
`−i(S†L − LS)`, which is large, and it breaks trace conservation (`tr Xdot = −tr(ΓX) ≠ 0`). The XR
form with this `S_tot` is therefore not a closure.

**Standing of route F after episode 4.**

- What is now derived: the pp-channel content of the rank-2 flow through `γ_3 = A_3(γ_1γ_1γ_1)` is a
  sandwich `i(S†X_d − X_dS)` with `S = K_{W^T} − K_V N`; its Hermitian part is the effective pair
  Hamiltonian `K_h + K_{W^T} + K_V − ½(K_V N + N K_V)` and its dissipative part is `Γ = i[N, K_V]`.
- What is now excluded: a Hermitian generator in pp pair space alone (any polynomial `K`) as an
  accurate closure. Half of the `λ_3 = 0` flow is particle-hole and lies outside that algebra.
- What follows: the ph channel has its own Gram matrix (`G`) and its own XR/commutator structure.
  The exact flow is a sum of channel sandwiches; `Phi_gl` should be identified as a product class in
  G-matrix space the same way `Phi_ggg` was in pair space. A two-channel effective generator is the
  natural object, and the `Q`/`G` row of the ledger — never controlled before — is exactly where it
  would act. That is episode 5.

---

## Episode 5 (2026-08-12): the particle-hole channel, and the whole `λ_3 = 0` flow in closed algebraic form

Tools `tools/phi_identify_ph.py`, `tools/phi_gl_sparse.py`, `tools/closed_form_sep.py`; outputs in `results/`.

**Symbols.** Pair (pp) space as before: `X = Mat_2`, `N = γ_1⊕γ_1`, `K_h`, `K_V`, `K_{W^T}`, `L = Mat(λ_2)`,
`X_d = X − L`. One-body `C[p,q] = Σ_rtu V[r,t,p,u] λ_2[r,t,q,u]` (the rank-1 collision matrix; a partial
trace of `V·λ_2`, not a matrix product), `K_C = C⊕C`, `K_{γC}`, `K_{Cγ}`. **ph space** (`M²×M²`, pair indices
ordered (ket, bra)): `ph(T)[(p,q),(v,w)] = T[p,v,q,w]`, `V_d = ph(V)`, `L_d = ph(λ_2)`; `Sw` swaps the two
entries of a pair index; dressings `γ_1⊗I`, `γ_1^T⊗I`, `I⊗γ_1`, `I⊗γ_1^T`. `back(A)` maps an ph matrix to a
4-index tensor, antisymmetrises over `(p,v)` and `(q,w)`, and returns `Mat_2/2`. The exchange pairing
`ph_x(T)[(p,w),(v,q)]` equals `−ph(T)` exactly by antisymmetry, so one pairing suffices.

**5a. Term-by-term classification of the 72 Wick patterns of `Phi[A_3(γ_1λ_2)]`** (coefficients common
across six states; per-state tests are vacuous because the class is underdetermined per state):

| pattern type (indices `V` shares with `λ_2`) | count | representation |
|---|---|---|
| one ket AND one bra | 32 | ph products `A·V_d·B·L_d·C†`, `A·L_d·B·V_d·C†` with dressings incl. `Sw` |
| both ket, or both bra | 8 | pp products `N·K_V·L`, `L·K_V·N` |
| one index; `V`'s other two closed on `γ_1` (W-type) | 16 | pp products `K_{W^T}·L`, `L·K_{W^T}` |
| three indices (`C`-type), `γ_1` in the other slot | 16 | `C⊗γ_1 + γ_1⊗C = ½{N, K_C} − K_{γC} − K_{Cγ}` |

With all four kinds the common residual for `Phi_gl` is **4e-14**; without the `C`-type cores 13%,
without the `Sw` dressings 60%, pp-only 56%, ph-only 33%.

**5b. Sparse identity** (IRLS/L1 selection on six states, re-fit, **verified on two held-out states at
8e-15**), coefficients ±2 and ±½:

    Phi[A_3(γ_1λ_2)] = i(S† L − L S) + i(Q − Q†) + Phi_ph
        S      = K_{W^T} − K_V N                                    (same S as in Phi_ggg)
        Q      = ½(N K_C − K_{γC}) + ½(K_C N − K_{Cγ})              (= C⊗γ_1 + γ_1⊗C, antisymmetrised)
        Phi_ph = 2i[(P_1 + P_2 − P_3 − P_4) − h.c.]  (after back-mapping)
        P_1 = Sw (γ_1^T⊗I) V_d L_d,   P_4 = Sw (I⊗γ_1) V_d L_d,
        P_2 = L_d V_d (Sw (γ_1⊗I))†,  P_3 = L_d V_d (Sw (I⊗γ_1^T))†

so the ph generator is `(γ_1^T⊗I − I⊗γ_1)·V_d`: the **occupation difference `(n_p − n_q)` times `V`**,
the RPA/GW structure — as `K_V(I − N)`, i.e. `(1 − n_p − n_v)V`, was the ladder structure in pp.

**5c. The whole `λ_3 = 0` rank-2 flow, closed form, verified directly at 1e-15 on eight states:**

    D + Phi[separable_3]  =  i[K_h + K_V, X]  +  i(S† X − X S)  +  i(Q − Q†)  +  Phi_ph

Note the pp sandwich now acts on the full `X` (Phi_ggg supplied the `X_d` part, Phi_gl the `L` part), which
resolves episode 4d: GMF-x failed because it lacked `i(Q − Q†) + Phi_ph`, not because the sandwich was on
the wrong object. The pieces are individually LARGER than the total and cancel strongly (H4 quench t=1:
sandwich 1.44, ph 1.36, total 0.75): the ladder and ring channels oppose each other, so neither channel
alone is a closure.

**Reading for route F.** The `λ_3 = 0` flow is a sum of three Lie-algebraic pieces plus one non-product term:
- pp channel: `i(S_pp† X − X S_pp)`, `S_pp = K_h + K_{W^T} + K_V(I − N)`, dissipative part `Γ_pp = i[N, K_V]`;
- ph channel: a two-sided action on `L_d = ph(λ_2)` with generator `(n_p − n_q)V_d`;
- `i(Q − Q†)`: the rank-1 collision matrix `C` tensored with `γ_1` — this is the term that makes the
  contraction identity `c(Phi) = (N−2)c(i[K_V,X])` work, and it is neither channel's product.
Every piece is polynomial in `(γ_1, γ_2, V)`, so every piece passes the Lipschitz and extensivity screens by
construction. The `λ_3` remainder is 5–20% of `|Phi|` and is what no separable structure carries.

**What follows (episode 6).** With the flow in channel form, "choose K for D, Q and G" becomes concrete:
the pp sandwich has the XR form on `X` (D-positivity structural if used alone); the ph piece acts on `ph(λ_2)`
and hence on the G-matrix `G = (γ_1⊗I)Sw − ph(γ_2)`; the question is whether each channel's action can be
written as a sandwich on ITS OWN Gram matrix (`X` for D, `G` for G, and the hole analogue for Q) with the
three sandwiches consistent on the shared `γ_2`. If yes, D, Q and G positivity are all structural for the
`λ_3 = 0` flow, and the `λ_3` remainder is the only thing left to model.
