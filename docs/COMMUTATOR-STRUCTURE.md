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
