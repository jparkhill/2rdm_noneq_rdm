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
The flow is split as `Xdot = drift + Phi[gamma_3]`, with `drift` the `gamma_2`-carrying terms of the
rank-2 MKNO term list (`gamma_3 := 0` in the list; a term-list operation, not a moment negation) and
`Phi` the collision part, linear in `gamma_3`. Ground truth `qdiss.eom.fullci`; the wick term list
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
