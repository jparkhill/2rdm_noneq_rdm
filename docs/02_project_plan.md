# Project: a perturbative 2-RDM EOM that stays N-representable

## The generalisation to search for

At rank 1 the invariant-region argument (docs/01, §4.2) is:

> the N-representable set is $\{0\le\gamma\le 1\}$; on its boundary the flux points inward,
> because every gain term carries a particle factor $\gamma$ and every loss term a hole
> factor $(1-\gamma)$ — and those factors are *the MK contractions themselves*.

The rank-2 statement we want is the exact analogue. Let $\mathcal{C}$ be the
(approximate) N-representable cone — the DQG conditions:

$$^2D\succeq0,\qquad {}^2Q\succeq0,\qquad {}^2G\succeq0$$

(optionally $T_1,T_2$). **Boundary-flux criterion:** for every $^2D\in\partial\mathcal{C}$
with null vector $|v\rangle$ (i.e. $\langle v|M[^2D]|v\rangle=0$ for $M\in\{D,Q,G\}$),

$$\boxed{\;\langle v|\,\dot M\big[{}^2D\big]\,|v\rangle\;\ge\;0\;}$$

Given a closure (e.g. $\lambda_3=0$, or a parameterised family of reconstruction
functionals), $\dot M$ is an explicit polynomial in $\gamma$ and $\lambda_2$, so this is a
**decidable-in-principle polynomial positivity statement** — the right shape for
(a) Positivstellensatz / SOS certificates numerically, and (b) a Lean proof of the
general theorem once a candidate is found.

## Why this is a good agentic search problem

- **Cheap refutation.** Any candidate EOM can be falsified in minutes by propagating a
  random initial 2-RDM and watching the smallest DQG eigenvalue go negative. Verification
  is objective and fast — the loop can run unattended over many hypotheses.
- **Large but structured hypothesis space.** Free knobs: the cumulant reconstruction
  functional for $\lambda_3$; where the secular approximation is applied; whether the
  dissipator is written in a manifestly Lindblad (completely positive) form on the
  two-particle space before contraction; whether the rates are symmetrised so KMS holds
  term-by-term; contraction-consistency constraints tying $^2D$ to $\gamma$.
- **Strong prior from rank 1.** Manifest CP-form + term-by-term KMS is what made rank 1
  work. The rank-1 dissipator can be written with jump operators; if the rank-2 closure
  can be arranged so the *closed* equation is still the contraction of a genuine Lindblad
  flow on 2-particle space, positivity is structural rather than accidental.

## Pipeline

1. **Symbolic derivation.** Re-implement the MK generalized-normal-ordering engine
   (`MKNO.nb`) in Python/SymPy so derivations are scriptable, diffable, and machine-searchable.
   Validate by reproducing the notebook's rank-1 result exactly.
2. **Numerical falsifier.** Small-basis propagator + DQG eigenvalue monitor + Lyapunov check.
3. **Search.** Enumerate closures; auto-derive, auto-falsify, keep survivors.
4. **Lean.** Formalise the boundary-flux theorem for survivors. Rank-1 case first as a
   warm-up (it is short and true), then the rank-2 candidate.

## Status
- [x] Rank-1 EOM re-derived and verified (FD to $10^{-10}$, $0\le\gamma\le1$ invariant)
- [x] Confirmed against published C++ (`TCL.h` builds its dissipator from `HoleRho = I - Rho`)
- [ ] Symbolic MK engine in Python
- [ ] 2-RDM falsifier
- [ ] Lean formalisation
