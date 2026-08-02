# The oracle is live, and it says the engine was run with the cumulants off

Companion to `docs/07_mkno_primary_sources.md` (literature) — this note records
what the *running* Wolfram Engine says about `MKNO.nb` itself.

## 1. Setup

Wolfram Engine 14.1 is installed on the agent computer and licensed through an
on-demand entitlement (see `tools/wl`; billed against Wolfram cloud credits,
~10/kernel-hour, entitlement expires after 7 days).

- `docs/MKNO_inputs.{wl,m}` — all 346 Input cells of `MKNO.nb` exported as plain
  Wolfram source in notebook order, with cell indices.
- `tools/mkno_replay.wl` — replays cells in order under a per-cell time limit.
  260 ok / 40 error / 1 timeout before the run died at cell 1549 on a
  `RecursionLimit` in unrelated scratch work.
- `tools/mkno_oracle.wl` — loads only the definition cells and re-runs a
  targeted derivation.

## 2. The engine implements the FULL generalized Wick theorem

`AllStepsOfProductRule` genuinely produces cumulant contractions. Measured, for
the product of a 1-body and a 2-body normal-ordered operator:

```
AllStepsOfProductRule[1-body, 2-body]  ->  30 terms, 15 of which carry lambda
after applying the notebook's NoCumulants rule -> 15 terms, 0 with lambda
```

**Exactly half the terms are cumulant terms.** This matches the textbook GWT
(`docs/07`, verified numerically in `tools/verify_gwt.py` to 4e-16):

$$\{a^p_q\}\{a^r_s\} = \{a^{pr}_{qs}\} + \eta^r_q\{a^p_s\} - \gamma^p_s\{a^r_q\}
  + \gamma^p_s\eta^r_q + \boxed{\lambda^{pr}_{qs}}$$

So the engine is not naive. It knows about $\lambda$.

## 3. …but the derivations were run with the cumulants deleted

Cell 145 defines a rewrite rule whose only job is to delete any term containing
a $\lambda$ factor:

```
NoCumulants = {any1___, {{"Sign", a_}, any5___, {{λ, r_}, p___}, any8___}, any4___}
              -> {any1, any4};
```

with a companion `InvertAwayVacNormalOperatorsNoCum`. `NoCumulants` is applied in
cells 153, 204, 230, 244, 248, 250, 516, 519, 521, 1584.

Consequence, from `docs/07` §(a): the defining property of generalized normal
ordering, $\langle\{\ldots\}\rangle = 0$, then fails — and it fails by *exactly*
$\lambda_2$ (verified numerically, deviation 1.6e-16). What the notebook
implements after `NoCumulants` is the **quasi-normal-order theory** of Kong,
Nooijen & Mukherjee (JCP 132, 234107, 2010): one-body contractions only, built
on $\gamma$ alone. That is a real, self-consistent algebra — but it is *exact*
only for quasi-free (fermionic-Gaussian) references, i.e. exactly when
$\lambda_k \equiv 0$.

**This explains why rank 1 worked and rank 2 did not, without appealing to any
symmetry breaking.** A mean-field thermal state has a non-idempotent $\gamma$ and
all cumulants exactly zero (verified). The rank-1 Fermi–Dirac programme therefore
sat *inside* the engine's domain of exactness. A correlated 2-RDM does not.

## 4. Why that is the positivity mechanism

From `docs/07`, verified numerically here:

- $D = \tfrac12\,\gamma\wedge\gamma + \lambda_2$ and
  $Q = \tfrac12\,\eta\wedge\eta + \lambda_2$ — the *same* $\lambda_2$, same sign
  (even-rank cumulants are particle-hole invariant).
- $\tfrac12\,\gamma\wedge\gamma \succeq 0$ automatically whenever
  $0 \preceq \gamma \preceq I$ (300 random draws, no failures).

Therefore **$\lambda_2$ is the only place $D \succeq 0$ or $Q \succeq 0$ can
fail.** This is also the rigorous form of the particle-hole conjecture: D and Q
are not two conditions that happen to travel together, they are two constraints
on one object.

Now put those together with §3. The trajectory generates $\lambda_2 \neq 0$
(it must — that is what correlation *is*), while the generator was derived with
every $\lambda_2$ term deleted. So the EOM **produces** the one quantity that can
violate positivity and **contains no term that responds to it**. There is no
feedback channel by which a growing $\lambda_2$ can be checked. That is a
mechanism, and it is testable: restore the $\lambda_2$ terms and re-run the
boundary-flux test.

Status: hypothesis, consistent with everything measured so far. Not yet proven
to be *the* cause of the historical non-positivity.

## 5. A hard blocker on reproducing cell 1624

The 96-term 1-2-2 triple commutator (cells 1611–1624) is built with

- `AllStepsOfProductRule1PV`  ("1PV" = one-particle vacuum)
- `OneParticleVaccuum1`, `OneParticleVaccuum2`

**None of these three is defined anywhere in `MKNO.nb`, or in any other repo in
`refs/`.** They appear only at their four use sites (cells 1614, 1615, 1618,
1619, 1620, 1621). They were defined in a live kernel session that was never
saved.

So cell 1624's output **cannot be regenerated from the notebook alone**. The
96-term list we have is a stored output with no reproducible provenance. Its
name says what it did — restrict contractions to a one-particle vacuum, i.e.
drop the cumulants — which is consistent with the list containing only
$\gamma/\delta$ factors, and with its exact particle-hole antisymmetry
(`qdiss/symbolic/phtest.py`, 64/64, phase $-1$).

Options: (a) reconstruct `AllStepsOfProductRule1PV` from its behaviour, using
the stored 96-term output as the specification; (b) skip it and derive the
rank-2 EOM with the *full* GWT, keeping $\lambda_2$, which is where §4 says we
need to go anyway. (b) is the better use of the engine.
