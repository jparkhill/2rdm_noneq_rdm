# 09 — Source material for the three open issues

Written 2026-08-02 for JP, who offered to help derive. Everything below is a
pointer to a specific line, cell, or equation — no summaries.

---

## Issue A — `NoCumulants` deletes every λ term

**The rule itself.** `MKNO.nb` cell 145, exported verbatim to
`docs/MKNO_inputs.wl` (search for `NoCumulants =`):

```mathematica
NoCumulants = {any1___,
               {{"Sign", a_}, any5___, {{λ, r_}, p___}, any8___},
               any4___} -> {any1, any4};
```

Read it literally: match any term containing a factor whose head is `λ` of any
rank `r_`, and rewrite the term list **without that term**. Not "set λ to zero
and keep the structure" — the whole monomial is dropped.

**Where it is applied** (cell numbers are positions in `docs/MKNO_inputs.wl`):

| cell | context |
|---|---|
| 153 | `tripcomm = AssociativeSimplify[All111Terms //. {NoCumulants, NormalVaccuumRule}]` |
| 230 | `SimplifiedL1Terms = FullSimplify[RecognizableForm[L1Terms //. {NoCumulants, NormalVaccuumRule}]]` |
| 244 | `com12 = v1r1 //. {QSpaceProjector, NoCumulants}` |
| 248 | inside `InvertAwayVacNormalOperatorsNoCum`: `//. {NormalVaccuumRule, NoCumulants}` (2 sites) |
| 250 | `com21 = AllStepsOfProductRulePProjector[com1, v2] //. {NoCumulants, NormalVaccuumRule}` |
| 516 | `com12 = v1r1 //. {NoCumulants}` |
| 519 | same as 248, second copy (2 sites) |
| 521 | same as 250, second copy |
| 1584 | `RecognizableForm[tpp //. {NoCumulants, NormalVaccuumRule}]` — the 1-2-2 product |

**Measured effect**, under the live kernel:

```
AllStepsOfProductRule[1-body, 2-body]  ->  30 terms, 15 carry λ
                    //. NoCumulants    ->  15 terms,  0 carry λ
```

So `AllStepsOfProductRule` is correct — it implements the full generalized Wick
theorem. The loss happens downstream, deliberately, by rewrite.

**What it reduces the algebra to.** One-body contractions only, built on γ.
That is the *quasi-normal-order theory* of Kong, Nooijen & Mukherjee,
J. Chem. Phys. **132**, 234107 (2010). Named and characterised in
`refs/papers/mkno/margocsy_thesis_elte2022.pdf`, **PDF page 27**: quasi-normal
ordering "does not satisfy ⟨{S}⟩ = 0 for an MR wave function."

Verified numerically: ⟨{a^pq_rs}⟩ = λ^pq_rs exactly, to 1.6e-16.
`tools/verify_gwt.py`, log at `results/verify_gwt.log`.

---

## Issue B — the missing functions (this is where I need you)

Three symbols are **used** but **defined nowhere** — not in `MKNO.nb`, not in
any repo under `refs/`:

| symbol | used in cells | occurrences in the raw .nb |
|---|---|---|
| `AllStepsOfProductRule1PV` | 1614, 1618, 1619, 1620 | 4, none followed by `:=` |
| `OneParticleVaccuum1` | 1615, 1621 | 2, none followed by `:=` |
| `OneParticleVaccuum2` | 1615, 1621 | 2, none followed by `:=` |

The chain that used them, verbatim:

```mathematica
(* 1611 *) SubstituteDeltas = False;
           fac1 = NRta1[{"p"},{"q"}];
           fac2 = NRta2[{"r","t"},{"s","u"}];
           fac3 = NRta2[{"v","x"},{"w","y"}];
           temp1 = AllStepsOfProductRule[fac1, fac2];
(* 1614 *) rv1v2 = AllStepsOfProductRule1PV[temp1, fac3];
(* 1615 *) FullSimplify[RecognizableForm[rv1v2 //. {OneParticleVaccuum1, OneParticleVaccuum2}]]
(* 1618 *) v1rv2 = AllStepsOfProductRule1PV[AllStepsOfProductRule[fac2,fac1], fac3];
(* 1619 *) v2v1r = AllStepsOfProductRule1PV[fac3, AllStepsOfProductRule[fac2,fac1]];
(* 1620 *) v2rv1 = AllStepsOfProductRule1PV[fac3, AllStepsOfProductRule[fac1,fac2]];
(* 1621 *) AllTerms = Join[rv1v2, v2v1r, MultiplyScalar[v1rv2,-1], MultiplyScalar[v2rv1,-1]];
           CorrelatedForm = Expand[Simplify[RecognizableForm[
                              AllTerms //. {OneParticleVaccuum1, OneParticleVaccuum2}]]]
(* 1622 *) sAllTerms = AssociativeSimplify[AllTerms]
(* 1623 *) <stored output: the 96-term list, γ and δ factors only>
(* 1624 *) Length[AssociativeSimplify[AllTerms]]   (* reported 96 *)
```

**The stored output of 1623 contains only γ and δ factors — no λ, no η.** So
whatever `OneParticleVaccuum1/2` did, the result is a one-particle-vacuum
reduction with cumulants already gone.

**The two questions only you can answer:**

1. Was `…1PV` a *variant of the product rule* (different contraction set) or the
   ordinary rule followed by a projection? The naming suggests it contracts
   against a one-particle vacuum rather than the correlated γ vacuum.
2. Were `OneParticleVaccuum1/2` rewrite rules of the same shape as
   `NormalVaccuumRule` (cell 145) — i.e. pattern → replacement on the term list?
   If so I can reconstruct them from the 1623 output as a spec.

If neither is recoverable, we route around it: derive the same commutator with
plain `AllStepsOfProductRule` and cumulants kept. That is what is running now
(`tools/mkno_full_cumulant.wl`), and it is what the notebook itself does for the
2-1-1 case at cells 1630–1638, so it is not an unprecedented substitution.

---

## Issue C — what to derive *against*

These are the published product rules, in the same index convention, to check
term by term.

**Rank-1 × rank-1** — `margocsy_thesis_elte2022.pdf`, **PDF page 29**, Eq. (2.1.52):

    {a^p_q}{a^r_s} = {a^pr_qs} + η^r_q {a^p_s} − γ^p_s {a^r_q} + γ^p_s η^r_q + λ^pr_qs

Verified as an exact Fock-space operator identity, error 4.4e-16. Deleting the
λ term breaks it by 1.4e-1.

**Rank-2 × rank-1** — same page, Eq. (2.1.53). **19 terms, 8 carrying λ₂ and 1
carrying λ₃.** This is the strongest derivation-to-code regression test in the
open literature and is the single most useful thing to diff our engine against.

**The general rule and its one exception** — same section:

> product of two MRGNO operators is evaluated by considering all possible
> contractions between them, each producing a K-body cumulant (γ₁ = λ₁
> understood). The one exception is the pair-contraction between an
> annihilation operator of {S₁} and a creation operator of {S₂}, which yields η₁.

**Sign rule:** (−1)^(π_U + π_L − n_η), upper and lower index parities counted
separately, n_η = number of η₁ factors. Worth diffing against whatever
`AllStepsOfProductRule` does internally.

**Cumulant expansions** — same PDF, page 18 Eq. (2.1.24)–(2.1.26), inverted at
page 19 Eq. (2.1.29).

**Particle–hole conjugation** — page 19–20:

    η₂ = λ₂ + ½ η₁∧η₁
    η₃ = −λ₃ + λ₂∧η₁ + ⅙ η₁^∧3

general rule: even-rank cumulants are p–h invariant, odd-rank (k≥3) flip sign.
Both verified numerically.

**Free consistency tests for any cumulant implementation** — page 19,
Eq. (2.1.31)ff, which the thesis explicitly offers as implementation checks:

    Σ_r λ^{pr}_{qr} = δ^p_q (n_q² − n_q)          [natural-orbital basis]
    Σ_t λ^{prt}_{qst} = λ^{pr}_{qs} (n_p+n_r+n_q+n_s−2)

Note the RHS of the first: n² − n ≤ 0, so partial traces of λ₂ are negative
semidefinite with magnitude equal to the occupation-number variance.

**Secondary confirmations of the same rules**, if you want independent eyes:
`2205.01178v1.pdf` (Evangelista, Wick&d) PDF page 2, Eqs. (3) and (6);
`hanauer_thesis_mainz2013.pdf` §3.2.2 Eqs. (3.49)–(3.56);
`1602.05667v1.pdf` (Li & Evangelista) Eqs. (19)–(21);
`hergert_frontiers2020_tour.pdf` Eqs. (25)–(26) — nuclear convention, same signs.

---

## Why it matters for positivity (one paragraph, so the target is clear)

    D = ½ γ∧γ + λ₂        Q = ½ η∧η + λ₂

Same λ₂, same sign, in both. And ½γ∧γ ⪰ 0 automatically whenever 0 ⪯ γ ⪯ I
(proved by projecting γ⊗γ ⪰ 0 onto the antisymmetric subspace; 300 random
draws, no violations). **Therefore λ₂ is the only place D ⪰ 0 or Q ⪰ 0 can
fail.** A trajectory with correlation generates λ₂; a generator derived with
every λ₂ term deleted has no channel that responds to it. That is the mechanism
to confirm or kill, and confirming it needs the full-cumulant term list.
