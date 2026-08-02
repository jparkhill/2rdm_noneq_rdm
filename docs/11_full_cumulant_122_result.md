# 11 — The full-cumulant 1-2-2 triple commutator

Completed 2026-08-02. Derivation: `tools/mkno_stage.wl` (six staged Wolfram
products, checkpointed to `results/stage_*.m`). Evaluation and reduction:
`tools/term_audit.py` / `tools/orbit_reduce.py` (pure Python, no kernel).
Orbit-reduced answer: `results/commutator_122_orbits.txt`.

## The object

    [[ {a^p_q} , {a^{rt}_{su}} ] , {a^{vx}_{wy}} ]

built as `rv1v2 + v2v1r - v1rv2 - v2rv1`, with **every cumulant contraction
retained** — no `NoCumulants`.

## Result

| | terms | carrying λ |
|---|---|---|
| four triple products, each | 8,241 | 6,026 |
| joined commutator, raw | 32,964 | 24,104 |
| **after signed cancellation** | **1,312** | **560** |
| same, with `NoCumulants` applied | 752 | 0 |
| **DELTA — deleted by the old derivations** | **560** | |
| **orbit-reduced (independent content)** | **130** | **60** |
| same, `NoCumulants` | 70 | 0 |

λ census of the 1,312: **752 with no cumulant, 528 carrying λ₂, 32 carrying λ₃.**

**43% of the equation of motion was being deleted.** In reduced form: 60 of 130
independent terms.

Deleting λ before or after cancellation gives the identical 752 — λ-terms only
cancel against λ-terms — so the delta is unambiguous.

## Cancellation is where the real reduction lives

Within a *single* product, like-term merging gains 0% (docs/10). In the
*joined commutator* it is decisive: 32,964 → 12,599 distinct → 1,312 surviving.
**11,287 monomials cancel exactly.** That is why the earlier plan to simplify
the intermediates was worthless and this step is not.

Literal-identity cancellation and antisymmetry-aware canonicalization both give
1,312, i.e. the engine's index conventions are already consistent — the
cancellation is genuine cross-product cancellation, not bookkeeping.

## Symmetry: clean, once deltas are unpacked

With packed δ₂ factors expanded back to δ₁·δ₁, the final commutator is exactly
antisymmetric under all four external transpositions — `r<->t`, `s<->u`,
`v<->x`, `w<->y` — with **zero violations**. The 32 violations reported per
product in docs/10 were entirely a packing artefact and do not survive to here.

Orbit reduction under that 16-element group: 1,312 → **130 representatives**,
round trip verified exact (zero missing, zero extra, zero wrong coefficients).

130 terms is small enough to check by hand against Margócsy Eq. (2.1.53).

## Caveat on provenance

`AllStepsOfProductRule1PV` and `OneParticleVaccuum1/2` remain undefined
(docs/09), so this used the general `AllStepsOfProductRule` throughout — the
same substitution the notebook itself makes for the 2-1-1 chain at cells
1630-1638. This derivation is therefore **not** term-comparable to the stored
96-term cell-1624 list, which was built with the missing functions and carries
no cumulants at all.
