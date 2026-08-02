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

## Independent cross-check — two engines agree exactly

The Wolfram combine (`AssociativeSimplify` on the joined list) and the pure-Python
signed cancellation (`tools/term_audit.py`) were run independently on the same
raw joined commutator. They agree **term for term and coefficient for
coefficient**:

    Wolfram combine : 1312 terms, 560 carry lambda
    Python cancel   : 1312 terms, 560 carry lambda
    IDENTICAL       : True   (0 missing, 0 extra, 0 wrong coefficients)
    NoCumulants     : 752 vs 752, identical

Wolfram's own lambda census matches too: 752 with no cumulant, 528 with a single
lambda2, 32 with a single lambda3. Reproduce with `tools/xcheck_combine.py`.

So the 560-term delta does not rest on my parser: two independent
implementations of the cancellation reach the same equation.


## Deriving *in* the reduced representation — measured, and it does not work as a wrapper

The after-the-fact orbit reduction above is exact. The open question was whether
the engine could **derive** in that representation, never building the images.

Tested directly (`tools/mkno_repproduct.wl`): reduce the `r x v1` intermediate
from 30 terms to its 13 orbit representatives, then multiply *those* by `v2`.

| | terms |
|---|---|
| input, full | 30 |
| input, representatives | 13  (2.3x smaller) |
| output, from full input | 8,241 |
| output, from representatives | 7,711  (**1.07x smaller**) |

**Two things go wrong.**

1. **The saving does not propagate.** A 2.3x smaller input gives a 6% smaller
   output. The redundancy in the *product* is dominated by `v2`'s own symmetry
   (`v<->x`, `w<->y`), which is still fully expanded, and by the contraction
   combinatorics — not by the input's `r<->t`/`s<->u` images. Reducing one
   factor buys almost nothing.

2. **The coefficients come out wrong.** Antisymmetrizing the representative
   product gives exactly the right term *set* — 0 missing, 0 extra out of 8,241
   — but 7,904 coefficients disagree, with ratios 1, 3/2, 2, 9/4, 5/2, 8/3, …
   Non-integer ratios are the signature: the engine merges terms whose orbits
   have different stabilizers, and the merged output no longer records which
   representative each term came from, which is exactly the provenance the
   orbit-stabilizer weight needs.

**Conclusion.** Feeding representatives to the existing product rule is not the
route. Doing this properly means the contraction engine itself operating on
canonical antisymmetric tensors — carrying index symmetry inside the tensor and
never expanding it into separate terms. That is what Wick&d
(`refs/papers/mkno/2205.01178v1.pdf`) does, and it is an engine rewrite, not a
wrapper.

Nothing above depends on it: the 130-term reduced answer is obtained
after the fact and is verified exact.
