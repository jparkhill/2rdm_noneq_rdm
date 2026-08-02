# 12 — The full-cumulant derivation is invalid at rank 5 (engine defect)

Found 2026-08-02 while trying to certify the 1-2-2 result numerically before
using it. Recorded here because it invalidates part of docs/11.

## What was checked, and how

`tools/gno.py` builds Mukherjee–Kutzelnigg generalized-normal-ordered operators
`{a^P_Q}` and cumulants `λ_k` of **any rank** in an explicit Fock space. It is
self-validated, not assumed:

| check | result |
|---|---|
| `⟨{a_k}⟩ = 0`, k = 1,2,3 (all index tuples) | ≤ 1.1e-16 |
| `⟨{a_k}⟩ = 0`, k = 4,5 | ≤ 7.8e-18 |
| `λ_k = 0` at a single determinant, k = 2..5 | exactly 0 |
| `λ_2`, `λ_3` vs the independent formulas in `tools/verify_gwt.py` | ≤ 4.2e-17 |

`tools/verify_stage_operator.py` then tests a derived stage as an **operator
identity** — both sides applied to a random Fock vector — so the residual
`{a_k}` coefficients are tested too, not just the scalar part.

## Result

| stage | what it is | verdict |
|---|---|---|
| `r1v1` = `{a^p_q}·{a^rt_su}` | rank 1 × rank 2 | **exact**, 8.9e-16 |
| `rv1v2` = `r1v1·{a^vx_wy}` | up to rank 3 × rank 2 | **WRONG**, error 1.06 against a typical magnitude of 0.67 |

The closed 368-term EOM likewise disagrees with brute-force evaluation of
`⟨[[{a^p_q},{a^rt_su}],{a^vx_wy}]⟩` on ~13% of index tuples, by an amount the
size of the answer. Failures are confined to tuples with **repeated external
indices**; with all ten distinct both sides vanish, so that case proves nothing.

The 96-term cumulant-free part (= cell 1624) is **not** implicated: it matches
exactly, and at a single-determinant reference both sides vanish identically.

## Cause

`MKNO_inputs.wl` cell 62 defines the recursion that converts a tilde operator
into normal-ordered form:

    TildeSubstitution[ui_, li_] := If[Length[ui]==1, ta1[ui,li],
      If[Length[ui]==2, ta2[...], If[Length[ui]==3, ta3[...],
      If[Length[ui]==4, ta4[...],
        Print["Error, Please implement the corresponding recursive
               decomposition of tilde-a"]]]]]

and `NRTildeSubstitution` in the same cell has the identical cutoff. `Print`
returns `Null`, which is then fed to `DistributeOperatorProduct` — so the
failure is silent in the result and only visible as a message.

Rank 5 **is** implemented — `ta5`, `Na5`, `NRta5`, `Gd5` are all defined in cell
125, and cell 146 extends `InversionOfVacNormal` and `InversionOfGamma` to rank
5. But cell 146 extends them alongside `InvertAwayVacNormalOperatorsNoCum` and
`InvertAwayGammasNoCum`, i.e. **on the NoCumulants branch**. The full-cumulant
path through `AllTildesToNormals` → `TildeSubstitution` was never extended.

That explains the historical record exactly: the notebook's own derivations ran
through `NoCumulants`, where rank 5 works, so the gap was never exercised. A
1-2-2 product reaches rank 1+2+2 = 5, so the moment cumulants are retained the
derivation walks straight into it.

`AllTildesToNormals` also contains `If[TildeRank > 3, Print["Error! Please
recurse out the tilde operators of rank 4"]]` — a warning, not a stop.

## What survives from docs/11

- Cell 1624 reproduced exactly (96 terms, signs included) — stands.
- The channel classification (scalar / ta1 / ta2 / ta3) — stands as a statement
  about the term list, but the term list itself is not yet trustworthy.
- **The 1312 / 752 / 560 counts and the 130 orbit representatives do NOT stand**
  until the rank-5 case is implemented and the derivation re-run.
- The δ-packing finding and the orbit-reduction machinery are independent of
  this and stand.

## Fix

Extend `TildeSubstitution`, `ASubstitution` and `NRTildeSubstitution` to rank 5
using the already-defined `ta5` / `Na5` / `NRta5`, then re-derive. Rank 5 is
sufficient: a 1-2-2 product cannot exceed it.
