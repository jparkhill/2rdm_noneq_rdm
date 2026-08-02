# 10 — Is the derived term list adequately factorized?

Measured 2026-08-02 on `results/stage_rv1v2.m` and `results/stage_v1rv2.m`
(the two finished triple products of the full-cumulant 1-2-2 derivation,
8241 terms each). Tool: `tools/term_audit.py` — pure Python, no Wolfram kernel.

Three different things get called "simplification". They are not the same and
only one of them is currently being done.

## 1. Like-term merging — nothing left to gain (0%)

`AssociativeSimplify` = `RecognizableForm` -> `Simplify` -> `Expand` ->
`ListForm` (cell 141). All it can do is collect literally identical monomials.

**Measured: 8241 stored terms -> 8241 distinct. Zero merges**, even after
canonicalizing the antisymmetry of every γ_k, λ_k, η_k and {a_k}.

Consequence: the "simplify the 30-term intermediates before multiplying"
fallback staged in `tools/mkno_stage.wl` (`r1v1s`/`v1r1s`/`*s`) would have
bought **nothing**. It is not the lever. Do not spend a kernel on it.

## 2. Orbit factorization — 8x available, currently unexploited

The three factors are

    r  = {a^p_q}            no internal symmetry
    v1 = {a^{rt}_{su}}      antisym r<->t , antisym s<->u
    v2 = {a^{vx}_{wy}}      antisym v<->x , antisym w<->y

so the product must be antisymmetric under each of those four transpositions —
a group of order 2^4 = 16. **Verified: it is** (modulo §3).

    stored terms                     8241
    independent orbits               1008   (rv1v2) / 1019 (v1rv2)
    redundancy factor                ~8.2x

So ~87% of the term list is bookkeeping: eight signed copies of the same
content. The independent object is ~1000 terms, not 8241.

This is exact, not a truncation. Store one representative per orbit and
antisymmetrize on output. It shrinks every downstream step — the combine, the
spectral attribution, and the p-h audit — by the same factor, and it makes
term-by-term comparison against published rules tractable by hand.

## 3. A real defect: the packed delta is not canonical

32 terms in **each** stage break the `v<->x` antisymmetry. Every one of them,
in both stages, has the same shape signature: it contains **both** a `D1` and a
`D2`.

    16 x  D1 * D2 * G1 * ta1
     8 x  D1 * D2 * G1 * G1
     4 x  D1 * D2 * L2
     4 x  D1 * D2 * ta2

Cause. Cell 145 contains a rewrite that merges *a pair* of δ1 factors into a
packed δ2:

    {{d,1},{p},{q}}, ..., {{d,1},{r},{s}}  ->  {{d,2},{p,r},{q,s}}

When a term carries **three or more** deltas the rule pairs up whichever two it
matches first, so the leftover δ1 is chosen by term ordering rather than by any
canonical criterion. Two mathematically identical terms can therefore end up as
(δ1 on v)(δ2 on x) in one place and (δ1 on x)(δ2 on v) in another, and they no
longer match. The image term simply is not in the list under that key — the
audit shows the partner slot empty (coefficient 0), not wrong-signed.

This looks like a representation defect, not a physics error: the count is
identical (32) in two independently built products and confined to one shape
family. But it silently breaks term matching, symmetry checks and any
derivation-to-code diff — including the `//61{{Sign,...}}` stamps that tie
`RhoDotSec` back to the notebook. It affects the **old** derivations too.

Fix: canonicalize deltas before packing (sort the δ1 list by index, then merge
left to right), or do not pack at all and carry δ1 factors.

### Convention note

`{{d,2},{p,r},{q,s}}` is a **plain product** δ^p_q δ^r_s, not an antisymmetrized
pair. It must be canonicalized by sorting the (upper,lower) **pairs**. Treating
it as independently antisymmetric in its two index groups wrongly merges 159
distinct terms in `stage_rv1v2` — verified both ways.
