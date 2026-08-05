# Prior art: Mueller, Luoma and Schaefer (arXiv:2405.05093)

"Hierarchical approach to quantum many-body systems in non-Markovian environments", Kai Mueller,
Kimmo Luoma, Christian Schaefer. Submitted 8 May 2024, latest revision 17 July 2026. Found by
John after we derived the partial-trace condition, which gave the right search terms.

Referred to throughout the code and notes as **Mueller**.

## What they do

BBGKY combined with Hierarchical Equations of Motion, for open many-body systems coupled to
structured photonic and phononic baths. They propagate two-body reduced density matrices
`F_12^(n,m)` carrying HEOM indices `(n,m)` that track bath memory, closing the BBGKY chain with a
reconstruction of `F_123`. Applications: spin squeezing, superradiance, molecules in organic
crystals, Hubbard-type electronic systems.

Normalisation is ours: `Tr F_1 = N`, `Tr F_12 = N(N-1)`, `Tr F_123 = N(N-1)(N-2)`. No conversion
needed.

Their closure neglects non-Gaussian 3-body correlations, `Gamma_123 = 0`, plus 4-body. They are
explicit that this is not the usual cumulant expansion: "unlike the widely used cumulant
expansion, our approximation still allows for non-Gaussian single-particle states, since only
non-Gaussian correlations *between* particles are neglected."

Their Eq. 14, in the closed-system limit where the HEOM indices drop and `Tr(F_1) = N`, collapses
after collecting like terms to

    F~_123 = (N-2)/N [ F_12 F_3 + F_13 F_2 + F_1 F_23 ] - 2(N-1)(N-2)/N^2 F_1 F_2 F_3

Two structural facts follow, checked by hand and then numerically in
`tools/muller_vs_ours.py`:

- `Tr_3(F~_123) = (N-2) F_12` **exactly**, using `Tr F_3 = N` and `Tr_3 F_13 = (N-1) F_1`. So the
  reconstruction is contraction-consistent *before* antisymmetrisation. Their Eq. 44/45 exists to
  repair what antisymmetrisation breaks, not the reconstruction itself.
- The factor `(N-2)` makes it vanish identically at `N = 2`, where the exact `gamma_3` is zero. So
  **Mueller is exact for two electrons**, which `lambda_3 = 0` is not: we measured the latter wrong
  by 4.6e-2 there, and traced that to `lambda_3` genuinely not vanishing at N=2.

They then antisymmetrise, `F^AS = Lambda- F~ Lambda- / norm` with
`Lambda-_123 = (1-P_12)(1-P_13-P_23)`, which is six times the full antisymmetriser; and restore
contraction consistency by Eq. 44, `F_12^Delta = F_12 - Tr_3(F_123^R)/(N-2)`, and Eq. 45,
`F_123^AS = F_123^R + (N-2) perp-Delta F_123`.

Positivity is handled by a numerical patch: a dynamical purification scheme keeping the smallest
eigenvalue of `F_12` above `-eps_p` with `eps_p = 1e-3`. They are candid that "not all
antisymmetric three-particle density matrices can be obtained as a reduced state of an
antisymmetric N-particle density matrix... In practice, these conditions are hard to check."

## Where this overlaps our work

The overlap is real and on the point that matters. Their statement that "enforcing contraction
consistency of the reconstructed 3-RDM `Tr_3(F_123) = F_12` and implementing a purification
scheme can significantly increase the stability" is the same finding as our Table 2 in
`docs/CLOSURE-CONDITIONS.md`: of all the physical conditions one might impose, contraction
consistency and positivity are the only two that bite, because hermiticity, antisymmetry,
particle-number conservation at both ranks and energy conservation are automatic.

## Where the two differ

| | Mueller | ours |
|---|---|---|
| trace consistency | numerical projection after the fact, Eq. 44/45 | closed-form partial-trace inverse `B[I,.]/(M-4)`, exact target derived |
| positivity | purification patch, floor `min eig F_12 >= -1e-3` | `CP at a point == Nagumo == X_KK >= 0`; for `N <= 3` this degenerates to the LINEAR equality `X_KK = 0`, because `rank gamma_2 <= C(M,N-2)` is conserved |
| rank | fixed at 3 | order-general single-trace rule `sum lambda_k = anti(E, lambda_{k-1}) + lambda_{k/2} * lambda_{k/2}` verified at k = 2, 3, 4 |
| exact at N=2 | yes, by the `(N-2)` prefactor | `gamma_3 = 0` yes; `lambda_3 = 0` no |
| open systems | HEOM, non-Markovian baths, detailed balance | none yet |
| spectra | not addressed | pole comparisons against TDHF and FCI |

## Deviations in our implementation of Mueller

Recorded in `qdiss/muller.py` rather than hidden. Their antisymmetriser is `Lambda-` on both
sides with an unstated normalisation; we apply the full antisymmetric projector and fix the scale
by requiring `Tr gamma_3 = N(N-1)(N-2)`, which is unambiguous. Their perpendicular projection
Eq. 43, used to build the Eq. 45 correction, was not present in the text we obtained, so the
correction is built with our closed-form partial-trace inverse, which produces a tensor with any
prescribed partial trace. Same purpose, different construction.

Also: this reading came from a summarising fetch of the arXiv HTML, not the LaTeX source. The
prose and the equations quoted above are reliable; anything depending on the precise form of
Eq. 43 is not, and should be checked against the source before being built upon.

## Open problems after reading it

1. **A stationary state.** Neither has one. Their accuracy "only begins to deviate slightly at
   long times when explicit 3-body correlations accumulate"; ours gives `|lambda_2-dot| = 1.37e-02`
   at the exact ground state. The same failure in two languages, and the condition John
   identified as the most defensible.
2. **Structural positivity versus patching.** They floor an eigenvalue; we can characterise the
   condition exactly for `N <= 3` and even impose it linearly, and the result was violently
   unstable. Neither route yet gives positivity together with usable dynamics.
3. **Why contraction consistency stabilises them and not us.** Their tested regime is
   `U/J <~ 1`, weakly correlated. Whether the condition is genuinely stabilising or merely
   sufficient at weak coupling is open, and testable by pushing their scheme to larger `U/J`.
4. **Reconstruction accuracy in the correlated regime.** They are bounded to `U/J <~ 1`; our
   reconstructions are far smaller in magnitude than the true cumulant. Nobody has one that
   works when correlation is strong.
5. **The hierarchy upward.** Our even-rank chain rule is a conjecture from three ranks, with
   sharp predictions at k=5 (no chain) and k=6 (`lambda_3 * lambda_3`). They stop at rank 3.
6. **The dissipation/ETH question is now directly testable in their framework.** They have
   non-Markovian baths and detailed balance and never ask whether the dissipation stabilises by
   driving the density dominant over the correlations. That is a well-posed experiment rather
   than a programme.
