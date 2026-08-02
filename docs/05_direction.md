# Direction, as of 2026-08-02

Agreed with JP this session. Supersedes the "search over closures" framing in
`02_project_plan.md` §"Why this is a good agentic search problem" — there is no
menu of closures to enumerate; the term list is *determined* by the algebra once
the cumulant truncation, the normal ordering, and the reference are fixed. The
work is derive-then-analyse, not sample-then-test.

## The strategy (JP's, retained)

Rank 1 worked like this: build an EOM that is positive and thermalizing, start it
anywhere, and let it **relax** to the physical state. The ground state was not
imposed; it was the attractor. The Fermi–Dirac fixed point and the invariance of
0 ≤ n ≤ 1 came from the same MK hole factors.

**We keep that hope at rank 2.** Do not solve for a stationary ²D and then ask
whether it is stable. Build an EOM whose flow cannot leave the N-representable
cone, and let relaxation find the ground state, exactly as before. Ground-state
*stationarity* is a later question (and a different one — see
`03_literature_and_positivity.md` §1, P1 vs P2). Positivity first.

## What "full analysis of the EOM" means concretely

For each derived EOM, in order, cheapest first. Each step is a gate: fail and we
stop and localise, rather than propagating and inspecting the wreckage.

1. **Derivation integrity.** p-h covariance of the term list (γ ↔ η with phase
   −1), and term-by-term agreement with the Mathematica oracle. Already passes
   on the un-secularized 1-2-2 commutator (cell 1624, 96 terms, all 64 monomials
   match). Not yet run on the secular cells — that needs the extractor fix, and
   those are the ones that became C++.

2. **Boundary flux at HF — the decisive local test.** HF is an extreme point of
   the DQG cone whose kernel is the ov ⊕ vv sector (measured: 22 of 28 geminals
   for H4, 46 of 91 for H2O). So positivity at the first order the theory acts
   is exactly

       P_ker ²Ḋ P_ker ⪰ 0,   i.e. the ov/vv block of ²Ḋ is PSD.

   A finite Hermitian matrix. Linear in the term coefficients, so if it is
   infeasible even with the coefficients free, that is a **theorem** that no
   theory of that structural form is representable at HF — not a failed sample.

3. **Poles from the same object.** The linearisation L = D**F**(²D_HF) whose
   kernel block is tested in (2) has, as its spectrum, the excitation energies.
   Compare against ADC(2) (not CIS(D) — CIS(D) is the non-iterative first
   iteration of partitioned ADC(2), so it has no pole structure of its own).
   One matrix answers positivity and poles together.

4. **Vector-field test against FCI.** Not trajectories — trajectories confound
   error accumulation with model error. Take exact ²D(t) from FCI, feed *that*
   into the approximate EOM, compare ²Ḋ_approx with ²Ḋ_exact pointwise. Errors
   do not accumulate; the residual can be decomposed by tensor structure; and
   the exact λ₃ can be extracted along the way to test the closure in isolation.

5. **Global cone invariance.** Only if (2) passes. This is where the spectral /
   square-root parameterisation and, if needed, an SOS certificate belong.

6. **Relaxation.** Only if (5) holds. Propagate from a screened reference and
   check the attractor.

## Standing constraints

- **Electronic Hamiltonians only** (PySCF, or BAGEL if PySCF cannot reach it).
  No Lipkin, no model Hamiltonians. Default systems: H4/STO-3G linear and
  square, H2O/STO-3G, short H-chains — all with instant FCI references.
- **Screen the reference before blaming the generator.** The MP2 2-RDM is badly
  non-representable (H2O: 69 of 91 Q eigenvalues negative). `Mp2NOGuess` was an
  available initial condition in EE2, so some historical failures may be
  initial-condition failures. `scripts/check_reference_representability.py`.
- **Diagonalise as a diagnostic, propagate the matrix.** Near the boundary the
  null space is degenerate; individual geminal occupations are not trackable
  through a degeneracy, and the 1/(λ_h − λ_g) rotations blow up. The condition
  is on the kernel *block*.
- **No fitted coefficients.** Free c_i have exactly one legitimate use: proving
  infeasibility of a whole structural form.

## Immediate queue

- [ ] Fix `tools/extract_terms.py` to survive rank-2 tensors (currently drops
      rank-1 factors after a V, so cells 417/426/430/448 come back as bare V).
- [ ] p-h test on the secular cells. Current suspect for the symmetry break.
- [ ] Activate Wolfram Engine; use `MKNO.nb` as the derivation oracle.
- [ ] Build the HF-linearisation object (gate 2 + gate 3).
- [ ] Resolve the `RhoDot` (`+=`) vs `RhoDotGood` (`-=`) sign on `Eta*I4*Eta` —
      the ov/vv-block term whose sign gate (2) says is decisive.
