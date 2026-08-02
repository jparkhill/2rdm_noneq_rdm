# Positivity from principle, diagnosed in the spectral representation

Supersedes the "find where the secular approximation broke the symmetry"
line of attack. JP's correction, recorded verbatim in substance:

> I was really pulling at strings on the secular thing. I was merely HOPING
> it would be positive in analogy with secular Redfield — I had no reason to
> believe it would be.

Agreed, and the analogy is worth dissecting because dissecting it tells us
what the rank-2 condition actually is.

---

## 1. Why secular Redfield is positive, and why that does not transfer

Secularization is not a positivity-preserving operation. It happens that at
rank 1, applying it produces a generator that is *already in GKSL form*:

    ṅ_p = Σ_q [ W_{p←q} n_q (1 − n_p) − W_{q←p} n_p (1 − n_q) ]

Positivity here is structural and has nothing to do with the word "secular".
Every **loss** channel out of p carries the factor n_p, and every **gain**
channel into p carries (1 − n_p). So on the boundary the offending flux is
multiplied by zero:

    n_p = 0 ⟹ ṅ_p = Σ_q W_{p←q} n_q ≥ 0
    n_p = 1 ⟹ ṅ_p = −Σ_q W_{q←p} (1 − n_q) ≤ 0

Secularization was the *route* to that form in the Redfield case; it is not
the *reason* the form is positive. The reason is: **the generator is a sum of
jump terms whose rates are gated by the occupation of the state they move
population out of.** Any derivation that lands on that structure is positive.
Any that doesn't, isn't — secular or otherwise.

## 2. The variable whose positivity we need is not the variable the blocking
   factors are built from

This is, I think, the actual structural reason the rank-1 success did not
transfer, and it only becomes visible in the spectral representation.

MK generalized normal ordering is taken with respect to **γ**. Its elementary
contractions are

    ⟨a†_p a_q⟩ = γ_qp ,   ⟨a_p a†_q⟩ = δ_pq − γ_pq = η_pq

so every blocking factor the algebra can manufacture is a **one-particle**
hole factor η. But the object whose positivity we need at rank 2 is ²D, and
its spectral decomposition is over **geminals**:

    ²D = Σ_g λ_g |g⟩⟨g| ,   |g⟩ a two-particle state (Coleman natural geminals)

The positivity variables are the λ_g. Nothing in the algebra guarantees that a
generator whose blocking factors are (1 − n_p) produces loss channels gated by
λ_g.

**At rank 1 these two objects coincided** — the blocking factor and the
positivity variable were the same numbers n_p. At rank 2 they are different
objects living in different spaces. That coincidence, not any deep property of
MK ordering, is what made rank 1 work, and it is exactly what is absent at
rank 2.

### 2.1 What the pair-hole factor would have to be

The rank-2 analogue of η = 1 − γ is the **two-particle hole matrix**

    η₂ = (I ∧ I) − ²D   ... which, up to the γ-linear terms, is ²Q.

So "the generator should carry pair-hole blocking factors" and "the generator
should be p-h covariant so that Q is preserved alongside D" are **the same
statement**. ²Q ⪰ 0 *is* the pair-hole occupation being non-negative. This
collapses two separate items on the list into one, and it says where to look:

> A positivity-preserving rank-2 EOM must have its loss channel for geminal g
> gated by λ_g, and its gain channel gated by the corresponding eigenvalue of
> ²Q. If the derived EOM does not have that structure, that is the defect —
> stated as a property of the generator, not as an observed crash.

This is a conjecture about the required form. It is checkable term by term.

## 3. The diagnosis procedure (JP's proposal, made concrete)

Given a derived ²Ḋ = F(γ, ²D) with no secular approximation and no bath:

1. **Spectral transform.** ²D = U Λ U†. Then exactly
   λ̇_g = ⟨g| F |g⟩ and U̇ = U X, X_{gh} = ⟨g|F|h⟩/(λ_h − λ_g).
   The eigenvalue equation carries **no energy denominators** — every
   denominator lives in the eigenvector rotation, which positivity ignores.
   This is why the spectral view isolates the culprit: nothing can hide behind
   a near-singular denominator.

2. **Per-term attribution.** F is a sum of named contractions. Evaluate
   ⟨g| F_i |g⟩ for each term i separately, on states with λ_g → 0. Any term
   with a strictly negative limit is a **named** culprit: "this contraction,
   in this channel," not "the propagation went negative at t = 3.7."

3. **Degeneracy.** Near ∂C the null space is large, so the correct object is
   the block P_ker F P_ker and the condition is that this matrix is PSD — not
   that its diagonal entries are non-negative. Per-geminal occupations are
   meaningless through a degeneracy. Propagate the matrix; diagonalise only to
   diagnose.

4. **Sort into GKSL-compatible and remainder.** Attempt to write
   F = −i[H_eff, ²D] + Σ_k ( L_k ²D L_k† − ½{L_k†L_k, ²D} ) + R
   on the two-particle space. Everything absorbed into the first two pieces is
   positive by construction. **R is the entire problem**, and it is now an
   explicit object rather than a symptom.

5. **Order-count the remainder.** A second-order theory is *entitled* to
   violate N-representability at third order. If ‖R‖ = O(λ³) the theory is
   honest and the fix is a resummation/reordering. If ‖R‖ = O(λ²) the closure
   is structurally wrong and no amount of reorganisation saves it.

Step 5 is the falsifier, and it needs no exact reference.

## 4. Consequences for the queue

- **Dropped:** the hunt for where secularization broke p-h symmetry. The
  finding stands as archaeology (`ef87177`: 32 of 48 monomials in
  `RhoDotSecFull` have no p-h partner at all) but it is not on the critical
  path, because there was never a principled reason to expect the secular form
  to be positive.
- **Kept and promoted:** derive the closed, non-secular EOM; spectrally
  decompose; attribute per term; attempt GKSL sorting; order-count the
  remainder.
- **Still first:** the derivation itself must be trustworthy, hence the
  Mathematica oracle and a Python engine validated against it.
- **Unchanged constraints:** electronic Hamiltonians only; screen the reference
  for DQG-representability before blaming the generator (MP2's 2-RDM is badly
  non-representable — see `docs/04`); no fitted coefficients except to prove
  infeasibility.
