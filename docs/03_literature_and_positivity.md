# Literature dossier, and three different things called "positivity"

Sources mirrored in `refs/papers/lit/` (PDF + extracted text, 16 papers).
Findings below were gathered by a research pass and then checked against the
primary text where a claim is load-bearing; claims still unverified against
primary text are marked **[unverified]**.

---

## 0. Correction to an earlier claim of mine

I previously wrote that Tohyama & Schuck (arXiv:1003.0246) show "the pure rank-2
closure fails the stability condition," and implied this bears on our
N-representability problem. **That was an over-reading.** Checked against the
primary text (§III.C, p.595 ff. of the extracted text): their stability
condition is

> E = ⟨0|e^{-iF̂} Ĥ e^{iF̂}|0⟩ with F̂ Hermitian one-body; the condition that
> |0⟩ "corresponds to a minimum in the energy surface" is positive
> semidefiniteness of the stability matrix of second variations.

That is a **variational** statement about a stationary point of the energy — the
same object as an RPA/Hessian stability analysis. It is not a statement about
whether ²D remains N-representable along a trajectory. JP flagged this; he is
correct. The two are logically independent (see §1).

Their result is still worth having, but it belongs to the *ground-state*
workstream, not the positivity one. Their own finding is also weaker than the
abstract suggests: for Lipkin, all eigenvalues associated with the *physical*
operators (J₊, J₋, J_z) are positive definite; the negative eigenvalues occur
"for unphysical states" — non-collective phonons — and are "at most 1% of the
maximum positive eigenvalue."

---

## 1. Three distinct notions of positivity

Much of the confusion in this literature (and in my own earlier summary) comes
from three unrelated conditions sharing a word.

| | object | type of statement | when it fails |
|---|---|---|---|
| **(P1) N-representability** | ²D, ²Q, ²G ⪰ 0 (DQG cone `C`) | *pointwise* constraint on the state, at each time | trajectory leaves `C` |
| **(P2) Ground-state stability** | stability matrix of ∂²E ⪰ 0 | *local* property of a stationary point | stationary point is a saddle, not a minimum |
| **(P3) Strength-function positivity** | spectral weights of the response ≥ 0 | property of the *linear response* spectrum | ERPA matrix non-Hermitian ⇒ complex/negative-weight poles |

**These are independent.** A flow can keep `C` forward-invariant while its fixed
point is a saddle (invariance says nothing about interior dynamics). A flow can
have an asymptotically stable fixed point and still leave `C` transiently
(stability is local). (P3) is a statement about a matrix spectrum, not about a
state at all.

This project's problem is **(P1)**. Schuck's result is **(P2)**. TDDM's
"practically positive definite though the non-Hermitian Hamiltonian matrix
obtained from TDDM does not guarantee it" (nucl-th/0305088) is **(P3)**.

### 1.1 The one place they genuinely couple

(P1) and (P2) coincide **iff the stationary point lies on ∂C and the unstable
manifold is transverse to ∂C.** Then an unstable direction pointing out of the
cone *is* a representability violation.

This is not a hypothetical for us. **The HF reference is on ∂C, deeply.** For a
single determinant, γ is idempotent and

  ²D_{pq,rs} = γ_{pr}γ_{qs} − γ_{ps}γ_{qr}

has support only on occ⊗occ; ²Q only on virt⊗virt; ²G has its own kernel. So
²D_HF is not merely a boundary point of the DQG cone, it is an extreme point
with a very large kernel. **Any perturbation theory built on HF starts on the
boundary.** The first-order response must be tangent or inward in every null
direction, or the theory is non-representable at the first order it acts.

This is the exact rank-2 analogue of the rank-1 argument in `01_derivation_1rdm.md`
§4.2 (nₚ = 0 and nₚ = 1 are the boundary; the MK hole factors make the flux
vanish there). It is also, plausibly, why (P1) and (P2) kept getting tangled in
practice: about HF they are evaluated at the same point.

### 1.2 Why this is the cheap test, not an SOS program

I earlier proposed a Positivstellensatz/SOS feasibility program over free term
coefficients. That is still the general form, but it is the wrong *first* move —
the degree blows up and it answers a global question we do not yet need.

The local version at the HF point is nearly free:

- Let `L = DF(²D_HF)` be the linearization of the candidate EOM at HF.
- **Boundary flux (P1, necessary condition):** for every v ∈ ker(²D_HF),
  ⟨v| L[δ] |v⟩ ≥ 0, and likewise for ker(²Q_HF), ker(²G_HF). Since ²Ḋ is linear
  in the perturbation at first order, and the EOM is linear in the free
  coefficients c_i, **these are linear inequalities in the c_i.**
- **Pole structure (c):** the eigenvalues of the same `L`.

One object, both questions, and the positivity half is a linear feasibility
problem rather than a semialgebraic one. Infeasibility here is already a
theorem — no closure of that form is representable at first order about HF.

Only if this passes does the global cone-invariance question (and SOS) arise.

### 1.3 Consequence for the reference state

The MP2 2-RDM is a known non-N-representable object. EE2's `CorrelateInitial()`
offers exactly the knobs JP recalls turning — `GML` (adiabatic Gell-Mann–Low),
`FD` (Fermi–Dirac at `Temp`), `Mp2NOGuess` (MP2 natural orbitals, `=2` for
diagonal only), and `InitialCondition==1` (HOMO→LUMO). If the reference itself
sits outside `C`, no property of the generator can rescue the trajectory. So the
reference must be screened for representability *before* the generator is
blamed. That is a one-line check we currently do not perform.

---

## 2. Pole-structure target

**ADC(2), not CIS(D).** Confirmed against primary abstracts:

- Head-Gordon, Rico, Oumi & Lee, *Chem. Phys. Lett.* **219**, 21 (1994) —
  CIS(D) is "a second-order perturbation expansion of the coupled-cluster
  excited state method … in which CIS is zeroth order, and the first-order
  correction vanishes." Non-iterative, size-consistent, N⁵.
- Head-Gordon, Oumi & Maurice, *Mol. Phys.* **96**, 593 (1999) — CIS(D₀)
  (quasi-degenerate) and CIS(D∞) (iterative).
- Schirmer, *Phys. Rev. A* **26**, 2395 (1982) — ADC(n) constructed so the
  secular problem reproduces the diagrammatic expansion of the polarization
  propagator through order n; Hermiticity and size-consistency structural.
- Trofimov & Schirmer, *J. Phys. B* **28**, 2299 (1995) — working equations.
  Order structure: 1p1h/1p1h through 2nd, 1p1h/2p2h through 1st, 2p2h/2p2h 0th.
- Hättig, *Adv. Quantum Chem.* **50**, 37 (2005) — **CIS(D∞) ≡ partitioned
  ADC(2)**; CIS(D) is its single non-iterative first iteration seeded with the
  CIS eigenvalue. Also: ADC(2) is the Hermitian analogue of CC2. **[unverified
  against primary text — abstract-level only]**
- Dreuw & Wormit, *WIREs Comput. Mol. Sci.* **5**, 82 (2015) — review.

So CIS(D) has no pole structure of its own: it shifts a CIS root rather than
diagonalizing anything with second-order poles, and degrades for quasi-degenerate
states. ADC(2) is the Hermitian second-order reference to match.

**Caveat on the reference point.** The poles of a linearized EOM depend on where
you linearize. Linearize at HF → RPA-like; linearize at a correlated stationary
²D → ERPA/SCRPA-like, with ground-state correlation in the poles. ADC(2) is
built on an MP-partitioned reference. So the pole comparison is only meaningful
once the reference is pinned down — which is JP's point (d), and is why he was
iterating over ground-state choices in the original code.

### 2.1 A test that does *not* work

Chatterjee & Pernal, *J. Chem. Phys.* **137**, 204109 (2012): ERPA from Rowe's
EOM with matrix elements written entirely in terms of the ground-state 1- and
2-RDM. Verified from the abstract: "Applying the ERPA equations with the exact
2-RDM to a hydrogen molecule reveals that the resulting ¹Σg⁺ excitation energies
are not exact" — the *excitation-operator manifold* is incomplete, not the RDM.
Their ERPA2 adds selected doubles and is exact for two-electron singlets; the
APSG (geminal) 2-RDM gives excellent singlets for Be, LiH, and CCSD-quality
curves for N₂.

Consequence: "feed in the exact 2-RDM, check the poles" is **not** a valid
correctness test. Pole error and RDM error are separable failures.

This is also the closest prior art to the geminal idea: APSG 2-RDMs are used
inside ERPA. I found no work propagating a 2-RDM in a geminal parameterization
in real time.

---

## 3. N-representability during propagation — state of the art

- Lackner, Březinová, Sato, Ishikawa & Burgdörfer, *Phys. Rev. A* **91**, 023412
  (2015), arXiv:1411.0495 — closed 2-RDM EOM with a 3-RDM reconstruction
  preserving norm, energy, spin. "Approximately enforcing N-representability
  during time evolution is essential for achieving stable solutions."
  Enforcement is **purification each step**, not a structurally safe EOM.
- Lackner et al., *Phys. Rev. A* **95**, 033414 (2017), arXiv:1611.00888 — HHG.
- Donsa, Lackner, Burgdörfer, Bonitz, Kloss, Rubio & Březinová,
  *Phys. Rev. Research* **5**, 033022 (2023), arXiv:2303.04576 — 1D Fermi-Hubbard
  quench testbed against exact references.
- **Pescoller, Eder & Březinová, arXiv:2412.13566 (2024)** — projective
  purification: restores a chosen set of N-representability conditions,
  maintains contraction consistency across orders, and preserves all conserved
  quantities, minimally invasively. This is the principled replacement for EE2's
  diagonal clamping + `Posify()`, and because the correction is minimal by
  construction, **its size becomes a meaningful error metric**.
- Mazziotti, *Phys. Rev. E* **65**, 026704 (2002) — alternating eigenvalue
  corrections in particle/hole/particle-hole metrics; *Chem. Rev.* **112**, 244
  (2012) — DQG cone and SDP framing; *PRL* **97**, 143002 (2006) — ACSE.

**Nobody has a structurally positivity-preserving rank-2 EOM.** Every working
code purifies.

### 3.1 The D-vs-QG obstruction

Worth stating plainly because it kills the naive square-root idea:
²D = LL† buys **D** for free and says nothing about **Q** and **G**, which are
positivity conditions on *different linear images* of the same object. It is a
cone-projection problem, not a factorization problem — which is exactly why the
Vienna group projects rather than factorizes. Likewise, writing the EOM in
Lindblad form on the two-particle space buys D-positivity, not N-representability;
the latter is strictly stronger. I found nothing closing that gap.

---

## 4. Particle–hole symmetry ↔ simultaneous D and Q preservation

**No paper states or proves this.** The adjacent facts are standard — Q is the
two-hole conjugate of D (Coleman, *Rev. Mod. Phys.* **35**, 668 (1963);
Garrod & Percus, *J. Math. Phys.* **5**, 1756 (1964)) — so "a ph-conjugation-
covariant term list preserves both" is the natural corollary, but it is
unwritten. The only prior art in spirit is JP's own rank-1 result
(*J. Chem. Phys.* **142**, 134113 (2015)), where η = 1 − γ makes 0 ≤ γ ≤ 1
invariant and yields Fermi–Dirac rather than Boltzmann.

Treat as an **original conjecture of this project** and a Lean target.
**G is a different linear image and needs its own argument.**

Supporting datum, from this repo: `qdiss/symbolic/phtest.py` shows the 96-term
1-2-2 triple commutator (MKNO cell 1624) is *exactly* antisymmetric under
γ ↔ η, phase −1, all 64 monomials. So the un-secularized derivation already has
the covariance the conjecture would need.

---

## 5. TDDM ≠ MK cumulant truncation

TDDM closes by approximating the 3-RDM as an **antisymmetrized product of
2-RDMs** (arXiv:1312.3703), later patched for a normalization effect
(arXiv:1708.00569), and tested further in arXiv:1904.06780. That is a different
closure from λ₃ = 0. Running both on the same electronic Hamiltonian would
separate MK-specific failure from generic rank-2 truncation failure.

Also from this literature, cheap and worth adopting as standing checks:

- **Spurious-mode test** (Tohyama & Schuck, arXiv:nucl-th/0305087): conditions
  for the translational spurious mode and its double phonon to sit at exactly
  zero — no truncation of the single-particle space, and retention of *all*
  components of the one-body (and, for the double mode, two-body) amplitudes.
  Same conditions give the Kohn theorem and the continuity equations. They also
  flag the non-Hermiticity problem of HF-based ERPA and propose a Hermitian
  density-matrix formulation. This is a (P3) test and a good regression check.
- **Do not prune the term list by norm** (arXiv:1708.07219): small-norm pp/hh
  and non-2p2h amplitudes matter at strong coupling.

Overview to read first if reading one thing: Schuck, Delion, Dukelsky, Jemai,
Litvinova, Röpke & Tohyama, *Phys. Rep.* **929**, 1 (2021), arXiv:2009.00591.

Small-amplitude limit, for the ground-state workstream: STDDM (arXiv:nucl-th/0311016)
is "a version of extended RPA theories with ground-state correlations";
STDDM*-b's one-body sector is equivalent to Self-Consistent RPA
(Schuck & Tohyama, *Phys. Rev. B* **93**, 165117 (2016), arXiv:1603.08276).

---

## 6. What this changes about the plan

1. **Positivity and ground-state stability are separate workstreams.** Do
   positivity first (§1.2), since that is the blocker.
2. **The first concrete computation is the HF boundary-flux condition** — linear
   inequalities in the free term coefficients, not an SOS program.
3. **Screen the reference state for representability** before blaming the
   generator.
4. **Falsify poles against ADC(2)**, at a specified reference point.
5. **Do not use "exact 2-RDM in ⇒ exact poles out"** as a correctness test.
6. **Model systems: electronic Hamiltonians via PySCF (or BAGEL).** Lipkin is
   not needed — the earlier argument for it rested on the (P2) reading now
   retracted.
