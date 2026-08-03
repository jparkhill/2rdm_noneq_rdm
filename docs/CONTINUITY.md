# Continuity — read this first

This file exists so the project survives the loss of any particular machine,
session, or assistant context. If you are a new agent or a returning human with
no memory of this work, everything you need is here or linked from here.

**Repo:** `https://github.com/jparkhill/2rdm_noneq_rdm.git`
**Principal:** John Parkhill (`github.com/jparkhill`), formerly Prof. of
Theoretical Chemistry, Notre Dame; now Austin, TX.

---

## 1. The problem in one paragraph

Find a **perturbative two-particle reduced-density-matrix (2-RDM) equation of
motion** — dissipative, second-order time-convolutionless (TCL2), derived with
**Mukherjee–Kutzelnigg generalized normal ordering** — that **does not violate
N-representability**. John derived exactly this at rank 1 in
*Nguyen, Nanguneri & Parkhill, "How electronic dynamics with Pauli exclusion
produces Fermi–Dirac statistics," J. Chem. Phys. **142**, 134113 (2015)*, where
the one-electron truncation thermalizes correctly to Fermi–Dirac. He attempted
the same at rank 2 (a method he called **EE2**), never published it, and
abandoned it — the 2-RDM kept leaving the N-representable set. This project
retries that, with the explicit intent that every step be machine-verifiable.

**Why this problem:** it is chosen to be *verifiable in nature*, so an agent can
check its own work. Verifiability is the point.

What is wanted is a **verifier, not a falsifier** — earlier drafts of this file
said the latter and it undersold the target. A candidate EOM has to be judged on
three independent criteria, two of which are positive statements:

1. **Positivity.** DQG eigenvalues stay `≥ 0` under propagation; on the boundary
   `boundary_flux() ≥ 0` (`qdiss/representability/dqg.py`).
2. **Stability, with the right fixed point.** Bounded trajectories, no runaway
   occupations, a monotone Lyapunov function, and thermalization to the correct
   distribution — Fermi–Dirac in the rank-1 case.
3. **Correct perturbative limit.** In the appropriate limit the EOM must reduce
   to a perturbatively corrected excited-state theory and reproduce known
   excitation energies. This is the sharpest constraint on the derivation,
   because it is sensitive to coefficients and not merely to signs: a term list
   can preserve positivity and still be the wrong physics.

A refutation is cheap and one negative eigenvalue is enough to kill *one
candidate*. Confirmation requires all three.

But note the scope: refuting a single closure is a data point, not a result about
the programme. docs/02 lists at least five independent knobs -- the reconstruction
functional for `lambda_3`, where the secular approximation is applied, whether the
dissipator is put in manifestly completely-positive form on the two-particle space
*before* contraction, whether rates are symmetrised so KMS holds term by term, and
contraction-consistency constraints tying `^2D` to `gamma`. **Any** EOM that
manipulates at most 4-index quantities and can be related to a perturbative
correlated correction is of interest. So showing that `lambda_3 = 0` is
non-positive establishes very little; a real no-go has to quantify over the class.

And the strong prior from rank 1 points the other way: manifest CP form plus
term-by-term KMS is *why* rank 1 worked. If the closed rank-2 equation can be
arranged to remain the contraction of a genuine Lindblad flow on the two-particle
space, positivity is structural -- a consequence of complete positivity -- rather
than something to be checked case by case. That is the direction with a mechanism
behind it, not a search.

## 2. Start here, in this order

1. `./scripts/bootstrap.sh` — installs deps, clones the upstream reference
   repos into `refs/`, regenerates derived artifacts, runs the self-check.
2. `docs/01_derivation_1rdm.md` — the rank-1 result, re-derived from scratch and
   verified numerically. **This is the solved case and the template.** The key
   section is §4.2, the invariant-region argument.
3. `docs/02_project_plan.md` — the rank-2 generalization being searched for, the
   boundary-flux criterion, and the pipeline.
4. `docs/00_code_inventory.md` — what is in each of the seven upstream repos and
   what it implies for the rebuild.

## 3. The one idea everything hangs on

At rank 1 the N-representable set is just `0 ≤ γ ≤ 1`. The MK hole contractions
`⟨a_p a†_p⟩ = 1 − γ_pp` put a factor `γ` on every gain term and `(1 − γ)` on
every loss term, so **the flux through the boundary vanishes identically on that
boundary** and the box is invariant. Representability is not enforced; it is
structural.

Numerically confirmed: keeping the hole contractions gives Fermi–Dirac to 1e-10
and `0 ≤ γ ≤ 1` throughout. Dropping them gives Boltzmann and drives one
spin-orbital to **2.81 electrons**. Same mechanism, both results.

**The rank-2 question is whether an analogous structural argument exists** for
the DQG cone: for `²D` on the boundary with null vector `|v⟩`, does
`⟨v| Ṁ[²D] |v⟩ ≥ 0` hold for `M ∈ {D, Q, G}`?

That criterion is **Nagumo's theorem** specialized to the PSD cone, and naming it
matters because it says what we do and do not need. Nagumo: a closed convex `K`
is forward-invariant under `ẋ = f(x)` iff at every boundary point `f(x)` lies in
the tangent cone — equivalently `⟨ν, f(x)⟩ ≤ 0` for every outward normal `ν`.
The outward normal to the PSD cone at a matrix with null vector `v` is `−vv†`,
which gives the criterion above. **It is a pointwise algebraic condition on the
right-hand side.** We therefore need only the EOM, never a propagator, and the
rank-1 argument above is the same theorem in the scalar case.

Two things follow. The primary attack is **structural**: ask whether every term
of `Ṁ` carries a factor that vanishes when `M(²D)v = 0`, which is an identity in
abstract indices and so independent of `n`. Refutation is **small-`n` numerics**:
build a boundary `²D` at `n = 4, 6` and look for a negative. SOS/Positivstellensatz
is a *late, optional* tool — the SDP is `O(n⁴)` in the RDM elements and does not
scale, so it is worth reaching for only after a structural reduction to a few
scalar invariants, or for a machine-checked certificate at fixed small `n`.

Two caveats. DQG is a **relaxation** — necessary but not sufficient for
N-representability, so invariance of DQG does not prove the 2-RDM is physical.
And at non-smooth boundary points where several constraints are simultaneously
active the normal cone has several generators, and the condition must hold for
all of them, including across `D`/`Q`/`G`; a check that tests one at a time can
pass while the dynamics still escapes.

## 3b. How the EOMs are derived (do not reinvent this)

The physical route, and the only one to use:

1. Write the **Heisenberg-picture EOM for an operator of the rank you want to
   propagate** -- for the 1-RDM that operator is `a^p_q`, for the 2-RDM it is
   `a^{pr}_{qs}`. Not a generic product of operators: the operator whose
   expectation value *is* the quantity being propagated.
2. **Take the expectation value.** That converts an operator equation into an
   equation for the RDM element. It also kills every residual `{a_k}`, since
   `<{a_k}> = 0` for `k >= 1` by construction of GNO.
3. **Return to the Schrodinger picture** for the quantity being propagated.

The closure enters at step 2 and nowhere else: taking the expectation of the
TCL2 dissipator raises RDMs of higher rank than the one being propagated, and
whatever removes them is the closure. At rank 1 the dissipator is quartic in
fermion operators, because it is a product of two one-body couplings `A_alpha`,
so the `gamma` equation picks up `<a+ a+ a a> = gamma_2`; `lambda_2 = 0` closes
it. See docs/01 sections 2 and 3.

The corollary, learned the hard way: do **not** start by building a generic
rank-1 x rank-2 x rank-2 operator product and then look for its scalar part.
That object is not any EOM. Start from the rank of the propagated quantity and
let the Heisenberg equation tell you which commutators appear.

## 4. Notation conventions John uses (match these)

- **η** is the blocking factor / one-hole RDM, substituted via
  `δ^a_b = η^a_b + γ^a_b`.
- **η is never propagated.** It is always reset from γ. Preserve this — it is a
  deliberate choice, not an oversight.
- `λ` is the cumulant; `λ₂ = 0` is the rank-1 closure.

## 5. What is in this repo

```
docs/     the written record — read 01 and 02 before writing code
qdiss/    the new implementation
  rdm1.py            rank-1 dissipative EOM (working, verified)
  symbolic/          MK generalized normal ordering engine (to build)
  dynamics/          propagators (to build)
  representability/  DQG monitors and projections (to build)
tools/    Mathematica .nb parsers — these work, do not rewrite them
  nbparse.py         non-recursive index-based cell parser
  nbdump.py          whole-notebook -> readable markdown
  mmaparse.py        box-expression -> Python structures
  extract_terms.py   derived tensor terms -> structured data
lean/     formalization (empty; rank-1 boundary-flux theorem is the warm-up)
tests/    pytest
scripts/  bootstrap.sh
refs/     upstream repos — NOT committed, recreated by bootstrap.sh
```

## 6. Upstream sources and what each is for

All under `github.com/jparkhill`. Cloned to `refs/` by the bootstrap script.
Full detail in `docs/00_code_inventory.md`.

| repo | why it matters |
|---|---|
| `RealtimePySCF` | **the mature EE2.** `lib/tdscf/TCL_EE2.h`, 3743 lines, five propagators. The only place a real 2-RDM (`r2`, spin-orbital, `4n²×4n²`) is propagated. |
| `TightBindingTCL2` | the 1-RDM ancestor. Five side-by-side attempts at the correlated `RhoDot` — a record of what was tried and rejected. |
| `MukherjeeKutzelniggNO` | `MKNO.nb`, the derivation engine. 3294 cells. Parsed by `tools/`. |
| `Symbolic-Many-Fermion-Codes` | `CC2To6LamDm.nb` (18850 cells) covers density matrices and "Programs for all DM blocks"; `ContLibrary/` is the C++ contraction runtime it emits against. The autocoder precedent. |
| `CorrelatedPolaron` | `Wick.py`. **Ordinary** normal ordering only — no cumulants. GNO exists only in MKNO.nb. |
| `SupportingInformation_PerovskiteDynamics` | carries a real production `TDSCF.prm`. Says `Stabilize 2`. |

`tdscf_pyscf` **does not exist publicly** — the full 30-repo listing was checked.
Private, renamed, or local-only.

## 7. Findings worth not rediscovering

- **`docs/12`'s stated cause is wrong.** It blamed cell 62's rank-4 cutoff:
  `TildeSubstitution` calls `Print` past rank 4, `Print` returns `Null`, and the
  `Null` multiplies silently into the product. Tested directly:
  `tools/rank5_patch.wl` extends all three cell-62 dispatchers to rank 5 using
  the `ta5`/`Na5`/`NRta5` cell 125 already defines, and rank 5 does go from
  `Null` to a 1670-term list. Re-deriving `rv1v2` with it produces a file
  **byte-identical** to the old one, containing **zero `Null`s**, failing by the
  same amount. So that path is never reached in this product.
- **The live lead is `Signature`.** The build emits
  `Signature::normal: Nonatomic expression expected at position 1 in
  Signature[upper$…]`. The code computes
  `OuterSign = Signature[upper]*Signature[lower]` with `upper = Join[u1,u2]`, so
  `Join` came back unevaluated, meaning it received a non-list. Both stage
  runners open with `Off[Join::heads]`, which suppresses precisely the warning
  that would have named the malformed operand. Re-run with that `Off` removed
  before theorizing further.
- **The failure is specific to three factors.** Both two-factor products verify
  exact (8.6e-16, 1.3e-15); all four three-factor orderings fail (1.09–1.26 on
  typical magnitudes 0.67–0.94). Same code, same oracle, so it is not a
  convention or version artifact. Published work is not implicated: it ran the
  `NoCumulants` branch, and cell 1624's 96 terms match.
- **A sandbox's durable computer is durable only at `/root`.** Wolfram Engine
  installed to `/usr/local` was lost when the container was recreated; the
  entitlement in `~/.wolfram_entitlement` survived. Install into the volume:
  Wolfram now lives at `/root/.local/wolfram/14.1` and the Python venv at
  `/root/.venv`. Same trap applies to pip installs and apt packages.
- **The sandbox cannot push.** `origin` is HTTPS and the credential died with
  the old container. Move commits out with `git bundle` and push over SSH from
  the host, or re-point the remote at SSH with a key mounted.

- **EE2 enforced representability three ways, all crude or disabled.**
  (a) diagonal clamping `if(Rho(i,i)<0) Rho(i,i)=0` at every RK4 substep;
  (b) `Posify()` eigenvalue projection in `TCL_prop.h`, gated on `Stabilize`,
  which production runs had set to `2`;
  (c) `TracedownConsistency()` — the ²D↔γ contraction-consistency constraint —
  **written, then commented out of the BBGKY loop**, with the note "Ensures the
  separable part = Rho at all times." That commented-out line is the most
  interesting artifact in the whole codebase.
- **`RhoDotSec` in `TightBindingTCL2/TCL_EE2.h` stamps every hand-coded line
  with the MKNO term it came from**, in the same `{Sign, {γ,1},{p},{q}, …}`
  format `tools/extract_terms.py` reads out of the notebook. This is a free
  derivation→code regression test. Use it.
- `RhoDotCorr` in `RealtimePySCF` has two term lists: a live `if(1)` block with
  24 terms and a dead `if(0)` block with ~120 including 50 direct `Out(p,q)`
  terms. Same structure, one a strict superset. Determine which is which before
  rebuilding.
- MKNO cell 1624 (1-2-2 triple commutator, no projectors) = **96 terms**, 64 of
  signature δ·δ·γ·γ·γ and 32 of δ·δ·δ·γ·γ. Cell 430's secular form matches the
  `RhoDotSec` code comments term-for-term.
- None of the C++ compiles standalone — all Q-Chem-coupled (`rem_read`,
  `FileMan`, `../cdman/`). `RealtimePySCF/lib/tdscf/tdscf.cpp` is the decoupled
  fork: same physics, plain C-ABI, needs only Armadillo. That is the template.
- The Python layer of `RealtimePySCF` is Python 2 + TensorFlow 0.x. Dead.

## 8. Status

- [x] Rank-1 EOM re-derived from scratch and verified (FD to 1e-10; `0≤γ≤1`
      invariant; Lyapunov monotone)
- [x] Confirmed against the published C++ (`TCL.h` builds its dissipator from
      `HoleRho = I - Rho`)
- [x] All upstream source collected and inventoried
- [x] MKNO.nb parsed; terms machine-extractable
- [x] Symbolic MK engine in Python (`qdiss/symbolic/mkno.py`) — products of GNO
      operators, gamma form and cumulant form, ~200 lines, no GWT rules
      transcribed: it reads the one definition `tools/gno.py` validates, forwards
      to expand a bare string into GNO operators and solved for its `m=0` term to
      go back, with ordinary anticommutation between
- [x] Python `A.B` agrees with the notebook's `r1v1` **term for term, all 30
      terms and signs**, in the cumulant basis. Two independent derivations in
      the same basis, so the convention is pinned and neither carries a private
      blind spot
- [x] Python `A.B.C` verifies as an exact operator identity (1e-15 at M=8 and
      M=10) where **every** Wolfram three-factor stage fails by the size of the
      answer
- [ ] 2-RDM falsifier with DQG monitor (`qdiss/representability/`)
- [ ] Closure search
- [ ] Lean formalization — rank-1 boundary-flux theorem first as warm-up

## 9. Working agreement

- Commit regularly. This repo is the durable record; a local machine is not.
- Prefer refutable claims. If a candidate EOM cannot be falsified numerically in
  minutes, the falsifier is not good enough yet.
- Verification is the deliverable as much as the derivation is.
