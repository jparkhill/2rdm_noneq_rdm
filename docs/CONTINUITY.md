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
check its own work. Any candidate EOM is refuted in minutes by propagating and
watching a DQG eigenvalue go negative. Verifiability is the point.

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
`⟨v| Ṁ[²D] |v⟩ ≥ 0` hold for `M ∈ {D, Q, G}`? That is a polynomial positivity
statement — tractable by SOS/Positivstellensatz numerically, and a Lean target
once a candidate survives.

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
- [ ] Symbolic MK engine in Python (`qdiss/symbolic/`)
- [ ] 2-RDM falsifier with DQG monitor (`qdiss/representability/`)
- [ ] Closure search
- [ ] Lean formalization — rank-1 boundary-flux theorem first as warm-up

## 9. Working agreement

- Commit regularly. This repo is the durable record; a local machine is not.
- Prefer refutable claims. If a candidate EOM cannot be falsified numerically in
  minutes, the falsifier is not good enough yet.
- Verification is the deliverable as much as the derivation is.
