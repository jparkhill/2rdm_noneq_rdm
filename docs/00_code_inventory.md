# Code inventory — EE2 / dissipative RDM dynamics

All sources cloned to `/root/workspace/refs/`. Every repo is a single squashed
commit (no branch history to mine). Collected 2026-08-02.

| repo | HEAD | date | what it is |
|---|---|---|---|
| TightBindingTCL2 | 776fe1d | 2016-01-28 | EE2 as of 2015 — the *tight-binding / standalone* fork |
| RealtimePySCF | e5edfca | 2017-09-11 | EE2 mature — fullest C++ + a decoupled C-ABI fork + Py2 layer |
| MukherjeeKutzelniggNO | 22f3b6d | 2017-12-21 | MKNO.nb — the GNO derivation engine (4.9 MB, 3294 cells) |
| Symbolic-Many-Fermion-Codes | b09e6b4 | 2016-11-10 | CC/λ/DM diagram derivation + autocoding, plus a C++ contraction runtime |
| CorrelatedPolaron | 716355d | 2016-11-10 | Python Wick engine + TCL2 system–bath propagators (2011–12) |
| Transient_Absorption_Paper_SupportingInformation | — | 2017 | GFP3 pump–probe inputs + spectra (data only) |
| SupportingInformation_PerovskiteDynamics | — | 2017 | perovskite Q-Chem input + **production `TDSCF.prm`** |

**`tdscf_pyscf` does not exist as a public repo.** Full listing of the 30 public
repos under github.com/jparkhill contains no such name — it is private, renamed,
or local-only.

---

## 1. TightBindingTCL2 — the 2015 EE2 (the one you linked)

`TCL_EE2.h`, 2048 lines, one `EE2` class. This is the **1-RDM-only** ancestor:
no `r2`, no BBGKY, no spin-orbital V. Its value is that it contains **four
side-by-side attempts at the same correlated `RhoDot`**, which is the historical
record of what you tried:

| routine | lines | character |
|---|---|---|
| `RhoDot` | 819–986 | *live*. 4 intermediates I1–I4, **diagonal/secular** (`int y=x;`), assembled as `Rho·I1·Eta + Eta·I2·Rho + Rho·I3·Rho + Eta·I4·Eta`, signs **+** |
| `RhoDotGood` | 987–1129 | same structure, `Vt` carries the antisymmetrized `2V−V^x`, and assembly signs are **−**. Comment split: I1/I2 "MP2 like", I3/I4 "excited state correlation" |
| `RhoDotYanai` | 1130–1498 | spin-resolved: separate same-spin (`Vs_ss`,`Vt_ss`) and opposite-spin (`Vs_os`,`Vt_os`) channels; explicit ρρρ cubic terms in the pp block |
| `RhoDotSec` | 1499–1683 | "most up-to-date secular term". RI-driven: `Vs`, `V13_1`, `V14_1`, `V23_2`, `V24_2` built as `kron(B,B)` over aux index, then 7 intermediates I1–I7. **Each line is annotated with the MKNO term it came from** |
| `RhoDotSecFull` | 1745–2040 | the unfactorized reference, marked "dead/nonworking/notmaintained" from line 1684 |

The `RhoDotSec` / `RhoDotSecFull` comment annotations are the Rosetta stone:

```
//61{{Sign, -1}, {{Vs,2},{p,u},{p,r}}, {{Vt,2},{p,s},{p,t}},
//   {{γ,1},{p},{p}}, {{γ,1},{u},{t}}, {{η,1},{p},{q}}, {{η,1},{s},{r}}}
Out1(p,q)+=(-1.)*Vs[...]*V24_2[...]*Rho(p,p)*Eta(p,q);
```

i.e. every hand-coded line is stamped with the machine-derived MKNO term, in
exactly the `{Sign, {tensor,rank},{upper},{lower}}` format my parser now reads
out of MKNO.nb. **This gives a direct, checkable derivation → code mapping.**

Other content: `TCLprop.cpp` is a standalone `main()` on a linear polyene
tight-binding H (`BuildLinearPolyene`), with initial conditions `FermiDirac`,
`Boltzmann`, `Gs`, `Es`, `Left`, `CoherentEven`, `CoherentLeft`, `Infinite` — a
ready-made numerical test bed. **The file has unresolved git conflict markers**
(two different `main()`s merged badly), as does `makefile` and `TCL_prop.h`.
`run` (8196 lines) is a captured log of a PYXAID-driven perovskite job.

## 2. RealtimePySCF — the mature EE2

`lib/tdscf/TCL_EE2.h`, **3743 lines**, banner `EE2 JAP 2015`, self-described
"Does cumulant decomposed BBGKY2". Five propagators dispatch from `step()`:

| flag | routine | what |
|---|---|---|
| default | `RhoDotCorr` | correlated 1-RDM EOM (nonlinear TDHF correction) |
| `2Born` | `RhoDot2B` | second Born |
| `TDTDA` | `tdtda`/`sotdtda` | time-dependent Tamm–Dancoff |
| `BBGKY` | `bbgky1`+`bbgky2` | **the 2-RDM propagation** |
| `StaticStep` | — | fixed MP partitioning, "based off what I learned about stability of BBGKY" |

**The 2-RDM path.** `r2` is a full spin-orbital 2-RDM, `cx_mat(4n²,4n²)`.
`Separabler2()` initialises it as the antisymmetrised Grassmann product ρ∧ρ.
`BBGKYstep` propagates ρ and `r2` with separate RK4s inside a shared MMUT
exponential; `Transform2` rotates `r2` into/out of the interaction frame (N⁵
sequential contraction). Energy from `CorrelatedEnergy(r2)` — reduced-Hamiltonian
trace with the (2N−1) one-body distribution factor.

**The representability machinery is all present and all disabled or crude:**
- `RhoDotCorr` opens by hard-clamping the diagonal
  (`if(real(Rho(i,i))<0) Rho(i,i)=0; if(abs(Rho(i,i))>1) Rho(i,i)=1;`) —
  rank-1 representability by fiat, at every RK4 substep.
- `TCL_prop.h::Posify()` (line 1176) — eigen-decompose ρ−ρ_b, clamp eigenvalues
  into [0,1], rescale electron/hole branches to match, re-add missing charge to
  the HOMO- or LUMO-like natural orbital. Gated by `Stabilize` (default 2).
  Called on `Rho` *and* `RhoM12` after every EE2 step.
- **`TracedownConsistency()`** (line 2455) — force-antisymmetrises `r2`,
  rescales it to the correct trace, extracts an effective *n*, and checks the
  cumulant traces to γ·η. **Commented out of the propagation loop** at line
  2173, with the note *"Ensures the separable part = Rho at all times."*
  This is the single most interesting artifact in the codebase: the ²D↔γ
  contraction-consistency constraint was implemented, then switched off.

**Two hand-derived term lists in `RhoDotCorr`.** The live `if(1)` block (3323)
has 24 terms in `Oute/Outr/Oute2/Outr2` with no direct `Out(p,q)` terms; the
dead `if(0)` block (3423) has ~120 terms including 50 direct `Out(p,q)`
contributions. Same structure, one a strict superset — a full derivation and a
reduced one, with the reduced one shipping.

`TCL_EE2.hfullb` (1804 lines) is a third variant: same class, `RhoDot` written
as a clean 5-index intermediate cascade (`tmp/tmp2/tmp3` vs `Tmp/Tmp2/Tmp3` for
the γ- and η-contracted channels) with an explicit inhomogeneity branch
(`g0flag==3` feeds only `rho0`) and an `IncMP2` option. Much the most readable
statement of the 1-RDM EOM in the tree.

**Everything above is Q-Chem-coupled** (`rem_read`, `FileMan`, `formpqr`,
`../cdman/rimp2grad.h`) — none of it compiles standalone.
`lib/tdscf/tdscf.cpp` (1033 lines) is the **decoupled fork**: same physics,
plain C-ABI entry points taking raw pointers, needing only Armadillo. That is
the realistic template for a rebuild.

Python layer (`tdscf/`, Py2 + TensorFlow 0.x, effectively dead) is nonetheless
the clearest spec: `tdcis.py::BBGKYstep/bbgky1/bbgky2` states the whole 2-RDM
step in ~80 lines of `np.einsum`, including the `BBGKY==2` branch that derives ρ
from ²D by `einsum("prqr->pq", rho2)/(n_e−1)`.

## 3. MukherjeeKutzelniggNO — the derivation engine

`MKNO.nb`, 3294 cells, 99 headings. Non-recursive index-based parser now works
(`tools/nbparse.py`, `tools/nbdump.py`, `tools/extract_terms.py`); full markdown
render at `docs/MKNO_dump.md`.

Structure: brute-force permutation implementation of MK generalized normal
ordering → cumulant decomposition of Γ to rank 4 → 1-particle triple commutator
→ accelerated >4-particle versions → a long series of secular approximations
(`Secular8`, `SecularTwo`, "the diagonal/secular approximation I like best",
spin adaptation, Armadillo autocoding) → `Secularizer (which actually works)` →
**`Derivation of BBGKY Heirarchy up to two`** and the 1-2-2 / 2-1-1 triple
commutators for correlated dynamics with and without projectors.

13 cells hold machine-derived term lists; term extraction works and yields clean
Python data (`docs/MKNO_terms.json`):

- cell **1624** (1-2-2 triple commutator, no projectors) → **96 terms**,
  signature δ·δ·γ·γ·γ ×64 and δ·δ·δ·γ·γ ×32
- cell **1578** ("If I do the diagrams by hand I get") → the hand check of 1624
- cells 417/426/430/448 → the secular approximations, in V·γ·γ·η·η form
- cell 430 is the one that matches the `RhoDotSec` code comments term-for-term

Margin note in the BBGKY section: *"comparing with Tohyama"* — Tohyama & Schuck
is the literature comparison point.

## 4. Symbolic-Many-Fermion-Codes — the autocoder

Notebooks that diagrammatically derive **and emit C++ for** coupled cluster of
arbitrary rank, plus λ-multipliers and density matrices:

- `CC2To6LamDm.nb` (18850 cells) — Multipliers / Algebra / Doubles…Hextuples /
  **Density Matrices** / "Programs for all DM blocks"
- `WholeCC6Update.nb` (4646 cells) — Algebra + Code
- `CC246multipliers.nb`, `CCDiag2410.nb` (contraction planning, `AutoProgramDsk`
  writes a populate-file readable from disk), `CCSDDisk.nb`, `BergmanPath.nb`

`ContLibrary/` is the C++ runtime the notebooks emit against (John Parkhill
2008): `MGC_GenTensor.h` (`OrbitalTypes`, `GenTensor`), `MGC_TensorSymmetry.h`
(`TensorSymmetry` — permutational index symmetry), `MGC_Contraction.h` 2151
lines (`ContractionMapSource`, `ExpandedContractionMap`, `Contraction`,
`Contractions` — symmetry-aware contraction planner with cached maps).

**This is the piece that matters for the rebuild**: a working precedent for
going derived-terms → executable contraction code, including the DM blocks.

## 5. CorrelatedPolaron — the Python Wick engine

`Wick.py`, 134 KB. Fermion Wick's theorem for **both genuine and Fermi vacua**
(`AntiCommuteGenuine`, `AntiCommuteFermi`), `NormalOrder` / `NormalOrderDiag`
with diagrammatic sign and 1/2^{n_eq} factor (`DiagrammaticSignAndFactor`),
term filters (`ConnectedPart`, `NoSelfContractPart`, `QPRankPart`, `ApplyTDA`,
`QPart` Q-space projector), `Secularize`, and numerical evaluation via
`numpy.tensordot` + `scipy.weave` inline C++.

Note: this engine is **ordinary** normal ordering, not MK generalized normal
ordering — no cumulants, no λ. The GNO lives only in MKNO.nb.

`Propagatable_TCL2_Herm.py` is the TCL2 system–bath propagator with explicit
switches (`QSpaceProjector`, `HermitianTheory`, `TammDancoff`, `OnlyPH`,
`ForceAntiHermitian`, `SecondOrderDipole`). `TensorNumerics.py` carries
`Params.Adiabatic` 0–4 (numerical/analytical ΔS integration, Markov).

Dead as code: Python 2, `scipy.weave` (removed in SciPy 1.0), `izip`,
`__builtin__`. Alive as specification.

## 6. Production parameters (the SI repos)

`SupportingInformation_PerovskiteDynamics/TDSCF.prm` is a real production run —
worth keeping as the reference for what settings were actually used:

```
Stabilize 2        ActiveSpace 1     ActiveCutoff 0.1   FirstOrb 262  LastOrb 347
UpdateFockMatrix 0 TCLOn 1  Lamb 0   DiagonalHole 0     Temp 300.0
dt 0.02  MaxIter 3000000  tolmax 1e-5  FockBuildThresh 1e-4  RateBuildThresh 1e-4
```

Note `Stabilize 2` — i.e. the full `Posify` positivity projection was on in
production. The published dynamics relied on it.

---

## What this means for the rebuild

1. The **derivation** is recoverable from MKNO.nb as structured data — done,
   term extraction works.
2. The **derivation→code mapping** is recoverable from the `//61{{Sign,…}}`
   comments in `RhoDotSec` — a free regression test for any new autocoder.
3. The **2-RDM propagation** exists and runs, but its N-representability was
   enforced by (a) diagonal clamping, (b) `Posify` eigenvalue projection, and
   (c) `TracedownConsistency`, which was written and then disabled. Those three
   are exactly the failure modes to attack.
4. The **rebuild target** is `tdscf.cpp`'s C-ABI shape or a pure-NumPy/PySCF
   port of `tdcis.py::BBGKYstep`; PySCF 2.14 is installed and working.
