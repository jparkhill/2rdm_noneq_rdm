# Standing requirements for the 2-RDM / cumulant-hierarchy work

**Read `docs/STATE-<latest date>.md` first** -- currently `docs/STATE-2026-08-09.md`. It records
what is established, what has been superseded, and the next experiment with its known defects.
Several committed results are special cases later shown to be non-general; the state document says
which, and reading commit messages alone will mislead.

These are program invariants, not preferences. Check work against them before reporting.
John edits this file to steer; if something here conflicts with an instinct, this file wins.

## The goal

Isolate **order-general** conditions on the k-cumulant which, imposed at each rank, give a
hierarchy of closures approximating fermionic dynamics at successively finer detail. Individual
good closures are not the goal. The target regime is **non-equilibrium and mixed states**, so a
condition that only makes sense for a pure ground state is of limited value.

## Before calling any condition order-general

**Test it at two adjacent ranks.** A statement verified only at rank 2 is a rank-2 statement.
Whenever the algebra is cheap enough to do at rank 3 as well, do both in the same pass rather
than reasoning narrowly about `lambda_2` and generalising by assertion.

## Screens to apply BEFORE investing in a closure

1. **Lipschitz.** Reject any closure containing an unregularised inverse of a `gamma_1`-dependent
   operator, or a denominator that can vanish on the physical domain. One-minute test; it
   disqualified both `blocked` (pole surface) and unregularised `pt` (rank-jumping pseudo-inverse).
   Survivors: polynomials, truncated Neumann series, level-shifted inverses with FIXED `sigma > 0`.
2. **Extensivity.** For a system separating into non-interacting fragments the cumulant must be
   block-additive with exactly zero cross-fragment elements, and the within-fragment part must not
   change when a distant spectator is added. Global normalisations fail this.
3. **Perturbative limit.** The reconstruction must carry `V`-dependence, or it cannot reach the
   coherent amplitude block that holds most of the exact cumulant's norm at low temperature.

## Reporting discipline

- **State the approximation.** Name the equation propagated, its state variables, what was
  reconstructed and what was truncated, before any number.
- **Define every symbol on first use**, in explicit index or matrix form, including operations
  like `Tr_1`. No bare letter standing for a group of terms without the next line expanding it.
- **Never infer the objective.** The metrics that matter: pole positions vs FCI and their errors;
  trace errors in the 1- and 2-RDM; positivity errors in the 1- and 2-RDM. Ask if unclear.
- **Perturb the wavefunction, never the RDM.** To probe a spectrum, apply a dipole operator to
  the state and compute the reduced density matrices from the result. Adding a kick directly to
  `gamma_2` produces a starting point that is not N-representable, and since the exact `gamma_2`
  of a pure state can sit exactly on the PSD boundary, the kick alone then dominates every
  positivity and trace metric.
- **Always include a mixed-state row.** Ground-state-only benchmarks are the hardest case for
  positivity and the least representative. Read positivity against the exact clearance
  (smallest eigenvalue of the exact matrix) at that state, not against zero.
- **Keep a clean path** derivation => code => simulation. Two derivation routes (plain Wick, and
  MKNO/GWT) are kept separate and validated against each other and against `qdiss.eom.fullci`.
  Never correct factors in a circle with no ground truth. Report expectation-engine bugs at once.

## Hard rules

- Admissible closures: **cumulant reconstructions derived by the MKNO path**, and
  **Mueller-2405-type** approximations. No outright negation of a moment (`gamma_3 = 0` etc).
- **Never test a 2-RDM EOM on a two-electron system.** Use `N >= 3`, and prefer geometries with
  no spatial symmetry.
- Two different Muellers: `muller2405` (arXiv:2405.05093, a 3-RDM reconstruction) and
  `muller_rdmft` (the RDMFT `sqrt(n_p n_v)` exchange functional). Always suffix.
- Do not say "isospectral" for occupation-conserving. The spectrum means the poles.
- Do not run `docker volume prune`; `qdiss-work` holds a harness volume.
