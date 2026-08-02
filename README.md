# 2rdm_noneq_rdm

Search for a **perturbative, dissipative two-particle reduced-density-matrix
equation of motion** — TCL2, derived with Mukherjee–Kutzelnigg generalized
normal ordering — that **does not violate N-representability**.

At rank 1 this works and is published: the one-electron truncation thermalizes
to Fermi–Dirac, because the MK hole contractions make `0 ≤ γ ≤ 1` an invariant
region ([J. Chem. Phys. **142**, 134113 (2015)][paper]). At rank 2 it was
attempted, called **EE2**, and abandoned unpublished — the 2-RDM left the
N-representable set. This repo retries it, with every step built to be
machine-checkable.

[paper]: https://doi.org/10.1063/1.4916822

## Quick start

```bash
git clone https://github.com/jparkhill/2rdm_noneq_rdm.git
cd 2rdm_noneq_rdm
./scripts/bootstrap.sh      # deps, upstream refs, derived artifacts, self-check
python3 verify_1rdm.py      # reproduces the rank-1 result
```

Then read **[`docs/CONTINUITY.md`](docs/CONTINUITY.md)** — it is written for
someone arriving with no prior context.

## The idea

The rank-1 N-representability conditions are just `0 ≤ n_p ≤ 1`. Generalized
normal ordering about the *correlated, non-idempotent* γ puts a factor `n` on
every gain term and a hole factor `(1 − n)` on every loss term:

```
ṅ_p = Σ_q [ W_{p←q} n_q (1 − n_p) − W_{q←p} n_p (1 − n_q) ]
```

so at `n_p = 0` the derivative is `≥ 0` and at `n_p = 1` it is `≤ 0`. **The flux
through the boundary of the representable set vanishes on that boundary.**
Representability is structural, not enforced. The same factors give Fermi–Dirac
rather than Boltzmann at the fixed point — one mechanism, both results.

Verified numerically: keeping the hole contractions gives Fermi–Dirac to `1e-10`
with `0 ≤ γ ≤ 1` throughout; dropping them drives a spin-orbital to **2.81
electrons**.

The rank-2 question is whether the analogue holds on the DQG cone — for `²D` on
the boundary with null vector `|v⟩`, whether `⟨v| Ṁ[²D] |v⟩ ≥ 0` for
`M ∈ {D, Q, G}`. That is a polynomial positivity statement: approachable by SOS
certificates numerically, and formalizable in Lean once a candidate survives.

## Layout

| path | contents |
|---|---|
| `docs/` | the written record — start with `CONTINUITY.md` |
| `qdiss/` | implementation: `rdm1.py` (working), `symbolic/`, `dynamics/`, `representability/` |
| `tools/` | Mathematica `.nb` parsers for the MKNO derivation notebook |
| `lean/` | formalization (boundary-flux theorem) |
| `scripts/` | `bootstrap.sh` — rebuilds everything from a bare clone |
| `refs/` | upstream source repos, **not committed**; recreated by bootstrap |

## Documents

- **[`docs/CONTINUITY.md`](docs/CONTINUITY.md)** — orientation, conventions,
  findings, status. Read first.
- [`docs/01_derivation_1rdm.md`](docs/01_derivation_1rdm.md) — the rank-1
  derivation and its numerical verification. The solved case and the template.
- [`docs/02_project_plan.md`](docs/02_project_plan.md) — the rank-2
  generalization, the boundary-flux criterion, the search pipeline.
- [`docs/00_code_inventory.md`](docs/00_code_inventory.md) — what is in each of
  the seven upstream repos.
