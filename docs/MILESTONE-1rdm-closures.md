# Milestone: 1-RDM equations of motion from cumulant reconstruction

Self-contained record of the 1-RDM line of work, so it can be picked up cold. Every claim
names its derivation method and the tool that produced it. Ground truth throughout is
`qdiss.eom.fullci`, exact Fock-space linear algebra.

## 1. What this is, and what it is not

The object is a reconstruction of the two-body CUMULANT, `lambda_2[gamma_1]`, substituted into
the **unapproximated** MKNO rank-1 equation of motion. That equation is a function of `gamma_1`
and `lambda_2` alone -- 2 one-body and 12 two-body terms, factors `D1 x14, G1 x18, L2 x4` --
and contains no `lambda_3`, which is why no closure of the higher rank is needed for it.

This is **not** RDMFT. Nothing minimises an energy functional and nothing is inserted as a
coupling matrix; TD-RDMFT is already developed and is not this niche. The RDMFT functionals
named below supply only the algebraic *shape* of the cumulant.

## 2. Verified infrastructure

`qdiss/eom/` keeps the derivation methods separate and holds each against FullCI alone:

- `hamiltonian.py` -- the single definition of `H`, including `TWO_BODY_PREFACTOR = 0.25`,
  imported by every producer so no derivation can be tuned into agreement with another.
- `wick.py` -- derivation A, ordinary vacuum Wick, bare moments.
- `gwt.py` -- derivation B, MKNO generalized normal ordering; exposes bare and cumulant forms.
- `fullci.py` -- ground truth; consumes no term list.
- `validate.py` -- each derivation against FullCI separately. Worst error `3.3e-15`.

The two derivations produce **identical** term lists (`|bare - GNO| = 0.0`), which validates
the GWT machinery, since thousands of intermediate terms cancel to the same expression.

`tools/physical_hamiltonian.py` builds spin-orbital `h` and antisymmetrised `V` from PySCF and
checks `build_H` against the PySCF FCI total energy: H2/6-31g `4.4e-16`, H4 chain/sto-3g
`2.7e-15`, HeH+/sto-3g `8.9e-16`.

## 3. The closure family (`qdiss/closures.py`)

All but `mp2` are diagonal in the natural-orbital pair basis and fixed by one symmetric
coefficient, with `c = n_p n_v - X`:

    lambda_2[pv,qw] = (1/2) c(n_p,n_v) ( d_pw d_vq - d_pq d_vw )

| name | `X` | note |
|---|---|---|
| `hf` | `n_p n_v` | `lambda_2 = 0`; reproduces TDHF `i[(h+W)^T, gamma_1]` at `2.7e-16` |
| `muller` | `sqrt(n_p n_v)` | Mueller / Buijse-Baerends shape |
| `power` | `(n_p n_v)^alpha` | Sharma et al., `alpha ~ 0.6` |
| `bbc` | Mueller, `c = 0` on the diagonal, root sign flipped for weakly-occupied pairs | BBC1 + Goedecker-Umrigar, **not** BBC3 |
| `blocked` | `A^f + f^A`, `A = gamma_1 - gamma_1^2` | ours; the only member satisfying the trace identity by construction |
| `mp2` | `T + T-dagger`, `T = V n n (1-n)(1-n)/Delta e` | not occupation-diagonal |

`blocked`'s coefficient, in the natural-orbital basis with `a_i = n_i(1-n_i)`, `T = sum a_j`,
`Q = sum_j a_j/(T - 2a_j)`:  `f_i = -a_i / [(1+Q)(T - 2a_i)]`. Both factors go as `a_i`, so each
vanishes as an orbital empties or fills.

`bbc`'s philosophy: Mueller patched by occupation block. The weak-weak sign flip removes
Mueller's spurious enhanced exchange among virtual-like orbitals; the `p = v` diagonal restored
to `n_p^2` removes the self-pair term. In cumulant language it is a sign change between blocks,
`c < 0` for occupied and mixed pairs, `c > 0` for virtual-virtual.

## 4. The partial-trace identity

Every physical `lambda_2` satisfies

    sum_x lambda_2[px,qx] = -(gamma_1 - gamma_1^2) = -gamma_1 eta_1

so the contraction of the cumulant **is** the particle-hole blocking product, and since
`gamma_1 eta_1 >= 0` exactly when `0 <= gamma_1 <= 1`, this is the Pauli bound in cumulant
variables. Verified at `1e-16` on exact states at (4,2), (6,3), (8,4).

`closures.trace_correct` enforces it on any reconstruction: for antisymmetric `W` with partial
trace `S`, `Corr(S) = (1/(M-2))(S^I + I^S) - [tr S/((M-2)(M-1))](I^I)` has partial trace exactly
`S`, so `W - Corr(S_W) + blocked` satisfies the identity while keeping every coherence of `W`.
Measured residuals fall from `1.4e-02 ... 4.4e-01` to `~3e-16`.

The trace-carrying part added is the `blocked` structure -- one particular solution among an
affine family. **The choice is not neutral**: it removes a spurious near-zero mode from `mp2`
and destroys Mueller's `0.86` pole on H2. Trying other particular solutions is an open knob.

## 5. Results

Poles at the exact FCI ground-state `gamma_1`. That reference is used because at a
determinantal `gamma_1` every occupation-diagonal reconstruction vanishes identically and the
`sqrt(n)` ones have divergent derivatives. It is the exact fixed point: the exact `lambda_2`
there gives `|gamma_1-dot| = 4.5e-16`.

**H2 / 6-31g**, exact `[0.3948 x3, 0.5626, 0.8599 x2]`, TDHF at its own HF point
`[0.3599 x3, 0.5363 x3]`:

    blocked   raw | trace 1.9e-16 | 0.3635 0.3635 0.3635 0.4667 0.4825 0.4825
    muller    raw | trace 1.8e-01 | 0.3701 0.3701 0.3701 0.5374 0.8660 0.8660
    bbc       raw | trace 1.8e-01 | 0.3694 0.3694 0.3694 0.5360 0.8615 0.8615
    mp2       tc  | trace 1.9e-16 | 0.3915 0.3915 0.3915 0.4014 0.4051 0.4053

TDHF misses the `0.86` manifold entirely, placing a spurious triple at `0.5363`; the
Mueller-family reconstructions recover it (`bbc` gives `0.8615` against `0.8599`).
Trace-corrected `mp2` gives the best low pole, `0.3915` against `0.3948` -- 0.8% against TDHF's
8.8%.

**H3 linear / sto-3g**, exact `[0, 0.0777 x2, 0.1251 x3]`, TDHF `[0.0122 ... 0.0589]`:

    bbc       tc  | trace 4.4e-16 | 0.0247 0.0774 0.0874 0.1300 0.1390 0.1610

`0.0774` against `0.0777`, and `0.1300` for a manifold TDHF lacks entirely.

**Stationarity predicts pole quality.** `|gamma_1-dot|` at the FCI density measures how badly a
reconstruction fails to hold the exact fixed point. On H2 `mp2_tc` is most stationary
(`0.0038`) and has the best low pole; on H3 `bbc_tc` is most stationary by tenfold (`0.0027`
against `0.027` for `lambda_2 = 0`) and has the best poles.

## 6. Why occupations stay frozen

The rank-1 EOM splits into 8 two-body terms carrying no `lambda_2` plus 2 one-body terms --
exactly TDHF, a commutator, which moves **no** occupation at any time (`7e-16` throughout a
propagation) -- against 4 terms carrying `lambda_2`, which supply all occupation motion. The
blocking factors are present and work perfectly; they simply cannot constrain `lambda_2`, which
is an independent tensor rather than a product of `gamma` and `eta`.

Any `lambda_2` that is a **matrix function of `gamma_1`** commutes with it, is diagonal in the
natural-orbital pair basis, and contributes nothing to the diagonal of `gamma_1-dot`. This is
the known TD-RDMFT statement that adiabatic occupation-number functionals leave occupations
frozen.

`mp2` is not occupation-diagonal, yet also freezes them (`8.9e-17`). `tools/cumulant_phase.py`
locates the reason -- splitting the amplitude into its two Hermitian combinations:

    in phase       0.25 (T + T-dagger)    max|ndot| 8.9e-17
    out of phase   0.25 i(T - T-dagger)   max|ndot| 3.3e-01
    exact lambda_2                        max|ndot| 2.6e-02

Both are equally pair-Hermitian and built from the same `T`, and only the out-of-phase one
transfers population. `mp2` picks the combination appearing in the static MP2 density. The
exact cumulant overlaps both (`0.052` in phase, `-0.049` out of phase).

So the no-go is **not** "adiabatic implies frozen": `i(T - T-dagger)` is adiabatic and moves
occupations strongly. It is that occupation transfer lives entirely in the component out of
phase with the interaction, and nothing instantaneous fixes that component's magnitude -- it
overshoots thirteenfold. The physical value is the cumulant's lag behind the density, i.e.
retardation, which requires propagating `lambda_2`.

## 7. Retracted along the way

- "A reconstruction must depend on `V` to move occupations, so `mp2` can." False; `mp2`
  freezes them. The earlier supporting measurement used a **complex random** `V`, Hermitian but
  not real-symmetric, which smuggles in an out-of-phase part no molecular `V` has. On real
  integrals a raw `V` perturbation gives `4.4e-17`.
- "Occupation-conserving means it cannot beat TDHF." False. The `lambda_2` contribution is
  purely off-diagonal, the generator differs by 19% (`|X_blocked - X_TDHF| = 0.50` against
  `|X_TDHF| = 2.64`), and the poles differ. Same occupations, different theory.
- Never say "isospectral" for preserving occupations. The spectrum means the poles.
- An energy expression written as `tr(h @ g1)` pairs indices wrongly; with
  `gamma_1[p,q] = <a+_p a_q>` it is `sum h[p,q] gamma_1[p,q]`. The wrong form showed TDHF
  failing to conserve its own energy functional, which is impossible and flagged the error.

## 8. Open items

- The trace-correction particular solution is a free choice; `blocked` destroys Mueller's
  `0.86` pole and a different choice may not.
- All reconstructions are far smaller than the true cumulant (`|lambda_2 exact| = 0.076` on H2
  against `|blocked| = 0.005`, `|mp2| = 0.049`). They are perturbative nudges.
- H3 is strongly correlated (`|gamma_1_HF - gamma_1_FCI| = 0.738`), so its pole comparison is
  reference-dominated.
- The raw `mp2` near-zero mode is not converged; the finite-difference Jacobian is unstable
  against small energy denominators. It vanishes under trace correction.
- PySCF 2.14 ships `pyscf.adc` with `radc_ee` / `uadc_ee` and `adc(2)`, `adc(2)-x`, `adc(3)`,
  available as a benchmark whenever wanted.

## 9. Reproducing

    python3 -m qdiss.eom.validate          derivations against FullCI
    python3 tools/physical_hamiltonian.py  PySCF Hamiltonians against FCI
    python3 tests/test_closures.py         antisymmetry, trace identity, TDHF equivalence
    python3 tools/poles_real.py            poles on H2 and H3 against exact
    python3 tools/trace_corrected.py       poles with and without trace correction
    python3 tools/stationarity.py          how stationary the FCI density is
    python3 tools/cumulant_phase.py        the in-phase / out-of-phase occupation result
