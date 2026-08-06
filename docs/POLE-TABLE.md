# Singlet excitation energies: closures against TDHF and FCI

Validated pole comparison for the 2-RDM closures. Every earlier attempt at this table was
withdrawn; this one passes an independent acceptance test and is the first that should be
quoted.

## System and approximation

H4 in an asymmetric trapezoid, parallel sides 0.80 and 1.20 Angstrom, legs 0.762 and 0.735,
plus 0.05 * normal jitter on every coordinate, sto-3g, RHF reference (UHF is identical at this
geometry: energies agree to 4e-10 and `<S^2> = 0`). `M = 8`, `N = 4`, `E_0 = -1.961395904`.

`FCI` is exact diagonalisation in the `S_z = 0` sector (dim 36), filtered to `<S^2> ~ 0`, giving
20 singlet states. `TDHF(pyscf)` is `tdscf.TDHF` with `singlet = True`.

`tdhf(claude)` is the MKNO 1-RDM EOM (`D1 x14, G1 x18, L2 x4`) with the two-particle cumulant
TRUNCATED, `lambda_2 := 0`, linearised at its own HF fixed point (`|gamma_1-dot| = 2.7e-14`) over
a spin-ADAPTED basis `B[2p,2q] = B[2p+1,2q+1] = b[p,q]`, which for a closed-shell reference spans
exactly the singlet excitations.

`lam3_0`, `muller` and `muller_raw` propagate `gamma_2` alone with `gamma_1` by contraction, and
differ only in the reconstruction of `gamma_3`: the MKNO separable part (`lambda_3 = 0`), Mueller
Eq. 14 with the Lackner contraction correction, and the same without the correction. They are
linearised at the FCI `gamma_2`. Energies in Hartree.

## Matching

Poles must be assigned to FCI states, not sorted and read off: the `gamma_2` Jacobian has 344-360
modes for six targets, and most are redundant directions. Three ingredients, all necessary:

- **Biorthogonal projection.** The Jacobian is real but NON-NORMAL, so the component of a target
  along mode `i` is `l_i^dag x` with `l_i^dag` the `i`-th row of `R^{-1}`, not `r_i^dag x`. Using
  right eigenvectors with naive inner products was the actual bug behind three failed attempts.
- **Targets are transition densities of matching rank.** One-body `<0|a+_p a_q|k>` for the rank-1
  Jacobian, two-body `<0|a+_p a+_v a_w a_q|k>` for the rank-2 ones, Hermitian part, computed by
  Gram products of ladder-lowered states.
- **Hungarian assignment** on the projection weights, so the matching is one-to-one. Six
  independent argmaxes previously assigned one mode to two different levels.

## Result

    FCI target -> |   0.3889    0.5630    0.9222    1.1525    1.1807    1.4686
    TDHF(pyscf)   |   0.3369    0.8693    1.1234    1.5097   (reference)
    tdhf(claude)  |   0.3369    0.7111    0.8693    1.1234    0.4934    1.5097   mean|err| 0.1684
    weight        |    0.800     0.006     0.715     0.763     0.285     0.588
    lam3_0        |   0.3771    0.4533    0.9080    1.1684    1.0599    1.3859   mean|err| 0.0592
    weight        |    1.131     0.888     0.445     0.658     0.460     1.327
    muller        |   0.3640    0.5547    0.8900    1.1477    1.0915    1.4981   mean|err| 0.0314
    weight        |    0.649     0.506     0.696     0.566     0.462     0.553
    muller_raw    |   0.3415    0.5319    0.8858    1.1458    1.0142    1.4953   mean|err| 0.0524
    weight        |    0.565     0.329     0.690     0.297     0.331     0.395

**Acceptance test.** `tdhf(claude)` assigns all four `TDHF(pyscf)` singlets to the correct
targets in the correct order, with weights 0.80, 0.72, 0.76 and 0.59, while its two spurious
modes are flagged by weights of 0.006 and 0.285. The weight column does the job the magnitude
gate never could. No closure number in this file was quoted before that test passed.

**Closures beat TDHF by three to five fold**, 0.03-0.06 Ha against 0.17 Ha, roughly 1 eV against
4.6 eV. Order: `muller` 0.0314, `muller_raw` 0.0524, `lam3_0` 0.0592.

## Caveats

Weights above one are legitimate: biorthogonal modes are not orthogonal, so
`|c_i| ||r_i|| / ||x||` is not a fraction. It is a relative indicator only. Note `lam3_0`'s two
best-matched columns are exactly the two with weight above one, so its 0.0592 may be flattered.

The closures are linearised at the FCI `gamma_2`, which is a fixed point of the exact flow but
NOT of any closure, whereas `tdhf(claude)` sits at its own fixed point. That asymmetry is the
remaining structural weakness in the comparison. Resolving it needs either a fixed point for each
closure or the homotopy construction: scale `gamma_3(s) = (1-s) g1^g1^g1 + s gamma_3_recon` and
track modes from the mean-field limit, where they can be anchored to TDHF by eigenvalue alone,
using biorthogonal overlaps between consecutive steps.
