"""Is the reference state we perturb about even N-representable?

Screens HF, MP2 and FCI 2-RDMs of small electronic systems against the DQG
conditions, in the natural-geminal (spectral) basis.  Motivation: EE2's
CorrelateInitial() offers an MP2-natural-orbital initial condition; if that
reference is already outside the cone, no property of the generator can
rescue the trajectory.
"""
import sys, numpy as np
sys.path.insert(0, '/root/workspace')
from pyscf import gto, scf, mp, fci
from pyscf.scf import addons
from qdiss.representability.dqg import spinorb_2rdm, dqg_report, contract_to_1rdm

np.set_printoptions(linewidth=160, precision=6, suppress=True)


def hf_2rdm_so(mf, nso, nocc_a, nocc_b):
    """<a_i^+ a_j^+ a_l a_k> for a single determinant: g_ik g_jl - g_il g_jk."""
    g = np.zeros((nso, nso))
    n = nso // 2
    for i in range(nocc_a):
        g[i, i] = 1.0
    for i in range(nocc_b):
        g[n + i, n + i] = 1.0
    return (np.einsum('ik,jl->ijkl', g, g) - np.einsum('il,jk->ijkl', g, g))


def run(name, atom, basis='sto-3g'):
    print(f"\n=== {name}  ({basis}) ===")
    mol = gto.M(atom=atom, basis=basis, verbose=0)
    mf = scf.RHF(mol).run()
    norb, nelec = mf.mo_coeff.shape[1], mol.nelectron
    nso = 2 * norb
    print(f"  norb {norb}  nelec {nelec}  E(RHF) {mf.e_tot:.8f}")

    # --- HF
    D_hf = hf_2rdm_so(mf, nso, nelec // 2, nelec // 2)
    dqg_report(D_hf, nelec, "HF (single determinant)")

    # --- MP2 (UMP2 on RHF orbitals so we get spin blocks without re-converging)
    umf = addons.convert_to_uhf(mf)
    pt = mp.UMP2(umf).run()
    dm2aa, dm2ab, dm2bb = pt.make_rdm2()
    D_mp2 = spinorb_2rdm(dm2aa, dm2ab, dm2bb)
    print(f"  E(MP2) {pt.e_tot:.8f}")
    dqg_report(D_mp2, nelec, "MP2")

    # --- FCI
    cis = fci.FCI(mol, mf.mo_coeff)
    efci, ci = cis.kernel()
    (d1a, d1b), (d2aa, d2ab, d2bb) = fci.direct_spin1.make_rdm12s(
        ci, norb, (nelec // 2, nelec // 2))
    D_fci = spinorb_2rdm(d2aa, d2ab, d2bb)
    print(f"  E(FCI) {efci:.8f}")
    # sanity: trace and contraction consistency
    tr = np.einsum('ijij->', D_fci)
    g = contract_to_1rdm(D_fci, nelec)
    print(f"    sanity: sum_ij D2[i,j,i,j] = {tr:.6f}  (expect N(N-1) = {nelec*(nelec-1)})"
          f"   tr gamma = {np.trace(g):.6f}  (expect {nelec})")
    dqg_report(D_fci, nelec, "FCI (exact)")


if __name__ == '__main__':
    run("H4 linear, 1.0 A", '; '.join(f'H 0 0 {i*1.0}' for i in range(4)))
    run("H4 square, 1.23 A", 'H 0 0 0; H 1.23 0 0; H 0 1.23 0; H 1.23 1.23 0')
    run("H2O, 1.0x r_e (stretched)", 'O 0 0 0; H 0 -1.2 0.9; H 0 1.2 0.9')
