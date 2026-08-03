#!/usr/bin/env python3
"""Does the 1-body EOM with lambda_2 = 0 reproduce TDHF? Checked against PySCF.

The 1-body EOM is  i d<a^p_q>/dt = <[H, a^p_q]>, whose two commutators were
derived and verified in tools/eom1.py:

    <[a^r_s,   a^p_q]> = d^p_s g^r_q - d^r_q g^p_s
    <[a^rt_su, a^p_q]> = d^p_s g2^rt_qu - d^p_u g2^rt_qs
                        - d^r_q g2^pt_su + d^t_q g2^pr_su

Closing with lambda_2 = 0, i.e. g2^{pr}_{qs} = g^p_q g^r_s - g^p_s g^r_q, should
give exactly TDHF. Rather than propagate and fit peaks out of an FFT, linearize
the closed EOM about the HF ground state and diagonalize: the eigenvalues must be
+/- i omega_n with omega_n the RPA/TDHF excitation energies PySCF reports.

Contractions are built from the symbolic term lists, not transcribed.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from pyscf import gto, scf, tdscf, ao2mo
from qdiss.dynamics.coherent2 import compile_terms, evaluate
from tools.eom1 import one as ONE1, two as TWO1

np.set_printoptions(precision=6, suppress=True)

# ---------- a minimal molecule ----------
mol = gto.M(atom="H 0 0 0; H 0 0 0.74", basis="6-31g", verbose=0)
mf = scf.RHF(mol).run()
print("molecule : H2 / 6-31G   E(RHF) = %.10f" % mf.e_tot)

td = tdscf.TDHF(mf).run(nstates=8)
tdt = tdscf.TDHF(mf); tdt.singlet=False; tdt.run(nstates=8)
w_s = np.array(sorted(w for w in td.e if w > 1e-8))
w_t = np.array(sorted(w for w in tdt.e if w > 1e-8))
w_ref = np.array(sorted(set(np.round(np.concatenate([w_s,w_t]),10))))
print("PySCF TDHF singlets:", w_s)
print("PySCF TDHF triplets:", w_t)
print()

# ---------- integrals in the MO spin-orbital basis ----------
nao = mf.mo_coeff.shape[1]
hcore_mo = mf.mo_coeff.T @ mf.get_hcore() @ mf.mo_coeff
eri_mo = ao2mo.restore(1, ao2mo.kernel(mol, mf.mo_coeff), nao)   # (pq|rs)

M = 2 * nao                       # spin orbitals: even = alpha, odd = beta
nocc = mol.nelectron


def spin_h():
    h = np.zeros((M, M), dtype=complex)
    for p in range(nao):
        for q in range(nao):
            for s in range(2):
                h[2 * p + s, 2 * q + s] = hcore_mo[p, q]
    return h


def spin_V():
    """V[r,t,s,u] for H2 = (1/2) sum V a^{rt}_{su},  a^{rt}_{su}=a+_r a+_t a_u a_s.

    <rt|su> in physicists' order is (rs|tu) in chemists' order, spin-diagonal on
    (r,s) and (t,u). Antisymmetrized as the operator already is.
    """
    V = np.zeros((M,) * 4, dtype=complex)
    for r in range(M):
        for t in range(M):
            for s in range(M):
                for u in range(M):
                    if (r % 2 == s % 2) and (t % 2 == u % 2):
                        V[r, t, s, u] += eri_mo[r // 2, s // 2, t // 2, u // 2]
                    if (r % 2 == u % 2) and (t % 2 == s % 2):
                        V[r, t, s, u] -= eri_mo[r // 2, u // 2, t // 2, s // 2]
    return V


h = spin_h()
V = spin_V()

# ---------- the closed 1-body EOM ----------
ONE_C = compile_terms(ONE1, "rs", out="pq")
TWO_C = compile_terms(TWO1, "rtsu", out="pq")


def g2_from_g1(g):
    """lambda_2 = 0 : g2[p,r,q,s] = g[p,q] g[r,s] - g[p,s] g[r,q]"""
    return (np.einsum("pq,rs->prqs", g, g) - np.einsum("ps,rq->prqs", g, g))


def rhs(g):
    gam = {1: g, 2: g2_from_g1(g)}
    out = evaluate(ONE_C, h, gam, M, out_rank=2)
    out = out + 0.25 * evaluate(TWO_C, V, gam, M, out_rank=2)
    return 1j * out


# ---------- HF reference density ----------
g0 = np.zeros((M, M), dtype=complex)
for i in range(nocc):
    g0[i, i] = 1.0
print("stationarity of the HF density : |rhs(g0)| = %.3e" % np.abs(rhs(g0)).max())

# ---------- linearize and diagonalize ----------
eps = 1e-6
n = M * M
J = np.zeros((n, n), dtype=complex)
for k in range(n):
    dg = np.zeros((M, M), dtype=complex)
    dg[k // M, k % M] = eps
    J[:, k] = ((rhs(g0 + dg) - rhs(g0 - dg)) / (2 * eps)).reshape(-1)

ev = np.linalg.eigvals(J)
w = np.array(sorted(abs(e.imag) for e in ev if abs(e.imag) > 1e-6))
uniq = []
for x in w:
    if not uniq or abs(x - uniq[-1]) > 1e-6:
        uniq.append(x)
print("EOM linear-response frequencies :", np.array(uniq))
print()

for wr in w_ref:
    close = [x for x in uniq if abs(x - wr) < 1e-5]
    print("  TDHF %.8f  ->  %s" % (wr, ("matched %.8f" % close[0]) if close else "NOT FOUND"))
