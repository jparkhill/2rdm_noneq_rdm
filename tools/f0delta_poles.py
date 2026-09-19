#!/usr/bin/env python3
r"""Route F, episode 7, step 4: poles of the F0 + delta closure against FCI.

Same system and recipe as docs/POLE-TABLE.md / tools/singlet_poles.py: H4 asymmetric trapezoid with the
0.05-jitter geometry (seed 20260805), sto-3g, M = 8, N = 4; FCI in the S_z = 0 sector filtered to <S^2> ~ 0;
targets = the six lowest singlet excitation energies; the closure is linearised at the FCI ground-state gamma_2
by finite differences over a real basis of Hermitian pair-space perturbations; modes are assigned to targets by
biorthogonal projection onto the two-body transition densities <0|a+_p a+_v a_w a_q|k> (Hermitian part) and a
one-to-one Hungarian assignment. Controls: `f0` (must reproduce POLE-TABLE lam3_0, 0.0592) and `muller`
(0.0314). New rows: `f0delta_D`, `f0delta`. Also reported: the largest real part among Jacobian eigenvalues
(a non-Hamiltonian closure may have damped/growing modes).

MIXED-STATE ROW. The same procedure at the thermal state kT = 0.5 over the S_z = 0 sector: targets are the
transition frequencies E_k - E_j between populated states, weighted by |w_j - w_k|; the six largest weights
(by |w_j - w_k| times the norm of the transition density) are matched.
"""
import sys
import time
import numpy as np
from scipy.optimize import linear_sum_assignment
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
import commutator_structure as cs
from commutator_k3 import from_mat2
from f0delta import F0Delta
from molecular_2rdm import molecule
from qdiss.eom import fullci

BASE = np.array([[0.00, 0.00, 0.0], [0.80, 0.00, 0.0], [1.05, 0.72, 0.0], [-0.15, 0.72, 0.0]])
geom = BASE + 0.05 * np.random.default_rng(20260805).normal(size=BASE.shape)
atom = "; ".join("H %.6f %.6f %.6f" % tuple(r) for r in geom)
mol, mf, M, h, V, enuc = molecule(atom, "sto-3g", 0, 0)
N = mol.nelectron
n = M // 2
EPS = 1e-6

sz0 = [d for d in range(1 << M)
       if bin(d & 0x5555555555).count("1") == N // 2 and bin(d & 0xAAAAAAAAAA).count("1") == N // 2]
idx = {d: i for i, d in enumerate(sz0)}
nd = len(sz0)


def op_matrix(apply_fn):
    Om = np.zeros((nd, nd), dtype=complex)
    b = np.zeros(1 << M, dtype=complex)
    for j, d in enumerate(sz0):
        b[:] = 0; b[d] = 1.0
        col = apply_fn(b)
        for d2, i in idx.items():
            Om[i, j] = col[d2]
    return Om


def Hfn(b):
    col = np.zeros_like(b)
    for p in range(M):
        for q in range(M):
            if abs(h[p, q]) > 1e-14:
                col += h[p, q] * fullci.ad(p, fullci.aa(q, b))
    for p, t, q, u in np.argwhere(np.abs(V) > 1e-14):
        col += 0.25 * V[p, t, q, u] * fullci.ad(p, fullci.ad(t, fullci.aa(u, fullci.aa(q, b))))
    return col


def S2fn(b):
    sp = np.zeros_like(b)
    for p in range(n):
        sp += fullci.ad(2 * p, fullci.aa(2 * p + 1, b))
    out = np.zeros_like(b)
    for p in range(n):
        out += fullci.ad(2 * p + 1, fullci.aa(2 * p, sp))
    return out


def embed(c):
    f = np.zeros(1 << M, dtype=complex); f[sz0] = c; return f


def trans2(ca, cb):
    """<a| a+_p a+_v a_w a_q |b> as a 4-index tensor, via the ladder of pair annihilators."""
    La = fullci._ladder(embed(ca), M, 2)
    Lb = fullci._ladder(embed(cb), M, 2)
    return (La.conj() @ Lb.T).reshape(M, M, M, M)


PAIRS = cs.pairs(M)
NP = len(PAIRS)
B2 = []
for a in range(NP):
    for b in range(a, NP):
        for ph in (1.0, 1j):
            Z = np.zeros((NP, NP), dtype=complex); Z[a, b] += ph; Z[b, a] += np.conj(ph)
            if np.abs(Z).max() < 1e-12:
                continue
            B2.append(from_mat2(Z, M))
n2 = len(B2)
nrm = np.array([float(np.real(np.tensordot(b.conj(), b, axes=4))) for b in B2])


def jacobian(cl, g2_0, mode):
    f0 = cl.g2dot(g2_0, mode)
    J = np.zeros((n2, n2))
    t0 = time.time()
    for k in range(n2):
        dn = (cl.g2dot(g2_0 + EPS * B2[k], mode) - f0) / EPS
        J[:, k] = np.array([float(np.real(np.tensordot(B2[m].conj(), dn, axes=4))) for m in range(n2)]) / nrm
    return J, time.time() - t0


def assign_bio(J, targetsH, targets, label):
    w, R = np.linalg.eig(J)
    Rinv = np.linalg.inv(R)
    keep = np.where(w.imag > 1e-6)[0]
    rn = np.array([np.linalg.norm(R[:, i]) for i in range(R.shape[1])])
    S = np.zeros((len(targets), len(keep)))
    for a in range(len(targets)):
        Hm = targetsH[a]
        x = np.array([float(np.real(np.tensordot(B2[m].conj(), Hm, axes=4))) for m in range(n2)]) / nrm
        nx = np.linalg.norm(x)
        c = Rinv @ x
        for b, i in enumerate(keep):
            S[a, b] = abs(c[i]) * rn[i] / (nx + 1e-30)
    r, cc = linear_sum_assignment(-S)
    row = np.full(len(targets), np.nan); sc = np.zeros(len(targets)); damp = np.zeros(len(targets))
    for a, b in zip(r, cc):
        row[a] = w[keep[b]].imag; sc[a] = S[a, b]; damp[a] = w[keep[b]].real
    err = np.nanmean(np.abs(row - targets))
    print("%-12s | %s  mean|err| %.4f   (%d modes; max|Re w| over all modes %.1e)"
          % (label, "  ".join("%8.4f" % x for x in row), err, len(keep), np.abs(w.real).max()))
    print("%-12s | %s  <- weight" % ("", "  ".join("%8.3f" % x for x in sc)))
    print("%-12s | %s  <- Re(w) of the assigned mode" % ("", "  ".join("%+8.1e" % x for x in damp)))
    return err


def main():
    Hs = op_matrix(Hfn); Hs = 0.5 * (Hs + Hs.conj().T)
    S2s = op_matrix(S2fn); S2s = 0.5 * (S2s + S2s.conj().T)
    es, vec = np.linalg.eigh(Hs)
    s2 = np.real(np.einsum("ik,ij,jk->k", vec.conj(), S2s, vec))
    sing = np.where(np.abs(s2) < 1e-6)[0]
    print("S_z=0 sector dim %d; singlets %d; E0 = %.9f" % (nd, len(sing), es[sing[0]] + enuc))
    exc = (es - es[sing[0]]).real
    sing_exc = np.array([exc[k] for k in sing[1:]])
    keep = np.argsort(sing_exc)[:6]
    targets = sing_exc[keep]
    tstate = [sing[1:][k] for k in keep]
    print("FCI singlet excitations (Ha): %s" % np.array2string(targets, precision=4))
    H2 = []
    for k in tstate:
        T = trans2(vec[:, sing[0]], vec[:, k])
        H2.append(0.5 * (T + np.transpose(T.conj(), (2, 3, 0, 1))))
    cl = F0Delta(h, V, M, N)
    g2_0 = fullci.rdm(embed(vec[:, sing[0]]), M, 2)
    print("\n### GROUND STATE   FCI target -> | %s" % "  ".join("%8.4f" % x for x in targets))
    for mode in ("f0", "muller", "f0delta", "f0delta_hard", "f0delta_hardonly"):
        J, dt = jacobian(cl, g2_0, mode)
        assign_bio(J, H2, targets, mode)
        print("%-12s | (Jacobian %.0fs)" % ("", dt))

    # ---------------- mixed-state row: thermal kT = 0.5 over the S_z = 0 sector
    kT = 0.5
    w = np.exp(-(es - es[0]) / kT); w /= w.sum()
    pop = np.where(w > 1e-3)[0]
    g2_T = sum(w[j] * fullci.rdm(embed(vec[:, j]), M, 2) for j in pop)
    cands = []
    for j in pop:
        for k in range(nd):
            if es[k] - es[j] > 1e-6 and abs(w[j] - w[k]) > 1e-4:
                T = trans2(vec[:, j], vec[:, k])
                Hm = 0.5 * (T + np.transpose(T.conj(), (2, 3, 0, 1)))
                cands.append((abs(w[j] - w[k]) * np.linalg.norm(Hm), es[k] - es[j], Hm, j, k))
    cands.sort(key=lambda t: -t[0])
    # distinct frequencies, six largest weights
    chosen, seen = [], []
    for c in cands:
        if all(abs(c[1] - s) > 1e-6 for s in seen):
            chosen.append(c); seen.append(c[1])
        if len(chosen) == 6:
            break
    tT = np.array([c[1] for c in chosen]); HT = [c[2] for c in chosen]
    print("\n### THERMAL kT=%.1f (populated states %d)   targets E_k - E_j -> | %s" % (kT, len(pop), "  ".join("%8.4f" % x for x in tT)))
    print("%-12s | %s  <- |w_j - w_k| |T|" % ("", "  ".join("%8.3f" % c[0] for c in chosen)))
    for mode in ("f0", "muller", "f0delta", "f0delta_hard", "f0delta_hardonly"):
        J, dt = jacobian(cl, g2_T, mode)
        assign_bio(J, HT, tT, mode)
        print("%-12s | (Jacobian %.0fs)" % ("", dt))


if __name__ == "__main__":
    main()
