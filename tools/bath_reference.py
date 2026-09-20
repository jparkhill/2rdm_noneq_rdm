#!/usr/bin/env python3
r"""Episode 10, step 1: the Born-Markov-secular (Davies/GKSL) bath dissipator on the full Fock space, as reference.

Following docs/BATH-DERIVATION.md. System H4 (M=8, N=4, S_z=0 sector... no: full N-sector, dim C(8,4)=70).
Coupling H_SB = A (x) (b + b^dag) summed over modes, A = x-dipole one-body operator; bath Ohmic,
J(w) = eta w exp(-w/w_c), thermal at kT. A1 Born, A2 Markov, A5 secular.

Under A5 with the coupling operator evolved under the EXACT H_S (this file; A3 is the closure we test later),
the GKSL form is
    D[rho] = sum_{w} gamma(w) ( A_w rho A_w^dag - 1/2 {A_w^dag A_w, rho} ),
    A_w = sum_{E_b - E_a = w} |a><a| A |b><b|      (jump lowering energy by w; w<0 raises),
    gamma(w) = 2 pi J(|w|) [ n_B(|w|) + theta(w) ]      detailed balance gamma(-w) = e^{-beta w} gamma(w).
Energy differences are binned with tolerance TOL (secular grouping). Fixed point: the Gibbs state of H_S in the
N-sector, by Davies' theorem -- CHECK 0 verifies this numerically (|D[rho_Gibbs]| ~ 0) and that the Lamb-shift is
omitted (it commutes with H_S, so it does not move the fixed point; noted, not included).

Outputs: the generator as a superoperator on the 70x70 sector density matrix; propagation from the HF
determinant to the fixed point; entropy production; and the exact RDM rates d/dt gamma_1, gamma_2 under D alone,
saved for CHECK 1-4 (rank-1 and rank-2 Wick contractions).
"""
import sys, time
import numpy as np, scipy.linalg as sla
sys.path.insert(0, "."); sys.path.insert(0, "tools")
import commutator_structure as cs
from molecular_2rdm import build_H_sector, molecule, sector
from qdiss.eom import fullci

ETA, WC, TOL = 0.05, 1.0, 1e-8


def J(w):
    return ETA * w * np.exp(-w / WC)


def gamma_rate(w, kT):
    aw = abs(w)
    if aw < 1e-12:
        # w -> 0 limit of 2 pi J(w) n_B(w) = 2 pi eta kT  (Ohmic); pure dephasing rate
        return 2 * np.pi * ETA * kT
    nB = 1.0 / (np.exp(aw / kT) - 1.0)
    return 2 * np.pi * J(aw) * (nB + (1.0 if w > 0 else 0.0))


class DaviesBath:
    def __init__(self, Hs, A_list, kT):
        """A_list: list of Hermitian one-body coupling operators (sector matrices), each with its own bath."""
        es, U = np.linalg.eigh(Hs)
        self.es, self.U, self.kT = es, U, kT
        n = len(es)
        self.jumps = []
        for A in A_list:
            self._add(A, es, U, n, kT)
        self._build(n)

    def _add(self, A, es, U, n, kT):
        Ae = U.conj().T @ A @ U                  # coupling in the energy basis
        n = len(es)
        gaps = es[None, :] - es[:, None]          # gaps[a,b] = E_b - E_a : |a><a| A |b><b| lowers energy by gaps[a,b]
        # secular grouping: unique gaps within TOL
        flat = gaps.reshape(-1)
        order = np.argsort(flat)
        groups = []; cur = [order[0]]
        for i in order[1:]:
            if abs(flat[i] - flat[cur[0]]) < TOL: cur.append(i)
            else: groups.append(cur); cur = [i]
        groups.append(cur)
        for g in groups:
            w = float(np.mean(flat[g]))
            Aw = np.zeros((n, n), dtype=complex)
            for i in g:
                a, b = divmod(i, n); Aw[a, b] = Ae[a, b]
            if np.abs(Aw).max() < 1e-14: continue
            self.jumps.append((gamma_rate(w, kT), Aw, w))

    def _build(self, n):
        self.n = n
        # vectorised form: sum_w g_w A_w rho A_w^dag = sum_{ab,cd} R[a,b,c,d] rho[b,d] with
        # R[a,b,c,d] = sum_w g_w Aw[a,b] conj(Aw[c,d]); since each (a,b) belongs to exactly one group,
        # R[a,b,c,d] = G[a,b,c,d] Ae[a,b] conj(Ae[c,d]) where G = g(gap) if gap[a,b]==gap[c,d] else 0.
        # Build the sparse structure: for each group the outer product of its entries. n=70 -> n^4 = 24M complex = 384MB; do it as
        # a (n^2 x n^2) matrix K[(a,b),(c,d)] acting on vec(rho): rho'[a,c] = sum_bd K[(ab),(cd)] rho[b,d].
        K = np.zeros((n * n, n * n), dtype=complex)
        Gsum = np.zeros((n, n), dtype=complex)                  # sum_w g_w A_w^dag A_w
        for g, Aw, w in self.jumps:
            v = Aw.reshape(-1)
            nz = np.nonzero(v)[0]
            K[np.ix_(nz, nz)] += g * np.outer(v[nz], v[nz].conj())
            Gsum += g * (Aw.conj().T @ Aw)
        self.K, self.Gsum = K, Gsum

    def D(self, rho_e):
        """Dissipator on rho in the ENERGY basis (vectorised)."""
        n = self.n
        # sum_w g A_w rho A_w^dag : out[a,c] = sum_{b,d} K[(a,b),(c,d)] rho[b,d]
        Kt = self.K.reshape(n, n, n, n)
        jump = np.einsum("abcd,bd->ac", Kt, rho_e)
        return jump - 0.5 * (self.Gsum @ rho_e + rho_e @ self.Gsum)

    def rhs(self, rho_e):
        return -1j * (np.diag(self.es) @ rho_e - rho_e @ np.diag(self.es)) + self.D(rho_e)


def main():
    kT = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
    atom, N, spin = cs.H4
    mol, mf, M, h, V, enuc = molecule(atom, "sto-3g", 0, spin); dets = sector(M, N)
    Hs = build_H_sector(h, V, M, dets)
    # coupling operators: the x-dipole plus a generic set of Hermitian one-body operators so that no symmetry
    # sector is left uncoupled (a single dipole has selection rules and relaxes to conditional Gibbs states)
    rng = np.random.default_rng(3)
    A_list = [cs.sector_operator(cs.dipole(mol, mf, M), M, dets)]
    for _ in range(int(sys.argv[3]) if len(sys.argv) > 3 else 3):
        B = rng.normal(size=(M, M)); B = 0.5 * (B + B.T) / np.sqrt(M)
        A_list.append(cs.sector_operator(B.astype(complex), M, dets))
    bath = DaviesBath(Hs, A_list, kT)
    es, U = bath.es, bath.U
    print("H4 M=%d N=%d sector dim %d | Ohmic eta=%.3g w_c=%.3g | kT=%.2f | %d coupling operators, %d secular jump groups" % (M, N, len(dets), ETA, WC, kT, len(A_list), len(bath.jumps)))
    # CHECK 0: Gibbs state is a fixed point; CP structure (rates >= 0)
    w = np.exp(-(es - es[0]) / kT); w /= w.sum(); rhoG = np.diag(w).astype(complex)
    print("CHECK 0  |D[rho_Gibbs]| = %.2e   min rate = %.2e   detailed balance max|g(-w)/g(w) - e^{-w/kT}| = %.1e"
          % (np.abs(bath.D(rhoG)).max(), min(g for g, _, _ in bath.jumps),
             max(abs(gamma_rate(-wv, kT) / gamma_rate(wv, kT) - np.exp(-wv / kT)) for _, _, wv in bath.jumps if wv > 1e-6)))
    # propagate from the HF determinant (in the sector basis -> energy basis)
    hf = sum(1 << p for p in range(N)); i0 = dets.index(hf)
    rho = np.zeros((len(dets),) * 2, dtype=complex); rho[i0, i0] = 1.0
    rho_e = U.conj().T @ rho @ U
    dt, T = 0.05, float(sys.argv[2]) if len(sys.argv) > 2 else 400.0
    nsteps = int(T / dt)
    print("propagating rho from the HF determinant, RK4 dt=%.2f a.u., t_end=%.0f a.u. = %.1f fs" % (dt, T, T * 0.02419))
    print("    %7s | %8s %8s %10s %10s | %8s | %8s" % ("t", "E", "S", "|rho-Gibbs|", "min eig rho", "tr", "|D|"))
    def S_of(r):
        x = np.linalg.eigvalsh(0.5 * (r + r.conj().T)); x = x[x > 1e-14]; return -np.sum(x * np.log(x))
    def rk4(f, y, dt):
        k1 = f(y); k2 = f(y + .5 * dt * k1); k3 = f(y + .5 * dt * k2); k4 = f(y + dt * k3); return y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    rec = []
    for k in range(nsteps + 1):
        if k % int(20 / dt) == 0 or k == nsteps:
            E = np.real(np.trace(np.diag(es) @ rho_e)); S = S_of(rho_e)
            rec.append((k * dt, E, S, np.linalg.norm(rho_e - rhoG), np.linalg.eigvalsh(0.5 * (rho_e + rho_e.conj().T)).min()))
            print("    %7.1f | %8.5f %8.5f %10.2e %+10.1e | %8.6f | %.1e" % (rec[-1] + (np.real(np.trace(rho_e)), np.linalg.norm(bath.D(rho_e)))))
        if k < nsteps: rho_e = rk4(bath.rhs, rho_e, dt)
    EG = float(w @ es); SG = -np.sum(w * np.log(w))
    print("Gibbs: E=%.5f S=%.5f" % (EG, SG))
    # save: the generator's action on RDMs at a few correlated states for CHECK 1-4
    np.save("results/bath_ref_kT%.2f.npy" % kT, dict(es=es, U=U, jumps=[(g, Aw, wv) for g, Aw, wv in bath.jumps], rec=rec, dets=dets), allow_pickle=True)


if __name__ == "__main__":
    main()
