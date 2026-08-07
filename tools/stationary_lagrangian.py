#!/usr/bin/env python3
"""Stationary densities by constrained energy minimisation: residual iteration with DIIS.

Minimise the energy subject to constraints imposed by Lagrange multipliers, exactly as a
traditional HF solver does. Three constraint sets, each a separate Lagrangian:

    hf_var    L = E - Tr(Xi (g^2 - g)) - mu (tr g - N)
              The textbook HF variational problem: minimise the energy over idempotent densities
              of trace N. Idempotency is what a real HF solver's multiplier enforces (as
              orthonormality of the orbitals), and it implies 0 <= g <= 1 automatically. This is
              the ACCEPTANCE GATE: for lam2_0 it must return the SCF density and energy.

    eom       L = E - Tr(Lambda R) - mu (tr g - N)
              Minimise the energy subject to the equation of motion being stationary, R[g] = 0.
              No idempotency. This is the one that generalises to correlated closures, where the
              density must be allowed fractional occupations.

    eom_idem  both constraints at once, for cross-checking.

Gradients of every block come from torch autograd on the Lagrangian itself, so the residual is
[dL/dg, dL/dLambda, dL/dXi, dL/dmu] with consistent signs by construction, and the term arising
because the generator depends on gamma is exact. The iteration is x <- x - alpha r with
alpha = 1e-3 and DIIS on the residual history. DIIS needs no descent direction: it solves for the
combination minimising the residual norm, so it converges to stationary points of the residual
whether they are minima or saddles, which is what a Lagrangian requires.

APPROXIMATION per row. Solved variable is gamma_1 alone; the equation whose stationary points are
sought is the MKNO 1-RDM EOM, qdiss.eom.gwt.terms_cumulant(1), validated against exact Fock-space
algebra qdiss.eom.fullci at 1e-15. The energy is the same functional, E = h.g + (1/4) V.g2 with
g2 = g1 ^ g1 + lambda_2.
    lam2_0   two-particle cumulant TRUNCATED, lambda_2 := 0. This is TDHF, and the CONTROL: the
             answer is the HF density, known independently from PySCF SCF.
    blocked  ours: lambda_2 RECONSTRUCTED as A ^ f + f ^ A with A = gamma_1 - gamma_1^2 and f
             fixed by demanding the rank-2 partial-trace identity exactly. Written here as matrix
             functions of gamma_1 rather than in the natural-orbital basis, so autograd sees the
             full dependence and there is no degenerate-eigenvalue blowup.
Started from the FCI 1-RDM.
"""
import sys
import numpy as np
import torch
sys.path.insert(0, "."); sys.path.insert(0, "tools")
from molecular_2rdm import build_H_sector, molecule, sector
from qdiss import closures
from qdiss.eom import fullci, gwt

CUM1 = gwt.terms_cumulant(1)
BASE = np.array([[0.00, 0.00, 0.0], [0.80, 0.00, 0.0],
                 [1.05, 0.72, 0.0], [-0.15, 0.72, 0.0]])
geom = BASE + 0.05 * np.random.default_rng(20260805).normal(size=BASE.shape)
atom = "; ".join("H %.6f %.6f %.6f" % tuple(r) for r in geom)
mol, mf, M, h, V, enuc = molecule(atom, "sto-3g", 0, 0)
N = mol.nelectron
dt64 = torch.complex128
ht = torch.tensor(h, dtype=dt64); Vt = torch.tensor(V, dtype=dt64)
EYEt = torch.eye(M, dtype=dt64)

HB = []
for i in range(M):
    E = np.zeros((M, M), dtype=complex); E[i, i] = 1.0; HB.append(E)
for i in range(M):
    for j in range(i + 1, M):
        E = np.zeros((M, M), dtype=complex); E[i, j] = 1.0; E[j, i] = 1.0; HB.append(E)
        E = np.zeros((M, M), dtype=complex); E[i, j] = 1j; E[j, i] = -1j; HB.append(E)
HBt = torch.stack([torch.tensor(b, dtype=dt64) for b in HB])
nB = HBt.shape[0]
HBn = torch.stack([b / torch.sqrt(torch.real((b.conj() * b).sum())) for b in HBt])


def herm(x):
    return torch.tensordot(x.to(dt64), HBn, dims=([0], [0]))


def coeffs(Gm):
    return torch.real(torch.tensordot(HBn.conj(), Gm, dims=([1, 2], [0, 1])))


def wedge_t(a, b):
    return (torch.einsum("pq,vw->pvqw", a, b) - torch.einsum("pw,vq->pvqw", a, b))


def lam2_blocked_t(g1):
    """The blocked cumulant as matrix functions of gamma_1 -- no eigendecomposition.

    Identical to closures.coefficient('blocked') fed through _from_c: in the natural-orbital
    basis that construction is c[p,v] = -2(a_p f_v + f_p a_v) with a = n - n^2, which assembles
    to exactly sym_wedge(A, f). Both A and f are matrix functions of gamma_1, so the wedge form
    is basis-free, autograd-transparent, and safe under the degenerate occupations that H4 has.
    """
    A = g1 - g1 @ g1
    T = torch.einsum("pp->", A)
    den = T * EYEt - 2.0 * A
    Dinv = torch.linalg.inv(den)
    Q = torch.einsum("pp->", A @ Dinv)
    f = -(A @ Dinv) / (1.0 + Q)
    return wedge_t(A, f) + wedge_t(f, A)


TERMS = [(s, [("".join(u) + "".join(l)) for nn, r, u, l in f], [(nn, r) for nn, r, u, l in f])
         for s, f in list(CUM1[0])]
TERMS2 = [(s, [("".join(u) + "".join(l)) for nn, r, u, l in f], [(nn, r) for nn, r, u, l in f])
          for s, f in list(CUM1[1])]


def lam2_t(g1, closure):
    if closure == "lam2_0":
        return torch.zeros((M,) * 4, dtype=dt64)
    return lam2_blocked_t(g1)


def rhs_t(g1, closure):
    tab = {("D", 1): EYEt, ("G", 1): g1, ("L", 2): lam2_t(g1, closure)}
    tot = torch.zeros((M, M), dtype=dt64)
    for s, subs, keys in TERMS:
        tot = tot + s * torch.einsum(",".join(["rs"] + subs) + "->pq",
                                     ht, *[tab[k] for k in keys])
    for s, subs, keys in TERMS2:
        tot = tot + s * 0.25 * torch.einsum(",".join(["rtsu"] + subs) + "->pq",
                                            Vt, *[tab[k] for k in keys])
    return 1j * tot


def energy_t(g1, closure):
    g2 = wedge_t(g1, g1) + lam2_t(g1, closure)
    return torch.real(torch.einsum("pq,pq->", ht, g1)
                      + 0.25 * torch.einsum("pvqw,pvqw->", Vt, g2))


# --- references -------------------------------------------------------------------------------
# molecule() returns h and V already in the converged SCF molecular-orbital basis, so the HF
# density is exactly the projector onto the N lowest spin orbitals. No SCF loop of our own, and
# the energy functional is checked against PySCF's total energy as a convention test.
g_hf = np.diag(np.array([1.0] * N + [0.0] * (M - N), dtype=complex))
E_HF_ours = float(energy_t(torch.tensor(g_hf, dtype=dt64), "lam2_0")) + enuc
Hs = build_H_sector(h, V, M, sector(M, N))
es_, vecs = np.linalg.eigh(Hs)
psi0 = np.zeros(1 << M, dtype=complex); psi0[sector(M, N)] = vecs[:, 0]
g_fci = fullci.rdm(psi0, M, 1)
E_HF, E_FCI = float(mf.e_tot), float(es_[0].real) + enuc
print("reference: E_HF(pyscf) = %.9f   E_HF(our functional at g_HF) = %.9f   diff %.2e"
      % (E_HF, E_HF_ours, abs(E_HF - E_HF_ours)))
print("           E_FCI = %.9f   |R[g_HF]| (lam2_0) = %.3e"
      % (E_FCI, float(coeffs(rhs_t(torch.tensor(g_hf, dtype=dt64), "lam2_0")).norm())))
l2b = closures.lambda2(g_fci, "blocked")
print("           blocked cumulant check |torch - numpy| = %.3e"
      % np.abs(lam2_blocked_t(torch.tensor(g_fci, dtype=dt64)).numpy() - l2b).max())
print()

MODES = {"hf_var": (False, True), "eom": (True, False), "eom_idem": (True, True)}
NITER = 4000


def solve(closure, mode, x0):
    use_R, use_I = MODES[mode]
    xg = x0.clone().requires_grad_(True)
    xL = torch.zeros(nB, dtype=torch.float64, requires_grad=True)
    xX = torch.zeros(nB, dtype=torch.float64, requires_grad=True)
    mu = torch.zeros(1, dtype=torch.float64, requires_grad=True)
    var = [xg] + ([xL] if use_R else []) + ([xX] if use_I else []) + [mu]
    alpha = 1.0e-3
    hist_x, hist_r = [], []
    print("=== closure %-8s mode %-9s ===" % (closure, mode))
    print("   %5s | %11s %11s %11s %11s | %13s"
          % ("iter", "|dL/dg|", "|R|", "|g^2-g|", "|tr-N|", "E - E_HF"))
    for it in range(NITER):
        G = herm(xg)
        Ls = energy_t(G, closure) - mu[0] * (torch.real(torch.einsum("pp->", G)) - N)
        if use_R:
            Ls = Ls - torch.real((herm(xL).conj() * rhs_t(G, closure)).sum())
        if use_I:
            Ls = Ls - torch.real((herm(xX).conj() * (G @ G - G)).sum())
        grads = torch.autograd.grad(Ls, var)
        r = torch.cat([g.detach().reshape(-1) for g in grads])
        sgn = torch.cat([(torch.ones if i == 0 else (lambda n: -torch.ones(n)))(g.numel())
                         if i == 0 else -torch.ones(g.numel(), dtype=torch.float64)
                         for i, g in enumerate(grads)])
        x = torch.cat([v.detach().reshape(-1) for v in var])
        hist_x.append(x.clone()); hist_r.append(r.clone()); hist_d = None
        if len(hist_x) > 10:
            hist_x.pop(0); hist_r.pop(0)
        nh = len(hist_x)
        cand = x - alpha * (sgn * r)
        if nh >= 2:
            B = torch.zeros(nh + 1, nh + 1, dtype=torch.float64)
            for a in range(nh):
                for b in range(nh):
                    B[a, b] = torch.dot(hist_r[a], hist_r[b])
            B[nh, :nh] = 1.0; B[:nh, nh] = 1.0
            rv = torch.zeros(nh + 1, dtype=torch.float64); rv[nh] = 1.0
            try:
                c = torch.linalg.solve(B + 1e-10 * torch.eye(nh + 1, dtype=torch.float64), rv)
                cand = (sum(c[a] * hist_x[a] for a in range(nh))
                        - alpha * (sgn * sum(c[a] * hist_r[a] for a in range(nh))))
            except Exception:
                pass
        with torch.no_grad():
            o = 0
            for v in var:
                v.copy_(cand[o:o + v.numel()].reshape(v.shape)); o += v.numel()
        if it % 1000 == 0 or it == NITER - 1:
            with torch.no_grad():
                Gn = herm(xg)
                print("   %5d | %11.3e %11.3e %11.3e %11.3e | %13.2e"
                      % (it, grads[0].norm(),
                         coeffs(rhs_t(Gn, closure)).norm(),
                         float(torch.abs(Gn @ Gn - Gn).max()),
                         abs(float(torch.real(torch.einsum("pp->", Gn))) - N),
                         float(energy_t(Gn, closure)) + enuc - E_HF))
    with torch.no_grad():
        Gn = herm(xg).numpy()
        occ = np.sort(np.linalg.eigvalsh(0.5 * (Gn + Gn.conj().T)))[::-1]
        print("   occupations: %s" % np.array2string(occ, precision=4, suppress_small=True))
        print("   |gamma_1 - gamma_HF| = %.4e     |gamma_1 - gamma_FCI| = %.4e"
              % (np.abs(Gn - g_hf).max(), np.abs(Gn - g_fci).max()))
    print()
    return coeffs(herm(xg).detach())


if __name__ == "__main__":
    x_fci = coeffs(torch.tensor(g_fci, dtype=dt64))
    # Is HF a STABLE fixed point of the iteration? Start from gamma_HF kicked by a small random
    # Hermitian perturbation. If the iteration returns to gamma_HF the formulation and solver are
    # sound and the earlier failure was basin selection (aufbau) rather than a defect.
    rk = np.random.default_rng(11)
    Kk = rk.normal(size=(M, M)) + 1j * rk.normal(size=(M, M))
    Kk = 0.5 * (Kk + Kk.conj().T); Kk *= 0.02 / np.abs(Kk).max()
    x_hfk = coeffs(torch.tensor(g_hf + Kk, dtype=dt64))
    print("start from gamma_HF + kick:  |start - gamma_HF| = %.3e" % np.abs(Kk).max())
    print()
    for closure in ("lam2_0", "blocked"):
        for mode in ("hf_var", "eom_idem"):
            solve(closure, mode, x_hfk)
            solve(closure, mode, x_fci)
