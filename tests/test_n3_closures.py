#!/usr/bin/env python3
"""N=3: does the closure choice decide whether the 2-RDM stays representable?

At N=3 gamma_3 is nonzero, so a closure does real work. Two endpoints:

    alpha = 0   gamma_3 = 0            exact at N=2, discards the 3-body term
    alpha = 1   lambda_3 = 0           gamma_3 = separable(gamma_1, lambda_2)

and the scan between them. John's hypothesis is that the cumulant
parameterization is what breaks positivity, because it hides the exact
cancellation that gamma_3 = 0 expresses at small N; if so, small alpha should
survive longer than alpha = 1.

The separable reconstruction is generated from mkno.gamma_to_cumulants -- the
same verified expansion that reproduced gamma_2 = lambda_2 + gamma_1 wedge
gamma_1 -- and is checked against explicit Fock space before being propagated.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from qdiss.symbolic import mkno
from qdiss.dynamics.coherent2 import compile_terms, evaluate, rk4
from qdiss.representability.dqg import contract_to_1rdm, build_Q, build_G
from tools.eom2_coherent import ONE, TWO
from gno import Fock

M, N = 6, 3
rng = np.random.default_rng(2026)
dim = 1 << M

# ---------- separable gamma_3, generated from the symbolic expansion ----------
UP3, LO3 = ("p", "v", "x"), ("q", "w", "y")
full3 = mkno.gamma_to_cumulants(UP3, LO3)
sep3 = [(s, f) for s, f in full3 if not any(n == "L" and r == 3 for n, r, u, l in f)]
print("gamma_3 expansion: %d terms total, %d separable (lambda_3 dropped)"
      % (len(full3), len(sep3)))


def compile_plain(terms, out):
    """Einsums over named tensors, no Hamiltonian operand."""
    comp = []
    for sign, facs in terms:
        subs, keys = [], []
        for name, rank, up, lo in facs:
            subs.append("".join(up) + "".join(lo))
            keys.append((name, rank))
        comp.append((sign, ",".join(subs) + "->" + out, keys))
    return comp


SEP3 = compile_plain(sep3, "pvxqwy")


def sep_gamma3(g1, l2):
    tens = {("G", 1): g1, ("L", 2): l2}
    out = np.zeros((M,) * 6, dtype=complex)
    for sign, sub, keys in SEP3:
        out = out + sign * np.einsum(sub, *[tens[k] for k in keys], optimize=True)
    return out


def lam2_of(g2, g1):
    return g2 - (np.einsum("pq,vw->pvqw", g1, g1)
                 - np.einsum("pw,vq->pvqw", g1, g1))


# ---------- validate the reconstruction against explicit Fock space ----------
psi0 = np.where(np.array([bin(s).count("1") for s in range(dim)]) == N,
                rng.normal(size=dim) + 1j * rng.normal(size=dim), 0.0)
psi0 = psi0 / np.linalg.norm(psi0)


def ad(p, v):
    out = np.zeros_like(v); bit = 1 << p
    for s in np.nonzero(v)[0]:
        if s & bit: continue
        out[s | bit] += ((-1.0) ** bin(s & (bit - 1)).count("1")) * v[s]
    return out


def aa(p, v):
    out = np.zeros_like(v); bit = 1 << p
    for s in np.nonzero(v)[0]:
        if not (s & bit): continue
        out[s ^ bit] += ((-1.0) ** bin(s & (bit - 1)).count("1")) * v[s]
    return out


def rdm2(psi):
    g = np.zeros((M,) * 4, dtype=complex)
    for p in range(M):
        for v in range(M):
            for q in range(M):
                for w in range(M):
                    g[p, v, q, w] = np.vdot(psi, ad(p, ad(v, aa(w, aa(q, psi)))))
    return g


def rdm3(psi):
    g = np.zeros((M,) * 6, dtype=complex)
    for p in range(M):
        for v in range(M):
            for x in range(M):
                for q in range(M):
                    for w in range(M):
                        for y in range(M):
                            g[p, v, x, q, w, y] = np.vdot(
                                psi, ad(p, ad(v, ad(x, aa(y, aa(w, aa(q, psi)))))))
    return g


g2_0 = rdm2(psi0)
g1_0 = contract_to_1rdm(g2_0, N)
l2_0 = lam2_of(g2_0, g1_0)
g3_true = rdm3(psi0)
g3_sep = sep_gamma3(g1_0, l2_0)
print("check: |gamma_3 true|      = %.4e" % np.abs(g3_true).max())
print("check: |separable part|    = %.4e" % np.abs(g3_sep).max())
print("check: |lambda_3| implied  = %.4e" % np.abs(g3_true - g3_sep).max())
print()

# ---------- Hamiltonian from PySCF-like random integrals (Hermitian) ----------
A = rng.normal(size=(M, M)) + 1j * rng.normal(size=(M, M))
h = (A + A.conj().T) / 2
B = rng.normal(size=(M,) * 4) + 1j * rng.normal(size=(M,) * 4)
V = B - np.transpose(B, (1, 0, 2, 3))
V = V - np.transpose(V, (0, 1, 3, 2))
V = (V + np.transpose(V.conj(), (2, 3, 0, 1))) / 2
h *= 0.5
V *= 0.5

ONE_C = compile_terms(ONE, "rs")
TWO_C = compile_terms(TWO, "rtsu")


def make_rhs(alpha):
    def rhs(g2):
        g1 = contract_to_1rdm(g2, N)
        l2 = lam2_of(g2, g1)
        g3 = alpha * sep_gamma3(g1, l2) if alpha else np.zeros((M,) * 6, dtype=complex)
        gam = {2: g2, 3: g3}
        out = evaluate(ONE_C, h, gam, M, out_rank=4)
        out = out + 0.25 * evaluate(TWO_C, V, gam, M, out_rank=4)
        return 1j * out
    return rhs


def Hpsi(v):
    out = np.zeros_like(v)
    for r in range(M):
        for s in range(M):
            out += h[r, s] * ad(r, aa(s, v))
    for r in range(M):
        for t in range(M):
            for s in range(M):
                for u in range(M):
                    c = V[r, t, s, u]
                    if c != 0:
                        out += 0.25 * c * ad(r, ad(t, aa(u, aa(s, v))))
    return out


def mineigs(g2):
    g1 = contract_to_1rdm(g2, N)
    res = {}
    for nm, T in (("D", g2), ("Q", build_Q(g2, g1)), ("G", build_G(g2, g1))):
        Mx = T.reshape(M * M, M * M)
        Mx = 0.5 * (Mx + Mx.conj().T)
        res[nm] = float(np.linalg.eigvalsh(Mx).min())
    return res


dt = 0.002
T = 1.0
nstep = int(round(T / dt))
alphas = [0.0, 0.5, 1.0]
states = {a: g2_0.copy() for a in alphas}
psi = psi0.copy()

print("N=%d M=%d  dt=%g  T=%g" % (N, M, dt, T))
print("%6s | %-22s | %-22s | %-22s" % ("t", "alpha=0 (gamma_3=0)",
                                       "alpha=0.5", "alpha=1 (lambda_3=0)"))
print("%6s | %7s %7s %6s | %7s %7s %6s | %7s %7s %6s"
      % ("", "err", "minD", "minQ", "err", "minD", "minQ", "err", "minD", "minQ"))

first_neg = {a: None for a in alphas}
for step in range(nstep + 1):
    if step % 100 == 0:
        gx = rdm2(psi)
        row = ["%6.2f |" % (step * dt)]
        for a in alphas:
            e = np.abs(states[a] - gx).max()
            m = mineigs(states[a])
            row.append(" %7.1e %+7.4f %+6.3f |" % (e, m["D"], m["Q"]))
        print("".join(row))
    for a in alphas:
        m = mineigs(states[a])
        if first_neg[a] is None and min(m["D"], m["Q"], m["G"]) < -1e-6:
            first_neg[a] = step * dt
        states[a] = rk4(make_rhs(a), states[a], dt)
    psi = rk4(lambda v: -1j * Hpsi(v), psi, dt)

print()
for a in alphas:
    print("  alpha=%-4s first DQG violation at t = %s"
          % (a, "%.3f" % first_neg[a] if first_neg[a] is not None else "never (T=%g)" % T))
