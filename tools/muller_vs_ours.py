#!/usr/bin/env python3
"""Mueller's 3-RDM reconstruction as a baseline against ours, static then propagated.

Derivation of the EOM: qdiss.eom.gwt / qdiss.eom.wick (identical term lists), validated against
qdiss.eom.fullci at 1e-15. Propagated state: gamma_2 ALONE, with gamma_1 recovered by
contraction, which is the scheme Mueller et al. use (they propagate F_12). Matrix elements from
PySCF, checked against FCI. Truth: exact wavefunction propagation.

Reconstructions compared, all producing gamma_3 from gamma_1 and gamma_2:

    muller      qdiss.muller.gamma3, closed-system limit of their Eq. 14 plus the Eq. 44/45
                contraction correction
    muller_raw  the same without the contraction correction, to isolate what the correction buys
    lam3_0      ours: gamma_3 = the separable part of the MKNO cumulant expansion, i.e. lambda_3 = 0

Static checks first: exactness at N = 2, error against the exact gamma_3, and the contraction
residual |Tr_3(gamma_3)/(N-2) - gamma_2|. Then propagation, tracking the error in gamma_2, the
trace, positivity of gamma_2, and the Pauli bound on gamma_1.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from poles_real import build as build_so
from qdiss import muller
from qdiss.eom import fullci, wick
from qdiss.symbolic import mkno

UP3, LO3 = ("p", "v", "x"), ("q", "w", "y")
SEP3 = [(s, ",".join("".join(u) + "".join(l) for n, r, u, l in f) + "->pvxqwy",
         [(n, r) for n, r, u, l in f])
        for s, f in mkno.gamma_to_cumulants(UP3, LO3)
        if not any(n == "L" and r == 3 for n, r, u, l in f)]


def wedge(a, b):
    return np.einsum("pq,vw->pvqw", a, b) - np.einsum("pw,vq->pvqw", a, b)


def separable(g1, g2, M):
    l2 = g2 - wedge(g1, g1)
    t = {("G", 1): g1, ("L", 2): l2}
    out = np.zeros((M,) * 6, dtype=complex)
    for s, sub, keys in SEP3:
        out = out + s * np.einsum(sub, *[t[k] for k in keys], optimize=True)
    return out


def recon(name, g1, g2, N, M):
    if name == "muller":
        return muller.gamma3(g1, g2, N, contraction_correct=True)
    if name == "muller_raw":
        return muller.gamma3(g1, g2, N, contraction_correct=False)
    if name == "lam3_0":
        return separable(g1, g2, M)
    raise ValueError("unknown reconstruction %r; moment negation is excluded" % name)


# N >= 3 only: at N = 2 the 2-RDM is the whole state, the rank-2 EOM is exact, and the test is
# vacuous. Moment negation (gamma_3 = 0) is excluded by the project ground rules.
NAMES = ("muller", "muller_raw", "lam3_0")
SYS = (("H 0 0 0; H 0 0 1.6; H 0 0 3.2", "sto-3g", 3, 0, 1, "H3 lin / sto-3g  (N=3)"),
       ("H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0", "sto-3g", 4, 0, 0, "H4 lin / sto-3g  (N=4)"))

for atom, basis, N, chg, spin, label in SYS:
    M, h, V, enuc = build_so(atom, basis, chg, spin)
    Hm = fullci.build_H(h, V, M)
    dim = 1 << M
    popc = np.array([bin(s).count("1") for s in range(dim)])
    sub = np.where(popc == N)[0]
    es, vec = np.linalg.eigh(Hm[np.ix_(sub, sub)])
    psi0 = np.zeros(dim, dtype=complex)
    psi0[sub] = vec[:, 0]
    rr = np.random.default_rng(5)
    kick = np.zeros(dim, dtype=complex)
    kick[sub] = rr.normal(size=len(sub)) + 1j * rr.normal(size=len(sub))
    psi0 = psi0 + 0.3 * kick / np.linalg.norm(kick)
    psi0 /= np.linalg.norm(psi0)
    g1e = fullci.rdm(psi0, M, 1)
    g2e = fullci.rdm(psi0, M, 2)
    g3e = fullci.rdm(psi0, M, 3)
    print("=== %s   M=%d ===   |gamma_3 exact| = %.5f" % (label, M, np.abs(g3e).max()))
    print("   %-11s | %13s %13s %13s" % ("recon", "|g3-g3exact|", "contract res", "tr g3"))
    for nm in NAMES:
        g3 = recon(nm, g1e, g2e, N, M)
        cr = (np.abs(np.einsum("pvxqwx->pvqw", g3) / (N - 2) - g2e).max()
              if N > 2 else float("nan"))
        print("   %-11s | %13.5f %13.3e %13.4f"
              % (nm, np.abs(g3 - g3e).max(), cr, np.einsum("pvxpvx->", g3).real))

    EYE = np.eye(M, dtype=complex)
    w, U = np.linalg.eigh(Hm)
    PAIRS = [(p, v) for p in range(M) for v in range(p + 1, M)]

    def to_pair(T):
        A = np.array([[2.0 * T[p, v, q, w] for (q, w) in PAIRS] for (p, v) in PAIRS])
        return 0.5 * (A + A.conj().T)

    def g2dot(g2, nm):
        g1 = np.einsum("pxqx->pq", g2) / (N - 1)
        g3 = recon(nm, g1, g2, N, M)
        return wick.gamma_dot(2, h, V, {("D", 1): EYE, ("G", 1): g1,
                                       ("G", 2): g2, ("G", 3): g3}, M)

    def rk4(f, y, dt):
        k1 = f(y); k2 = f(y + dt / 2 * k1); k3 = f(y + dt / 2 * k2); k4 = f(y + dt * k3)
        return y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

    st = {nm: g2e.copy() for nm in NAMES}
    dt, ns, ev = 0.004, 250, 50
    print()
    print("   propagation: |gamma_2 - exact| / min eig D / n_min(gamma_1) / tr gamma_2")
    print("   %6s | %s" % ("t", " | ".join("%-28s" % nm for nm in NAMES)))
    for step in range(ns + 1):
        if step % ev == 0:
            t = step * dt
            psit = U @ (np.exp(-1j * w * t) * (U.conj().T @ psi0))
            ex = fullci.rdm(psit, M, 2)
            cells = []
            for nm in NAMES:
                g2 = st[nm]
                if not np.isfinite(g2).all():
                    cells.append("%-28s" % "diverged")
                    continue
                g1 = np.einsum("pxqx->pq", g2) / (N - 1)
                cells.append("%7.4f %7.4f %6.3f %5.2f"
                             % (np.abs(g2 - ex).max(),
                                np.linalg.eigvalsh(to_pair(g2)).min(),
                                np.linalg.eigvalsh(0.5 * (g1 + g1.conj().T)).min(),
                                np.einsum("pvpv->", g2).real))
            print("   %6.2f | %s" % (t, " | ".join(cells)))
        if step < ns:
            for nm in NAMES:
                if np.isfinite(st[nm]).all():
                    try:
                        st[nm] = rk4(lambda y, k=nm: g2dot(y, k), st[nm], dt)
                    except Exception:
                        st[nm] = np.full((M,) * 4, np.nan, dtype=complex)
    print()
