#!/usr/bin/env python3
"""The CP question at fixed rho reduces to one PSD condition on the kernel block.

tools/cp_closure.py ran an SDP over a jump basis truncated to the first 40 of 224
matrix units, in row-major order, and reported the exact flow sitting 6.4% outside the
GKSL cone. That number is not trustworthy: the subspace was biased, and the only run
that came out at zero residual was gamma_3 = 0, which is a pure commutator and needs no
jump operators at all, so it exercised nothing of the dissipator.

No SDP is needed. Write K for ker rho and R for its orthogonal complement. On the kernel,

    <v| [K, rho] |v>          = 0     since rho|v> = 0 and <v|rho = 0
    <v| {F_b^dag F_a, rho} |v> = 0     for the same reason
    <v| F_a rho F_b^dag |v>            survives, and summed against c >= 0 is PSD

so the kernel block of any GKSL generator's output is PSD, and the K-R and R-R blocks are
freely reachable through the Hamiltonian part. Hence

    target is in the GKSL cone image at rho   <==>   X_KK is PSD

which is precisely the Nagumo condition for forward invariance of the PSD cone. CP at a
point adds nothing to positivity preservation.

Here rank gamma_2 <= 6, because a_v a_p psi lies in the one-electron sector, so on the
15-dimensional pair space the kernel is 9-dimensional and the condition has content.
Predictions: exact gives X_KK PSD; gamma_3 = 0 gives X_KK identically zero, since an
isospectral flow maps the kernel to itself; lambda_3 = 0 gives negative eigenvalues,
which is the D violation already on record, now identified as the CP obstruction itself.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from qdiss.symbolic import mkno
from tools.eom2_coherent import ONE, TWO

M, N = 6, 3
dim = 1 << M
PAIRS = [(p, v) for p in range(M) for v in range(p + 1, M)]
NP = len(PAIRS)
UP3, LO3 = ("p", "v", "x"), ("q", "w", "y")
SEP3 = [(s, ",".join("".join(u) + "".join(l) for n, r, u, l in f) + "->pvxqwy",
         [(n, r) for n, r, u, l in f])
        for s, f in mkno.gamma_to_cumulants(UP3, LO3)
        if not any(n == "L" and r == 3 for n, r, u, l in f)]


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


def rdms(psi):
    Q1 = np.array([aa(p, psi) for p in range(M)])
    Q2 = np.array([aa(v, aa(p, psi)) for p in range(M) for v in range(M)])
    Q3 = np.array([aa(x, aa(v, aa(p, psi)))
                   for p in range(M) for v in range(M) for x in range(M)])
    return (Q1.conj() @ Q1.T, (Q2.conj() @ Q2.T).reshape(M, M, M, M),
            (Q3.conj() @ Q3.T).reshape((M,) * 6))


def to_pair(T):
    return np.array([[2.0 * T[p, v, q, w] for (q, w) in PAIRS] for (p, v) in PAIRS])


def sep(g1, g2):
    l2 = g2 - (np.einsum("pq,vw->pvqw", g1, g1) - np.einsum("pw,vq->pvqw", g1, g1))
    t = {("G", 1): g1, ("L", 2): l2}
    out = np.zeros((M,) * 6, dtype=complex)
    for s, sub, keys in SEP3:
        out = out + s * np.einsum(sub, *[t[k] for k in keys], optimize=True)
    return out


def g2dot(h, V, g1, g2, g3):
    eye = np.eye(M, dtype=complex)
    gam = {1: g1, 2: g2, 3: g3}
    tot = np.zeros((M,) * 4, dtype=complex)
    for tl, Hop, pref, hsub in ((ONE, h, 1.0, "rs"), (TWO, V, 0.25, "rtsu")):
        for sign, facs in tl:
            subs = [hsub]; ops = [Hop]
            for n, r, u, l in facs:
                subs.append("".join(u) + "".join(l))
                ops.append(eye if n == "D" else gam[r])
            tot = tot + sign * pref * np.einsum(",".join(subs) + "->pvqw", *ops,
                                                optimize=True)
    return 1j * tot


popc = np.array([bin(s).count("1") for s in range(dim)])
print("%4s %6s | %28s | %28s | %20s"
      % ("H", "dim K", "exact: min/max eig X_KK", "lambda_3=0", "gamma_3=0"))
agg = {"exact": [], "l3": [], "g3": []}
for tr in range(8):
    r = np.random.default_rng(11 + tr)
    A_ = r.normal(size=(M, M)) + 1j * r.normal(size=(M, M))
    h = 0.25 * (A_ + A_.conj().T)
    Bt = r.normal(size=(M,) * 4) + 1j * r.normal(size=(M,) * 4)
    V = Bt - np.transpose(Bt, (1, 0, 2, 3))
    V = V - np.transpose(V, (0, 1, 3, 2))
    V = 0.25 * (V + np.transpose(V.conj(), (2, 3, 0, 1)))
    psi = np.where(popc == N, r.normal(size=dim) + 1j * r.normal(size=dim), 0.0)
    psi = psi / np.linalg.norm(psi)
    g1, g2, g3 = rdms(psi)
    rho = 0.5 * (to_pair(g2) + to_pair(g2).conj().T)
    w, U = np.linalg.eigh(rho)
    Kb = U[:, w < 1e-9]
    cells = []
    for key, g3u in (("exact", g3), ("l3", sep(g1, g2)),
                     ("g3", np.zeros((M,) * 6, dtype=complex))):
        X = to_pair(g2dot(h, V, g1, g2, g3u))
        X = 0.5 * (X + X.conj().T)
        XKK = Kb.conj().T @ X @ Kb
        e = np.linalg.eigvalsh(0.5 * (XKK + XKK.conj().T))
        agg[key].append(e.min()); cells.append((e.min(), e.max()))
    print("%4d %6d | %13.3e %13.3e | %13.3e %13.3e | %9.2e %9.2e"
          % (tr, Kb.shape[1], cells[0][0], cells[0][1], cells[1][0], cells[1][1],
             cells[2][0], cells[2][1]))
print()
for k, nm in (("exact", "exact"), ("l3", "lambda_3=0"), ("g3", "gamma_3=0")):
    a = np.array(agg[k])
    print("   %-11s worst min eig of X_KK = %+0.4e   %d of %d violate PSD"
          % (nm, a.min(), int((a < -1e-9).sum()), len(a)))
print()
print("X_KK PSD is necessary and sufficient for the flow to lie in the GKSL cone at rho,")
print("and identical to the Nagumo condition. A nonzero SDP residual for the exact flow")
print("would contradict unitarity, so 6.4% was an artifact of the truncated jump basis.")
