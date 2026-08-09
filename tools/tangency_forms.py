#!/usr/bin/env python3
"""Two checks on the tangency programme at rank 2 -- noting rank 2 is one rung, not the goal.

THE MACHINERY BEING TESTED. For the positive semidefinite cone and a boundary point with kernel
projector P_K, the tangent cone is {Y : P_K Y P_K >= 0}, so Nagumo's no-outward-flow condition is
    P_K (d/dt of the matrix) P_K  >=  0.
Two structural forms satisfy it for free, whatever A and L_i are:
    A X + X A^dag          compresses to EXACTLY zero, since X P_K = P_K X = 0  -> tangent
    sum_i L_i X L_i^dag    compresses to sum_i (P_K L_i) X (P_K L_i)^dag >= 0    -> inward
Their sum is Lindblad form. The same statement holds at every rank k with the k-tuple matrix in
place of the pair matrix; nothing below is special to rank 2 except the arithmetic.

CHECK 1. Does the gamma_2-dependent part of gamma_2-dot already have the Lyapunov shape
A X + X A^dag? Sharp test: such a map sends a RANK-1 input |Vec><Vec| to
(A Vec)<Vec| + |Vec>(A Vec)^dag, whose rank is at most 2 and whose range lies in
span{Vec, A Vec}. So feed rank-1 pair matrices and measure the output rank. Terms are taken from
qdiss.eom.wick.terms(2), the bare-moment rank-2 EOM validated against qdiss.eom.fullci at 1e-15,
partitioned by which moment they carry: ("G",2) terms form the map under test, ("G",3) terms are
the ones a closure would have to put into sum_i L_i X L_i^dag form.

CHECK 2. Where D's kernel is empty -- which the clearance data says is the generic mixed-state
situation -- is the tangency machinery still needed for the OTHER N-representability cones? The
three two-index conditions, all Gram matrices and so manifestly positive semidefinite, are

    D[(p,v),(q,w)] = <a^dag_p a^dag_v a_w a_q>      = <zeta_vp, zeta_wq>,  zeta_ab = a_a a_b |psi>
    Q[(p,v),(q,w)] = <a_p a_v a^dag_w a^dag_q>      = <eta_vp,  eta_wq >,  eta_ab  = a^dag_a a^dag_b |psi>
    G[(p,v),(q,w)] = <a^dag_p a_v a^dag_w a_q>      = <xi_vp,   xi_wq  >,  xi_ab   = a^dag_a a_b |psi>

`D` and `Q` are indexed by unordered pairs p<v (dimension C(M,2)); `G` by all ordered pairs
(dimension M^2). They are computed here DIRECTLY as Gram matrices of explicitly constructed
vectors, not from index gymnastics on gamma_2, so no convention error is possible. For a mixed
state rho = sum_k w_k |psi_k><psi_k| each is the w-weighted sum of the pure-state Grams, still
positive semidefinite.
"""
import sys
import numpy as np
sys.path.insert(0, "."); sys.path.insert(0, "tools")
import hylleraas as H
from molecular_2rdm import build_H_sector, sector
from qdiss.eom import fullci, wick

M, N, h, V, enuc = H.M, H.N, H.h, H.V, H.enuc
EYE = np.eye(M, dtype=complex)
PAIRS = [(p, v) for p in range(M) for v in range(p + 1, M)]
NP = len(PAIRS)

# ---------------- check 1 : Lyapunov shape of the gamma_2-dependent part ------------------------
one, two = wick.terms(2)
G2, G3 = ([], []), ([], [])
for which, lst in ((0, one), (1, two)):
    for s, f in lst:
        keys = [(n, r) for n, r, u, l in f]
        subs = ["".join(u) + "".join(l) for n, r, u, l in f]
        if ("G", 3) in keys:
            G3[which].append((s, subs, keys))
        elif ("G", 2) in keys:
            G2[which].append((s, subs, keys))


def build(terms, tab):
    tot = np.zeros((M,) * 4, dtype=complex)
    for s, subs, keys in terms[0]:
        tot += s * np.einsum(",".join(["rs"] + subs) + "->pvqw", h,
                             *[tab[k] for k in keys], optimize=True)
    for s, subs, keys in terms[1]:
        tot += s * 0.25 * np.einsum(",".join(["rtsu"] + subs) + "->pvqw", V,
                                    *[tab[k] for k in keys], optimize=True)
    return 1j * tot


def pair_to_tensor(Z):
    T = np.zeros((M,) * 4, dtype=complex)
    for i, (p, v) in enumerate(PAIRS):
        for j, (q, w) in enumerate(PAIRS):
            x = 0.5 * Z[i, j]
            T[p, v, q, w] = x; T[v, p, q, w] = -x
            T[p, v, w, q] = -x; T[v, p, w, q] = x
    return T


def tensor_to_pair(T):
    return np.array([[2.0 * T[p, v, q, w] for (q, w) in PAIRS] for (p, v) in PAIRS])


print("CHECK 1  does the gamma_2-dependent part of gamma_2-dot have the shape A X + X A^dag ?")
print("   term counts: %d one-body + %d two-body carry gamma_2;  %d + %d carry gamma_3"
      % (len(G2[0]), len(G2[1]), len(G3[0]), len(G3[1])))
rng = np.random.default_rng(0)
g1r = rng.normal(size=(M, M)) + 1j * rng.normal(size=(M, M)); g1r = 0.5 * (g1r + g1r.conj().T)
print("   rank-1 probe: feed X = |Vec><Vec|, a Lyapunov map must return rank <= 2")
print("      %6s %14s %14s" % ("trial", "rank of output", "sing. values 3..5"))
for t in range(4):
    vv = rng.normal(size=NP) + 1j * rng.normal(size=NP)
    vv /= np.linalg.norm(vv)
    Z = np.outer(vv, vv.conj())
    out = tensor_to_pair(build(G2, {("D", 1): EYE, ("G", 1): g1r, ("G", 2): pair_to_tensor(Z)}))
    sv = np.linalg.svd(out, compute_uv=False)
    print("      %6d %14d %14s" % (t, int((sv > 1e-10 * sv.max()).sum()),
                                   np.array2string(sv[2:5] / sv.max(), precision=3)))
print()

# ---------------- check 2 : D, Q, G kernels across states ---------------------------------------
dets = sector(M, N)
es, vec = np.linalg.eigh(build_H_sector(h, V, M, dets))


def grams(w):
    keep = np.where(w > 1e-13)[0]
    D = np.zeros((NP, NP), dtype=complex)
    Q = np.zeros((NP, NP), dtype=complex)
    G = np.zeros((M * M, M * M), dtype=complex)
    for k in keep:
        psi = np.zeros(1 << M, dtype=complex); psi[dets] = vec[:, k]
        zeta = [fullci.aa(p, fullci.aa(v, psi)) for (p, v) in PAIRS]
        eta = [fullci.ad(p, fullci.ad(v, psi)) for (p, v) in PAIRS]
        xi = [fullci.ad(p, fullci.aa(v, psi)) for p in range(M) for v in range(M)]
        Zm = np.array(zeta); Em = np.array(eta); Xm = np.array(xi)
        D += w[k] * (Zm.conj() @ Zm.T)
        Q += w[k] * (Em.conj() @ Em.T)
        G += w[k] * (Xm.conj() @ Xm.T)
    return D, Q, G


print("CHECK 2  smallest eigenvalue of D, Q, G -- is D's kernel empty while others are not?")
print("   %-18s %8s %12s %12s %12s" % ("state", "eff rank", "min eig D", "min eig Q", "min eig G"))
cases = [("pure ground", np.eye(len(dets))[0])]
for beta in (10.0, 5.0, 2.0):
    w = np.exp(-beta * (es - es[0])); w /= w.sum()
    cases.append(("thermal b=%.0f" % beta, w))
for K in (8, 16):
    w = np.zeros(len(dets)); w[35:35 + K] = 1.0 / K
    cases.append(("hot-K K=%d" % K, w))
for tag, w in cases:
    D, Q, G = grams(w)
    f = lambda A: np.linalg.eigvalsh(0.5 * (A + A.conj().T)).min()
    print("   %-18s %8.2f %12.3e %12.3e %12.3e"
          % (tag, 1.0 / (w ** 2).sum(), f(D), f(Q), f(G)))
