#!/usr/bin/env python3
"""Check the AGP geminal occupation formula, and whether anything beats it.

Claimed derivation for the equal-weight AGP condensate with K = M/2 pair levels
and n = N/2 pairs:

    a = <P+_k P_k>   = n/K                      (level k occupied)
    b = <P+_k P_l>   = n(K-n)/[K(K-1)]  (k != l, the coherence)
    <G+G>            = a + (K-1) b = n(K-n+1)/K
    eigenvalue       = 2<G+G> = N(M-N+2)/M

Verified three ways: a and b measured directly, <G+G> against the formula, and the
eigenvalue against 2<G+G>. Then random states are sampled to see whether any exceeds
the AGP value, which is evidence about whether it is the supremum -- evidence, not
proof.
"""
import numpy as np


def run(M, N, ntrial=300):
    dim = 1 << M
    K, n = M // 2, N // 2

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

    def Pdag(k, v):
        return ad(2 * k, ad(2 * k + 1, v))

    def P(k, v):
        return aa(2 * k + 1, aa(2 * k, v))

    vac = np.zeros(dim, dtype=complex); vac[0] = 1.0
    psi = vac
    for _ in range(n):
        acc = np.zeros_like(psi)
        for k in range(K):
            acc = acc + Pdag(k, psi)
        psi = acc
    psi = psi / np.linalg.norm(psi)

    a_meas = float(np.real(np.vdot(P(0, psi), P(0, psi))))
    b_meas = float(np.real(np.vdot(P(0, psi), P(1, psi)))) if K > 1 else 0.0
    a_pred = n / K
    b_pred = n * (K - n) / (K * (K - 1)) if K > 1 else 0.0

    Gpsi = np.zeros_like(psi)
    for k in range(K):
        Gpsi = Gpsi + P(k, psi) / np.sqrt(K)
    GdG = float(np.real(np.vdot(Gpsi, Gpsi)))
    GdG_pred = n * (K - n + 1) / K

    PAIRS = [(p, v) for p in range(M) for v in range(p + 1, M)]
    NP = len(PAIRS)

    def top_eig(state):
        kets = {b: aa(w, aa(q, state)) for b, (q, w) in enumerate(PAIRS)}
        D = np.zeros((NP, NP), dtype=complex)
        for i, (p, v) in enumerate(PAIRS):
            left = aa(v, aa(p, state))
            for j in range(NP):
                D[i, j] = 2.0 * np.vdot(left, kets[j])
        return float(np.linalg.eigvalsh(0.5 * (D + D.conj().T)).max())

    ev = top_eig(psi)
    formula = N * (M - N + 2) / M

    print("M=%d N=%d  (K=%d pair levels, n=%d pairs)" % (M, N, K, n))
    print("   a  measured %.6f   predicted n/K            %.6f" % (a_meas, a_pred))
    print("   b  measured %.6f   predicted n(K-n)/K(K-1)  %.6f" % (b_meas, b_pred))
    print("   <G+G> measured %.6f  predicted n(K-n+1)/K   %.6f" % (GdG, GdG_pred))
    print("   eigenvalue %.6f   2<G+G> %.6f   N(M-N+2)/M %.6f"
          % (ev, 2 * GdG, formula))

    popc = np.array([bin(s).count("1") for s in range(dim)])
    r = np.random.default_rng(3)
    best = 0.0
    for _ in range(ntrial):
        v = np.where(popc == N, r.normal(size=dim) + 1j * r.normal(size=dim), 0.0)
        best = max(best, top_eig(v / np.linalg.norm(v)))
    print("   best of %d random states: %.6f   %s"
          % (ntrial, best, "under AGP" if best < ev + 1e-9 else "EXCEEDS AGP"))
    print()


for M, N in ((4, 2), (6, 2), (8, 4), (10, 4)):
    run(M, N, ntrial=200 if M <= 8 else 60)
