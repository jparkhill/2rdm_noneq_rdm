#!/usr/bin/env python3
"""Does complete positivity in the HOLE representation characterize the
Q-inward generators?

Q is an affine function of the pair-space rho:  Q = C + A(rho), with A linear.
So the generator induced on Q is  Mh = A . L . A^{-1}, and the flow on Q is
CP-generated iff Mh is in GKSL form. The standard test is conditional complete
positivity: build the Choi matrix of Mh and require it to be PSD on the subspace
orthogonal to the maximally entangled vector.

If hole-CP holds exactly for the Q-inward generators and fails for the others,
that is the characterization, and it gives a constructive recipe -- build the
generator CP in both representations -- instead of a feasibility solve.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from qdiss.representability.dqg import contract_to_1rdm, build_Q, build_G

M, N = 6, 3
dim = 1 << M
eye = np.eye(M, dtype=complex)
PAIRS = [(p, v) for p in range(M) for v in range(p + 1, M)]
NP = len(PAIRS)
D2 = NP * NP


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


def to_pair(T):
    A = np.zeros((NP, NP), dtype=complex)
    for a, (p, v) in enumerate(PAIRS):
        for b, (q, w) in enumerate(PAIRS):
            A[a, b] = 2.0 * T[p, v, q, w]
    return A


def from_pair(A):
    T = np.zeros((M,) * 4, dtype=complex)
    for a, (p, v) in enumerate(PAIRS):
        for b, (q, w) in enumerate(PAIRS):
            x = 0.5 * A[a, b]
            T[p, v, q, w] = x
            T[v, p, q, w] = -x
            T[p, v, w, q] = -x
            T[v, p, w, q] = x
    return T


def qlin_tensor(dg2):
    """linear-in-gamma_2 part of Q (drops the constant delta.delta piece)."""
    dg1 = contract_to_1rdm(dg2, N)
    return (- np.einsum("ki,lj->ijkl", eye, dg1)
            + np.einsum("kj,li->ijkl", eye, dg1)
            + np.einsum("li,kj->ijkl", eye, dg1)
            - np.einsum("lj,ki->ijkl", eye, dg1)
            + np.einsum("klij->ijkl", dg2))


# ---- build A as a 225x225 matrix on vec(rho), column stacking ----
def vec(X):
    return X.reshape(-1, order="F")


def unvec(x):
    return x.reshape(NP, NP, order="F")


Amat = np.zeros((D2, D2), dtype=complex)
for k in range(D2):
    E = unvec(np.eye(D2, dtype=complex)[:, k])
    Amat[:, k] = vec(to_pair(qlin_tensor(from_pair(E))))
print("A (gamma_2 -> Q linear part) built; cond = %.3e" % np.linalg.cond(Amat))
Ainv = np.linalg.inv(Amat)


def superop(K, Ls):
    I = np.eye(NP, dtype=complex)
    S = -1j * (np.kron(I, K) - np.kron(K.T, I))
    for L in Ls:
        LdL = L.conj().T @ L
        S = S + (np.kron(L.conj(), L)
                 - 0.5 * np.kron(I, LdL) - 0.5 * np.kron(LdL.T, I))
    return S


def choi(S):
    """Choi matrix of the superoperator S (column-stacking vec)."""
    C = np.zeros((D2, D2), dtype=complex)
    for i in range(NP):
        for j in range(NP):
            E = np.zeros((NP, NP), dtype=complex)
            E[i, j] = 1.0
            out = unvec(S @ vec(E))
            C += np.kron(out, E)
    return C


Omega = np.zeros(D2, dtype=complex)
for i in range(NP):
    e = np.zeros(NP); e[i] = 1.0
    Omega += np.kron(e, e)
Omega /= np.linalg.norm(Omega)
P = np.eye(D2, dtype=complex) - np.outer(Omega, Omega.conj())


def ccp_min_eig(S):
    """conditional complete positivity: min eig of P Choi(S) P."""
    C = choi(S)
    Cp = P @ (0.5 * (C + C.conj().T)) @ P
    return float(np.linalg.eigvalsh(Cp).min())


# ---- reference state and null spaces ----
rng = np.random.default_rng(2026)
psi = np.where(np.array([bin(s).count("1") for s in range(dim)]) == N,
               rng.normal(size=dim) + 1j * rng.normal(size=dim), 0.0)
psi = psi / np.linalg.norm(psi)
g2 = rdm2(psi)
g1 = contract_to_1rdm(g2, N)
rho = to_pair(g2)


def mat(T):
    A = T.reshape(M * M, M * M)
    return 0.5 * (A + A.conj().T)


NULLQ = None
w, U = np.linalg.eigh(mat(build_Q(g2, g1)))
NULLQ = U[:, w < 1e-9]


def qflux(dg2):
    dot = qlin_tensor(dg2)
    return float(np.real(np.einsum("ai,ab,bi->i",
                                   NULLQ.conj(), mat(dot), NULLQ)).min())


print()
print("sanity: a GKSL generator must have ccp min eig >= 0 in the PARTICLE picture")
r0 = np.random.default_rng(1)
Ktest = 0.5 * (lambda A: A + A.conj().T)(r0.normal(size=(NP, NP)) + 1j * r0.normal(size=(NP, NP)))
Ltest = [0.3 * (r0.normal(size=(NP, NP)) + 1j * r0.normal(size=(NP, NP)))]
print("   particle-picture ccp min eig = %+.3e" % ccp_min_eig(superop(Ktest, Ltest)))
print()

print("%5s | %12s | %14s | %s" % ("run", "Q flux", "hole ccp min", "agree?"))
agree = 0
n = 24
rows = []
for run in range(n):
    r = np.random.default_rng(9000 + run)
    A = r.normal(size=(NP, NP)) + 1j * r.normal(size=(NP, NP))
    K = 0.5 * (A + A.conj().T)
    nL = 1 + run % 3
    Ls = [0.3 * (r.normal(size=(NP, NP)) + 1j * r.normal(size=(NP, NP)))
          for _ in range(nL)]
    S = superop(K, Ls)
    dg2 = from_pair(unvec(S @ vec(rho)))
    qf = qflux(dg2)
    Mh = Amat @ S @ Ainv
    hc = ccp_min_eig(Mh)
    qok = qf > -1e-9
    hok = hc > -1e-9
    ok = (qok == hok)
    agree += ok
    rows.append((run, qf, hc, qok, hok, ok))
    print("%5d | %12.4e | %14.4e | %s"
          % (run, qf, hc, "yes" if ok else "NO"))

print()
print("hole-CP agrees with Q-inwardness in %d/%d runs" % (agree, n))
nq = sum(1 for r in rows if r[3])
nh = sum(1 for r in rows if r[4])
print("  Q inward: %d   hole-CP: %d" % (nq, nh))
