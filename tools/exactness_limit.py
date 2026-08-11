#!/usr/bin/env python3
"""Does the EXACT state satisfy every constraint simultaneously, in the limit where the closure
becomes exact? Run this before concluding any constraint set is empty.

WHY. The constraint set under construction is {trace identity} and {connectedness} and {positive
generator}. If the exact solution fails any of them, the set is wrong no matter how appealing the
constraints are. The clean place to check is the top of the hierarchy: for an N-electron state,

    gamma_{N+1}[p_1..p_{N+1}, q_1..q_{N+1}] = <a+_p1 .. a+_p{N+1} a_q{N+1} .. a_q1>  ==  0

identically, because N+1 annihilation operators kill an N-particle state. So at k = N the 4k terms
carrying gamma_{k+1} vanish, the rank-N equation of motion CLOSES WITH NO APPROXIMATION, and it
must reproduce the von Neumann equation. Every constraint can then be tested against an object
that is known to be right.

WHAT IS CHECKED, on H3 bent / sto-3g, M = 6, N = 3, at rank k = 3 = N.

  1. gamma_4 == 0 from qdiss.eom.fullci. Confirms conventions before anything is built on them.
  2. The rank-3 term list from qdiss.eom.wick, evaluated with ("G",4) set to zero, reproduces
     qdiss.eom.fullci.rdm_dot exactly. This is the statement that the hierarchy closes at k = N.
  3. The exact lambda_3 satisfies the single-trace identity Tr_3 lambda_3 = anti(E, lambda_2),
     with E = gamma_1 - I/2 and Tr_3 T[pv,qw] = sum_y T[p v y, q w y].
  4. The exact lambda_3 is extensive: zero cross-fragment elements on a composite of
     non-interacting fragments.
  5. THE POSITIVE GENERATOR. The exact rank-3 flow, as a superoperator on
     Mat_3(T)[I,J] = 3! T[I+J] over sorted triples (dimension C(6,3) = 20), is tested for
     conditional complete positivity. A generator L is a GKSL generator iff it is
     Hermiticity-preserving, trace-annihilating, and
         Q Choi(L) Q  >=  0,     Choi(L) = sum_ij |i><j| (x) L[|i><j|],
         Q = 1 - |Omega><Omega|, |Omega> = (1/sqrt(n)) sum_i |i>|i>
     The eigenvalues of Q Choi(L) Q ARE the Kossakowski spectrum. For the exact flow at k = N the
     generator is the von Neumann commutator, so the Kossakowski spectrum must be zero to
     round-off. If it is not, the machinery is wrong and no closure result computed with it can be
     trusted.

THE CONTRAST, at rank k = 2 < N where a closure is actually needed: the same Kossakowski spectrum
for the linearised rank-2 flow under each lambda_3 reconstruction. The exact reduced flow at
k < N is NOT autonomous -- gamma_3 is independent data, not a function of gamma_2 -- so there is
no exact superoperator to compare against, and the reference is instead the drift alone, which is
the Lyapunov map A X + X A^dag and must also come out at zero.

Ground truth throughout: qdiss.eom.fullci.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from ansatz_hierarchy import fro, rank3_set
from molecular_2rdm import build_H_sector, molecule, sector
from qdiss import ansatz as az
from qdiss.eom import evaluate, fullci, wick

H3_BENT = "H 0.00 0.00 0.0; H 1.00 0.00 0.0; H 0.35 0.90 0.0"
N = 3

mol, mf, M, h, V, enuc = molecule(H3_BENT, "sto-3g", 0, 1)
dets = sector(M, N)
es, vec = np.linalg.eigh(build_H_sector(h, V, M, dets))
# A NON-STATIONARY state. The FCI ground state is an eigenstate, so d/dt gamma_k = 0 identically
# and every comparison against it is zero against zero -- a vacuous test. A normalised
# superposition is still a pure N-electron state, so gamma_{N+1} = 0 continues to hold exactly.
_rng = np.random.default_rng(7)
_c = vec[:, 0] + 0.4 * (_rng.normal(size=len(dets)) + 1j * _rng.normal(size=len(dets)))
_c /= np.linalg.norm(_c)
psi = np.zeros(1 << M, dtype=complex)
psi[dets] = _c
H = fullci.build_H(h, V, M)
g = {k: fullci.rdm(psi, M, k) for k in (1, 2, 3, 4)}
l2 = az.cumulant(2, g[2], {("G", 1): g[1]})
l3 = az.cumulant(3, g[3], {("G", 1): g[1], ("L", 2): l2})

print("=== exactness limit: H3 bent / sto-3g, M=%d N=%d, rank k = %d = N ===" % (M, N, N))
print("    E_FCI %.9f" % (es[0].real + enuc))
print()
print("1. gamma_{N+1} must vanish identically")
print("   max |gamma_4| = %.3e   (|gamma_3| = %.5f for scale)"
      % (np.abs(g[4]).max(), np.abs(g[3]).max()))

print()
print("2. the rank-N equation of motion closes with gamma_{N+1} = 0")
one3, two3 = wick.terms(3)
keys = sorted({(n, r) for lst in (one3, two3) for s, f in lst for n, r, u, l in f})
print("   rank-3 wick term list uses tensors: %s" % ", ".join("%s_%d" % k for k in keys))
tab = {("D", 1): np.eye(M, dtype=complex)}
for k in (1, 2, 3):
    tab[("G", k)] = g[k]
tab[("G", 4)] = np.zeros((M,) * 8, dtype=complex)
approx = wick.gamma_dot(3, h, V, tab, M)
truth = fullci.rdm_dot(psi, H, M, 3)
print("   |wick(gamma_4 = 0) - fullci d/dt gamma_3| = %.3e   (|d/dt gamma_3| = %.5f)"
      % (np.abs(approx - truth).max(), np.abs(truth).max()))
tab[("G", 4)] = g[4]
print("   |wick(gamma_4 = exact) - fullci|          = %.3e   (identical, since gamma_4 = 0)"
      % np.abs(wick.gamma_dot(3, h, V, tab, M) - truth).max())

print()
print("3. the exact lambda_3 satisfies the single-trace identity")
S = az.trace_target(3, g[1], l2)
print("   |Tr_3 lambda_3 - anti(E, lambda_2)| = %.3e   (|anti(E,lambda_2)| = %.5f)"
      % (az.trace_residual(l3, S, 3), np.abs(S).max()))

print()
print("4. the exact lambda_3 is extensive (checked in tools/ansatz_extensivity.py at 8.9e-16)")

# ---------------------------------------------------------------- the positive generator
TUP = az.tuples(M, 3)
n3 = len(TUP)


def mat3(T):
    return 6.0 * az.tensor_to_coord(T, M, 3)


def unmat3(A):
    return az.coord_to_tensor(A / 6.0, M, 3)


def superop(flow, n, to_mat, from_mat):
    """Matrix of a linear superoperator on n x n matrices, columns indexed by |i><j|."""
    L = np.zeros((n * n, n * n), dtype=complex)
    for i in range(n):
        for j in range(n):
            E = np.zeros((n, n), dtype=complex)
            E[i, j] = 1.0
            L[:, i * n + j] = to_mat(flow(from_mat(E))).reshape(-1)
    return L


def kossakowski(L, n):
    """Eigenvalues of Q Choi(L) Q, the Kossakowski spectrum. GKSL iff all are >= 0."""
    C = np.zeros((n * n, n * n), dtype=complex)
    for i in range(n):
        for j in range(n):
            LE = L[:, i * n + j].reshape(n, n)
            for k in range(n):
                for m in range(n):
                    C[i * n + k, j * n + m] = LE[k, m]
    om = np.zeros(n * n, dtype=complex)
    for i in range(n):
        om[i * n + i] = 1.0 / np.sqrt(n)
    Q = np.eye(n * n, dtype=complex) - np.outer(om, om.conj())
    K = Q @ C @ Q
    K = 0.5 * (K + K.conj().T)
    return np.linalg.eigvalsh(K), np.abs(C).max()


# gamma_4 is identically zero here, so its terms are dropped rather than contracting a 6^8 tensor
# once per superoperator column. Verified equivalent: check 2 above already showed that setting
# ("G",4) to zero reproduces fullci exactly.
_Z4 = np.zeros((M,) * 8, dtype=complex)
ONE3_C = [(sg, f) for sg, f in one3 if not any(n == "G" and r == 4 for n, r, u, l in f)]
TWO3_C = [(sg, f) for sg, f in two3 if not any(n == "G" and r == 4 for n, r, u, l in f)]


def flow3(T):
    """The rank-3 flow with the gamma_4 terms removed. Uses the validated assembler from
    qdiss.eom.evaluate, so the 1j and the TWO_BODY_PREFACTOR are imported, never re-chosen."""
    t = {("D", 1): np.eye(M, dtype=complex), ("G", 3): T}
    return evaluate.assemble(ONE3_C, TWO3_C, h, V, t, "pvxqwy", M)


print()
print("5. POSITIVE GENERATOR of the exact rank-N flow, on Mat_3 (dimension %d)" % n3)
print("   dropped %d of %d two-body terms carrying gamma_4; cross-check against the full"
      % (len(two3) - len(TWO3_C), len(two3)))
print("   assembler on the exact state: |flow3(gamma_3) - fullci| = %.3e"
      % np.abs(flow3(g[3]) - truth).max())
L3 = superop(flow3, n3, mat3, unmat3)
# machinery preconditions
X = np.random.default_rng(3).normal(size=(n3, n3)) + 1j * np.random.default_rng(4).normal(size=(n3, n3))
X = 0.5 * (X + X.conj().T)
Y = (L3 @ X.reshape(-1)).reshape(n3, n3)
print("   Hermiticity preserving : |L[X] - L[X]^dag| = %.3e" % np.abs(Y - Y.conj().T).max())
print("   trace annihilating     : |tr L[X]|         = %.3e" % abs(np.trace(Y)))
ev, cscale = kossakowski(L3, n3)
print("   Kossakowski spectrum   : min %.3e   max %.3e   (|Choi|_max %.3e)"
      % (ev.min(), ev.max(), cscale))
print("   -> the exact rank-N generator is the von Neumann commutator, so this must be zero;")
print("      any nonzero value here is a bug in the machinery, not physics.")

# ---------------------------------------------------------------- rank 2, where a closure is needed
PAIR = az.tuples(M, 2)
n2 = len(PAIR)
one2, two2 = wick.terms(2)


def mat2(T):
    return 2.0 * az.tensor_to_coord(T, M, 2)


def unmat2(A):
    return az.coord_to_tensor(A / 2.0, M, 2)


def flow2_drift(T):
    """The gamma_2-carrying terms alone: the Lyapunov drift, with gamma_3 held at zero."""
    t = {("D", 1): np.eye(M, dtype=complex), ("G", 1): g[1], ("G", 2): T,
         ("G", 3): np.zeros((M,) * 6, dtype=complex)}
    return wick.gamma_dot(2, h, V, t, M)


print()
print("6. rank 2, where a closure is required. Kossakowski spectrum of the linearised flow.")
print("   %-16s %13s %13s %13s" % ("generator", "min eig", "max eig", "|Choi|_max"))
ev, cs = kossakowski(superop(flow2_drift, n2, mat2, unmat2), n2)
print("   %-16s %13.3e %13.3e %13.3e" % ("drift only", ev.min(), ev.max(), cs))

g2m = mat2(g[2])
eps = 1e-5
for nm in ("null", "trace", "scaled", "scaled+trace"):
    def flow2(T, nm=nm):
        g1 = np.einsum("pxqx->pq", T) / (N - 1)
        lam2 = az.cumulant(2, T, {("G", 1): g1})
        g3 = az.separable(3, {("G", 1): g1, ("L", 2): lam2}) + rank3_set(g1, lam2, N, M)[nm]
        t = {("D", 1): np.eye(M, dtype=complex), ("G", 1): g1, ("G", 2): T, ("G", 3): g3}
        return wick.gamma_dot(2, h, V, t, M)

    f0 = mat2(flow2(g[2]))
    J = np.zeros((n2 * n2, n2 * n2), dtype=complex)
    for i in range(n2):
        for j in range(n2):
            E = np.zeros((n2, n2), dtype=complex)
            E[i, j] = 1.0
            J[:, i * n2 + j] = (mat2(flow2(unmat2(g2m + eps * E))) - f0).reshape(-1) / eps
    ev, cs = kossakowski(J, n2)
    print("   %-16s %13.3e %13.3e %13.3e" % (nm, ev.min(), ev.max(), cs))
print()
print("   Read these against the drift row, not against zero: the exact CLOSED flow is tangent")
print("   at the PSD boundary, not inward, so it is not CP-divisible with a nonzero gain part.")
