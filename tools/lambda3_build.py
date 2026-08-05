#!/usr/bin/env python3
"""First stab at a lambda_3 reconstruction: exact partial trace in closed form, then measure
how far the resulting 2-RDM flow is from positivity.

Derivation: qdiss.eom.gwt cumulant form, validated against qdiss.eom.fullci at 1e-15.
Propagated state for the rank-2 EOM is (gamma_1, lambda_2); lambda_3 is the object being
reconstructed.

The trace target is exact and already derived (tools/lambda3_factorization.py, verified 1e-16):

    sum_y lambda_3[pvy,qwy] = E lambda_2 + lambda_2 E ,   E = gamma_1 - 1/2

Inverting the partial trace. Let B[a,X] be the 9-term gamma_1-wedge-lambda_2 orbit evaluated
with the 2-index slot a and the 4-index slot X. By hand,

    sum_y B[I,X] = (M-4) X + (I^S_X) + (S_X^I) ,   S_X = sum_x X[.x,.x]

using sum_x (a^b)[px,qx] = a tr(b) - a b, with the four contraction terms collapsing to -4X
when a = I. Contracting that relation again fixes S_X from the target's own trace:

    (2M-6) S_X + tr(S_X) I = S_R      ==>   tr(S_X) = tr(S_R)/(3M-6),
                                            S_X = [S_R - tr(S_X) I]/(2M-6)

and then X = [R - (I^S_X) - (S_X^I)]/(M-4), so lambda_3 = B[I,X] has partial trace exactly R.
Every step is checked numerically here rather than trusted.

Positivity is then measured by the exact Nagumo condition for the PSD cone, which is also the
CP condition at a point: with K the kernel of gamma_2, the flow keeps gamma_2 >= 0 iff the
kernel block X_KK of gamma_2-dot is PSD. Benchmarks from tools/kernel_psd.py at M=6 N=3:
exact gamma_3 gives -4.2e-15, gamma_3 = 0 gives -4.1e-15, lambda_3 = 0 gives -1.73.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from qdiss.eom import fullci, gwt
from qdiss.eom.evaluate import assemble
from qdiss.eom.hamiltonian import random_V, random_h
from qdiss.symbolic import mkno

UP3, LO3 = ("p", "v", "x"), ("q", "w", "y")
ORB = {"A": [], "B": []}
for s, f in mkno.gamma_to_cumulants(UP3, LO3):
    if any(n == "L" and r == 3 for n, r, u, l in f):
        continue
    k = "B" if any(n == "L" and r == 2 for n, r, u, l in f) else "A"
    ORB[k].append((s, ",".join("".join(u) + "".join(l) for n, r, u, l in f) + "->pvxqwy",
                   [(n, r) for n, r, u, l in f]))
CUM2 = gwt.terms_cumulant(2)


def wedge(a, b):
    return np.einsum("pq,vw->pvqw", a, b) - np.einsum("pw,vq->pvqw", a, b)


def orbit(key, a, X, M):
    t = {("G", 1): a, ("L", 2): X}
    out = np.zeros((M,) * 6, dtype=complex)
    for s, sub, keys in ORB[key]:
        out = out + s * np.einsum(sub, *[t[k] for k in keys], optimize=True)
    return out


def ptrace(T):
    return np.einsum("pvyqwy->pvqw", T)


def ptrace2(X):
    return np.einsum("pxqx->pq", X)


def invert(R, M):
    """X with sum_y B[I,X] = R."""
    I = np.eye(M, dtype=complex)
    S_R = ptrace2(R)
    trS = np.trace(S_R) / (3 * M - 6)
    S_X = (S_R - trS * I) / (2 * M - 6)
    return (R - wedge(I, S_X) - wedge(S_X, I)) / (M - 4)


M, N = 6, 3
EYE = np.eye(M, dtype=complex)
rng = np.random.default_rng(11)
h, V = random_h(M, rng), random_V(M, rng)
psi = fullci.random_state(M, N, rng)
g1 = fullci.rdm(psi, M, 1)
g2 = fullci.rdm(psi, M, 2)
g3 = fullci.rdm(psi, M, 3)
l2 = g2 - wedge(g1, g1)
l3x = g3 - orbit("A", g1, l2, M) - orbit("B", g1, l2, M)

print("=== step 1: verify sum_y B[I,X] = (M-4)X + (I^S_X) + (S_X^I) ===")
r2 = np.random.default_rng(3)
Xt = r2.normal(size=(M,) * 4) + 1j * r2.normal(size=(M,) * 4)
Xt = Xt - np.transpose(Xt, (1, 0, 2, 3))
Xt = Xt - np.transpose(Xt, (0, 1, 3, 2))
S = ptrace2(Xt)
pred = (M - 4) * Xt + wedge(EYE, S) + wedge(S, EYE)
print("   residual %.3e   (scale %.4f)"
      % (np.abs(ptrace(orbit("B", EYE, Xt, M)) - pred).max(), np.abs(pred).max()))

print()
print("=== step 2: the closed-form trace-exact lambda_3 ===")
E = g1 - 0.5 * EYE
R = (np.einsum("py,yvqw->pvqw", E, l2) + np.einsum("vy,pyqw->pvqw", E, l2)
     + np.einsum("pvyw,yq->pvqw", l2, E) + np.einsum("pvqy,yw->pvqw", l2, E))
print("   |R - sum_y lambda_3 exact| = %.3e  (this is the verified identity)"
      % np.abs(R - ptrace(l3x)).max())
l3_tc = orbit("B", EYE, invert(R, M), M)
print("   |sum_y lambda_3_tc - R|    = %.3e" % np.abs(ptrace(l3_tc) - R).max())
au = np.abs(l3_tc + np.transpose(l3_tc, (1, 0, 2, 3, 4, 5))).max()
al = np.abs(l3_tc + np.transpose(l3_tc, (0, 1, 2, 4, 3, 5))).max()
print("   antisymmetry upper %.1e lower %.1e   |lambda_3_tc| %.4f (exact %.4f)"
      % (au, al, np.abs(l3_tc).max(), np.abs(l3x).max()))

print()
print("=== step 3: positivity of the resulting gamma_2 flow ===")
PAIRS = [(p, v) for p in range(M) for v in range(p + 1, M)]


def to_pair(T):
    A = np.array([[2.0 * T[p, v, q, w] for (q, w) in PAIRS] for (p, v) in PAIRS])
    return 0.5 * (A + A.conj().T)


rho = to_pair(g2)
lam, U = np.linalg.eigh(rho)
Kb = U[:, lam < 1e-9]


def kk(l3):
    t = {("D", 1): EYE, ("G", 1): g1, ("L", 2): l2, ("L", 3): l3}
    X = to_pair(assemble(CUM2[0], CUM2[1], h, V, t, "pvqw", M))
    B = Kb.conj().T @ X @ Kb
    return np.linalg.eigvalsh(0.5 * (B + B.conj().T)).min()


Z = np.zeros((M,) * 6, dtype=complex)
print("   kernel dim %d of %d" % (Kb.shape[1], len(PAIRS)))
print("   %-26s | %12s %12s" % ("lambda_3", "min X_KK", "trace resid"))
for nm, l3 in (("exact", l3x), ("zero (lambda_3 = 0)", Z),
               ("-(A+B)  (MK_gamma30)",
                -(orbit("A", g1, l2, M) + orbit("B", g1, l2, M))),
               ("closed-form trace-exact", l3_tc)):
    print("   %-26s | %12.4e %12.2e"
          % (nm, kk(l3), np.abs(ptrace(l3) - R).max()))
print()
print("A violation of order the lambda_3 = 0 benchmark means the closed form needs different")
print("structure; a small one means an SDP or cone correction is a cheap patch.")
