#!/usr/bin/env python3
"""Both constraints on lambda_3 are LINEAR, so no SDP is needed.

The exact flow has X_KK = 0 exactly, not merely PSD: rank gamma_2 <= C(M, N-2) is conserved,
so the kernel of gamma_2 is invariant and the kernel block of gamma_2-dot vanishes identically.
Measured at 4.2e-15 for the exact lambda_3. Since X_KK is linear in lambda_3, that is an affine
condition, and the partial-trace identity is affine too. Their intersection contains the exact
lambda_3, so it is non-empty, and a member can be found by linear least squares instead of a
cone program.

Construction: start from the closed-form trace-exact lambda_3 and add a combination of
TRACE-FREE antisymmetric basis tensors, which cannot disturb the trace, chosen to cancel X_KK.
The kernel is nontrivial only when C(M,N-2) < C(M,2), i.e. for two and three electrons, which
is exactly the regime of interest here.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from lambda3_build import (CUM2, EYE, Kb, M, N, ORB, R, V, assemble, g1, g2, h, invert,
                           l2, l3_tc, l3x, orbit, ptrace, to_pair, wedge)

NB = 150


def xkk(l3):
    t = {("D", 1): EYE, ("G", 1): g1, ("L", 2): l2, ("L", 3): l3}
    X = to_pair(assemble(CUM2[0], CUM2[1], h, V, t, "pvqw", M))
    B = Kb.conj().T @ X @ Kb
    return 0.5 * (B + B.conj().T)


def antisym6(T):
    import itertools
    out = np.zeros_like(T)
    for pu in itertools.permutations(range(3)):
        su = np.linalg.det(np.eye(3)[list(pu)])
        for pl in itertools.permutations(range(3)):
            sl = np.linalg.det(np.eye(3)[list(pl)])
            out += su * sl * np.transpose(T, tuple(pu) + tuple(3 + i for i in pl))
    return out / 36.0


print("kernel dim %d, so X_KK has %d real degrees of freedom"
      % (Kb.shape[1], Kb.shape[1] ** 2))
print("exact lambda_3 : |X_KK| = %.3e   (zero, not merely PSD)" % np.abs(xkk(l3x)).max())
print()

rng = np.random.default_rng(21)
basis, cols = [], []
X0 = xkk(l3_tc)
k = Kb.shape[1]
iu = np.triu_indices(k)


def flat(B):
    return np.concatenate([B[iu].real, B[iu].imag])


for i in range(NB):
    T = rng.normal(size=(M,) * 6) + 1j * rng.normal(size=(M,) * 6)
    T = antisym6(T)
    T = T - orbit("B", EYE, invert(ptrace(T), M), M)
    T = T / np.abs(T).max()
    basis.append(T)
    cols.append(flat(xkk(l3_tc + T) - X0))
A = np.stack(cols, axis=1)
sol, res, rank, sv = np.linalg.lstsq(A, -flat(X0), rcond=None)
l3_lin = l3_tc + sum(sol[i] * basis[i] for i in range(NB))

print("least squares over %d trace-free directions, rank %d of %d equations"
      % (NB, rank, A.shape[0]))
print()
print("%-28s | %12s %12s %10s" % ("lambda_3", "|X_KK|", "trace resid", "|lambda_3|"))
for nm, l3 in (("exact", l3x),
               ("zero", np.zeros((M,) * 6, dtype=complex)),
               ("-(A+B) (MK_gamma30)", -(orbit("A", g1, l2, M) + orbit("B", g1, l2, M))),
               ("closed-form trace-exact", l3_tc),
               ("+ linear positivity solve", l3_lin)):
    print("%-28s | %12.4e %12.2e %10.4f"
          % (nm, np.abs(xkk(l3)).max(), np.abs(ptrace(l3) - R).max(), np.abs(l3).max()))
print()
mn = np.linalg.eigvalsh(xkk(l3_lin)).min()
print("min eig X_KK for the linear solution : %.3e" % mn)
print("|lambda_3 solution - exact|          : %.4f  (|exact| %.4f)"
      % (np.abs(l3_lin - l3x).max(), np.abs(l3x).max()))
