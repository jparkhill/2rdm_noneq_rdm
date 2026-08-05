#!/usr/bin/env python3
"""Analytic confirmation that lambda_3 does not vanish for two electrons.

Take the two-configuration state |Psi> = c1|12> + c2|34> at M = 4. Then gamma_1 is
diagonal, diag(|c1|^2, |c1|^2, |c2|^2, |c2|^2), and the only surviving pieces of the
[123;123] element are the three terms with matched indices:

    A[123;123] = det gamma_1 on {1,2,3}                       = |c1|^4 |c2|^2
    lambda_2[23;23] = 0 - |c1|^2|c2|^2                        = -|c1|^2|c2|^2
    lambda_2[13;13] = -|c1|^2|c2|^2 ,  lambda_2[12;12] = +|c1|^2|c2|^2
    B[123;123] = |c1|^2 lam[23;23] + |c1|^2 lam[13;13] + |c2|^2 lam[12;12]
               = -2|c1|^4|c2|^2 + |c1|^2|c2|^4

    A + B = |c1|^2 |c2|^2 ( |c2|^2 - |c1|^2 )

which is the exact lambda_3 up to sign, since gamma_3 = 0. It is nonzero for every state
except the equal-weight one, so lambda_3 = 0 is not a property of two-electron systems and
no coefficient in the expansion can make it one.

Also included: the N = 1 control. There gamma_2 = 0 and gamma_1 is rank one, so
gamma_1 wedge gamma_1 vanishes identically and lambda_2 must come out exactly zero. If the
machinery got that wrong the rank-2 expansion would be suspect.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
from qdiss.eom import fullci
from qdiss.symbolic import mkno

UP3, LO3 = ("p", "v", "x"), ("q", "w", "y")
ORB = {"A": [], "B": []}
for s, f in mkno.gamma_to_cumulants(UP3, LO3):
    if any(n == "L" and r == 3 for n, r, u, l in f):
        continue
    k = "B" if any(n == "L" and r == 2 for n, r, u, l in f) else "A"
    ORB[k].append((s, ",".join("".join(u) + "".join(l) for n, r, u, l in f) + "->pvxqwy",
                   [(n, r) for n, r, u, l in f]))


def wedge(a, b):
    return np.einsum("pq,vw->pvqw", a, b) - np.einsum("pw,vq->pvqw", a, b)


def orbit(key, g1, l2, M):
    t = {("G", 1): g1, ("L", 2): l2}
    out = np.zeros((M,) * 6, dtype=complex)
    for s, sub, keys in ORB[key]:
        out = out + s * np.einsum(sub, *[t[k] for k in keys], optimize=True)
    return out


print("=== N = 1 control: lambda_2 must vanish exactly ===")
for M in (4, 6):
    rng = np.random.default_rng(3)
    psi = fullci.random_state(M, 1, rng)
    g1 = fullci.rdm(psi, M, 1)
    g2 = fullci.rdm(psi, M, 2)
    l2 = g2 - wedge(g1, g1)
    print("   M=%d  rank(gamma_1)=%d  |gamma_2|=%.3e  |gamma_1^gamma_1|=%.3e  |lambda_2|=%.3e"
          % (M, np.linalg.matrix_rank(g1, tol=1e-9), np.abs(g2).max(),
             np.abs(wedge(g1, g1)).max(), np.abs(l2).max()))
print()

print("=== N = 2, two-configuration state: measured A+B vs the closed form ===")
M = 4
dim = 1 << M
print("%8s %8s | %14s %14s | %12s"
      % ("|c1|^2", "|c2|^2", "measured A+B", "|c1c2|^2 dn", "|gamma_3|"))
for x in (0.1, 0.25, 0.5, 0.7, 0.9):
    c1, c2 = np.sqrt(x), np.sqrt(1 - x)
    psi = np.zeros(dim, dtype=complex)
    psi[0b0011] = c1
    psi[0b1100] = c2
    psi = psi / np.linalg.norm(psi)
    g1 = fullci.rdm(psi, M, 1)
    g2 = fullci.rdm(psi, M, 2)
    g3 = fullci.rdm(psi, M, 3)
    l2 = g2 - wedge(g1, g1)
    tot = orbit("A", g1, l2, M) + orbit("B", g1, l2, M)
    meas = tot[0, 1, 2, 0, 1, 2].real
    pred = (x * (1 - x) * ((1 - x) - x))
    print("%8.3f %8.3f | %14.9f %14.9f | %12.3e"
          % (x, 1 - x, meas, pred, np.abs(g3).max()))
print()
print("The [123;123] element matches |c1|^2|c2|^2(|c2|^2-|c1|^2) exactly, and vanishes only")
print("at equal weights. gamma_3 is zero throughout, so this IS -lambda_3.")
print()
print("Conclusion: lambda_3 is genuinely nonzero for two electrons. gamma_3 = 0 is exact")
print("there and lambda_3 = 0 is not, so MK_lambda_30 cannot reproduce N=2 dynamics.")
