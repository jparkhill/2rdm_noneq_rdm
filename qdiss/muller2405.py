"""The Mueller 3-RDM reconstruction, as a baseline.

Mueller, Luoma and Schaefer, "Hierarchical approach to quantum many-body systems in
non-Markovian environments", arXiv:2405.05093v3. Their scheme is BBGKY plus HEOM for open
systems; what is taken here is the closed-system limit of their 3-RDM reconstruction, so the
HEOM indices (n,m) are dropped and every superscripted F equals its unsuperscripted partner.

Normalisation matches ours: Tr F_1 = N, Tr F_12 = N(N-1), Tr F_123 = N(N-1)(N-2).

Their Eq. 14, with Tr(F_1) = N substituted and like terms collected, reduces to

    F~_123 = (N-2)/N [ F_12 F_3 + F_13 F_2 + F_1 F_23 ]
             - 2(N-1)(N-2)/N^2  F_1 F_2 F_3

By hand, Tr_3 of that is exactly (N-2) F_12, using Tr F_3 = N and Tr_3 F_13 = (N-1) F_1. So the
unantisymmetrised reconstruction is already contraction-consistent; their Eq. 44/45 exists to
repair what ANTISYMMETRISATION breaks, not the reconstruction itself.

Note (N-2) makes the whole thing vanish at N = 2, where the exact gamma_3 is identically zero.
So this reconstruction is exact for two electrons, which lambda_3 = 0 is not.

Two deviations from the paper are recorded rather than hidden. Their antisymmetriser is
Lambda- F~ Lambda- / norm with Lambda-_123 = (1-P_12)(1-P_13-P_23), which is 6 times the full
antisymmetriser; here the full antisymmetric projector is applied directly and the scale is
fixed by requiring Tr gamma_3 = N(N-1)(N-2), which is unambiguous. And their perpendicular
projection in Eq. 43, used to build the Eq. 45 correction, was not available in the text
obtained, so the correction is built with the closed-form partial-trace inverse derived in
tools/lambda3_build.py, which produces a tensor with any prescribed partial trace. That serves
the same purpose but is not their construction.
"""
import itertools
import numpy as np
from qdiss.symbolic import mkno

_UP3, _LO3 = ("p", "v", "x"), ("q", "w", "y")
_ORB_B = [(s, ",".join("".join(u) + "".join(l) for n, r, u, l in f) + "->pvxqwy",
           [(n, r) for n, r, u, l in f])
          for s, f in mkno.gamma_to_cumulants(_UP3, _LO3)
          if not any(n == "L" and r == 3 for n, r, u, l in f)
          and any(n == "L" and r == 2 for n, r, u, l in f)]


def wedge(a, b):
    return np.einsum("pq,vw->pvqw", a, b) - np.einsum("pw,vq->pvqw", a, b)


def antisymmetrize3(T):
    """Project a 6-index tensor onto antisymmetric in the upper and lower triples."""
    out = np.zeros_like(T)
    eye3 = np.eye(3)
    for pu in itertools.permutations(range(3)):
        su = np.linalg.det(eye3[list(pu)])
        for pl in itertools.permutations(range(3)):
            sl = np.linalg.det(eye3[list(pl)])
            out += su * sl * np.transpose(T, tuple(pu) + tuple(3 + i for i in pl))
    return out / 36.0


def _orbit_B(a, X, M):
    t = {("G", 1): a, ("L", 2): X}
    out = np.zeros((M,) * 6, dtype=complex)
    for s, sub, keys in _ORB_B:
        out = out + s * np.einsum(sub, *[t[k] for k in keys], optimize=True)
    return out


def trace_inverse(R, M):
    """A 6-index antisymmetric tensor whose partial trace over the third pair is R.

    Uses sum_y B[I,X] = (M-4) X + (I^S_X) + (S_X^I), verified in tools/lambda3_build.py.
    """
    I = np.eye(M, dtype=complex)
    S_R = np.einsum("pxqx->pq", R)
    trS = np.trace(S_R) / (3 * M - 6)
    S_X = (S_R - trS * I) / (2 * M - 6)
    X = (R - wedge(I, S_X) - wedge(S_X, I)) / (M - 4)
    return _orbit_B(I, X, M)


def gamma3(g1, g2, N, contraction_correct=True):
    """Mueller reconstruction of gamma_3 from gamma_1 and gamma_2."""
    M = g1.shape[0]
    if N <= 2:
        return np.zeros((M,) * 6, dtype=complex)
    prod = (np.einsum("pvqw,xy->pvxqwy", g2, g1)
            + np.einsum("pxqy,vw->pvxqwy", g2, g1)
            + np.einsum("pq,vxwy->pvxqwy", g1, g2))
    trip = np.einsum("pq,vw,xy->pvxqwy", g1, g1, g1)
    F = (N - 2) / N * prod - 2.0 * (N - 1) * (N - 2) / N ** 2 * trip
    F = antisymmetrize3(F)
    want = N * (N - 1) * (N - 2)
    got = np.einsum("pvxpvx->", F).real
    if abs(got) > 1e-14:
        F = F * (want / got)
    if not contraction_correct:
        return F
    delta = g2 - np.einsum("pvxqwx->pvqw", F) / (N - 2)
    return F + trace_inverse((N - 2) * delta, M)
