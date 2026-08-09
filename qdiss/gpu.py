"""GPU backend: the same symbolic term lists, contracted with torch einsums on device.

The symbolic layer is untouched. `qdiss.eom.wick` and `qdiss.eom.gwt` produce term lists of
(sign, factors) with factor = (name, rank, upper, lower); those are just strings and integers, so
the only thing that changes here is that the contraction is done by torch.einsum on tensors that
stay resident on the device. The integrals h and V are uploaded once.

Benchmarked at our actual shapes (contracting a 4-index with a 6-index, and the 36-term
antisymmetrisation), complex128 on an RTX 3090 against numpy:

       M |  einsum x  |  antisym x
       8 |    442     |    39.5
      10 |     40     |    45.6
      12 |     70     |    45.8
      14 |     79     |    47.5
      16 |    100     |    46.9

A 6-index complex128 tensor is 0.27 GB at M=16 and 1.02 GB at M=20, so 24 GB reaches M=16-20
against a practical CPU ceiling of M=8.

Every function here is validated against the numpy path in tools/gpu_validate.py; agreement is
required at complex128 round-off before any result from this module is used.
"""
import itertools
import numpy as np
import torch

from qdiss.eom.hamiltonian import TWO_BODY_PREFACTOR

_PERMS3 = [(pu, pl) for pu in itertools.permutations(range(3))
           for pl in itertools.permutations(range(3))]
_SGN3 = {p: int(round(np.linalg.det(np.eye(3)[list(p)]))) for p in
         itertools.permutations(range(3))}


def device(prefer="cuda"):
    return torch.device(prefer if torch.cuda.is_available() and prefer == "cuda" else "cpu")


def up(a, dev):
    """Upload a numpy array as resident complex128 on the device."""
    return torch.as_tensor(np.ascontiguousarray(a), dtype=torch.complex128, device=dev)


def down(t):
    return t.detach().cpu().numpy()


def contract(terms, Hop, h_subscript, prefactor, tensors, out_subscript):
    """Contract a symbolic term list. `tensors` maps (name, rank) to device tensors."""
    tot = None
    for sign, facs in terms:
        subs = [h_subscript]
        ops = [Hop]
        for name, rank, u, l in facs:
            key = ("D", 1) if name == "D" else (name, rank)
            subs.append("".join(u) + "".join(l))
            ops.append(tensors[key])
        piece = torch.einsum(",".join(subs) + "->" + out_subscript, *ops)
        piece = piece * (sign * prefactor)
        tot = piece if tot is None else tot + piece
    return tot


def assemble(one_terms, two_terms, h, V, tensors, out_subscript):
    """i<[H, O]> for the H of qdiss.eom.hamiltonian. The prefactor is imported, never chosen."""
    one = contract(one_terms, h, "rs", 1.0, tensors, out_subscript)
    two = contract(two_terms, V, "rtsu", TWO_BODY_PREFACTOR, tensors, out_subscript)
    return 1j * (one + two)


def antisymmetrize3(T):
    out = torch.zeros_like(T)
    for pu, pl in _PERMS3:
        out = out + (_SGN3[pu] * _SGN3[pl]) * T.permute(
            tuple(pu) + tuple(3 + i for i in pl))
    return out / 36.0


def wedge(a, b):
    return torch.einsum("pq,vw->pvqw", a, b) - torch.einsum("pw,vq->pvqw", a, b)


def muller2405_gamma3(g1, g2, N, M, contraction_correct=True, orbB=None):
    """Mueller Eq. 14 in the closed-system limit, optionally Lackner-corrected.

    orbB must be the precompiled rank-3 gamma-wedge-lambda_2 orbit term list, used for the
    partial-trace inverse; see qdiss.muller2405.
    """
    if N <= 2:
        return torch.zeros((M,) * 6, dtype=torch.complex128, device=g1.device)
    prod = (torch.einsum("pvqw,xy->pvxqwy", g2, g1)
            + torch.einsum("pxqy,vw->pvxqwy", g2, g1)
            + torch.einsum("pq,vxwy->pvxqwy", g1, g2))
    trip = torch.einsum("pq,vw,xy->pvxqwy", g1, g1, g1)
    F = (N - 2) / N * prod - 2.0 * (N - 1) * (N - 2) / N ** 2 * trip
    F = antisymmetrize3(F)
    want = N * (N - 1) * (N - 2)
    got = torch.einsum("pvxpvx->", F).real
    if abs(float(got)) > 1e-14:
        F = F * (want / got)
    if not contraction_correct or orbB is None:
        return F
    delta = g2 - torch.einsum("pvxqwx->pvqw", F) / (N - 2)
    return F + trace_inverse((N - 2) * delta, M, orbB)


def trace_inverse(R, M, orbB):
    """A 6-index antisymmetric tensor whose partial trace over the third pair is R.

    Uses sum_y B[I,X] = (M-4) X + (I^S_X) + (S_X^I), the same relation as Lackner's Eq. 43,
    whose (d_1 - 4) denominator this reproduces.
    """
    dev = R.device
    I = torch.eye(M, dtype=torch.complex128, device=dev)
    S_R = torch.einsum("pxqx->pq", R)
    trS = torch.einsum("pp->", S_R) / (3 * M - 6)
    S_X = (S_R - trS * I) / (2 * M - 6)
    X = (R - wedge(I, S_X) - wedge(S_X, I)) / (M - 4)
    out = None
    for sign, sub, keys in orbB:
        piece = sign * torch.einsum(sub, *[{("G", 1): I, ("L", 2): X}[k] for k in keys])
        out = piece if out is None else out + piece
    return out
