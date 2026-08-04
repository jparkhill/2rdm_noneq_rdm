"""Shared contraction of a symbolic term list. Infrastructure, not a derivation.

A term is (sign, factors) with factor = (name, rank, upper, lower) and name one of
'D' delta, 'G' gamma, 'L' lambda. This module only performs einsums; it makes no
statement about which terms are correct. Both derivations use it, so a fault here would
affect them equally and could not be caught by comparing them to each other -- which is
why both are held against fullci.py instead.
"""
import numpy as np
from .hamiltonian import TWO_BODY_PREFACTOR


def contract(terms, Hop, h_subscript, prefactor, tensors, out_subscript, M):
    tot = np.zeros((M,) * len(out_subscript), dtype=complex)
    for sign, facs in terms:
        subs = [h_subscript]
        ops = [Hop]
        for name, rank, up, lo in facs:
            key = ("D", 1) if name == "D" else (name, rank)
            subs.append("".join(up) + "".join(lo))
            ops.append(tensors[key])
        tot = tot + sign * prefactor * np.einsum(
            ",".join(subs) + "->" + out_subscript, *ops, optimize=True)
    return tot


def assemble(one_terms, two_terms, h, V, tensors, out_subscript, M):
    """i<[H, O]> for the H of hamiltonian.py. The prefactor is imported, never chosen."""
    one = contract(one_terms, h, "rs", 1.0, tensors, out_subscript, M)
    two = contract(two_terms, V, "rtsu", TWO_BODY_PREFACTOR, tensors, out_subscript, M)
    return 1j * (one + two)
