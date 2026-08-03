#!/usr/bin/env python3
"""Coherent 2-RDM dynamics, contracted straight out of the verified term lists.

The path is derivation -> code -> simulation with no hand transcription: the
einsum for every term is built from the symbolic factors themselves, so a term
list that passed the operator-identity check is the thing that propagates.

Conventions, fixed once and checked against explicit Fock space:
    a^P_Q      = a+_{p1}..a+_{pk} a_{qk}..a_{q1}
    gamma2[p,v,q,w] = <a^{pv}_{qw}>
    H          = sum_rs h[r,s] a^r_s + (1/2) sum_rtsu V[r,t,s,u] a^{rt}_{su}
    i d gamma2 / dt = <[H, a^{pv}_{qw}]>
"""
import numpy as np

OUT = "pvqw"


def einsum_for(term, ham_indices):
    """(sign, factors) -> (sign, subscript string, list of operand keys).

    Operand keys are ('H',), ('G', rank) or ('I',) for a delta, in the order the
    subscripts expect. Nothing about the physics is hard-coded here.
    """
    sign, facs = term
    subs, ops = [], []
    subs.append("".join(ham_indices))
    ops.append(("H",))
    for name, rank, up, lo in facs:
        if name == "D":
            subs.append(up[0] + lo[0])
            ops.append(("I",))
        elif name == "G":
            subs.append("".join(up) + "".join(lo))
            ops.append(("G", rank))
        else:
            raise ValueError("unexpected factor %r" % (name,))
    return sign, ",".join(subs) + "->" + OUT, ops


def compile_terms(terms, ham_indices):
    return [einsum_for(t, ham_indices) for t in terms]


def evaluate(compiled, H, gammas, M, dtype=complex):
    """Sum the compiled einsums. gammas maps rank -> tensor."""
    eye = np.eye(M, dtype=dtype)
    out = np.zeros((M,) * 4, dtype=dtype)
    for sign, sub, ops in compiled:
        arrays = []
        for key in ops:
            if key[0] == "H":
                arrays.append(H)
            elif key[0] == "I":
                arrays.append(eye)
            else:
                arrays.append(gammas[key[1]])
        out = out + sign * np.einsum(sub, *arrays, optimize=True)
    return out


class Coherent2:
    """i dgamma2/dt from a one-body h and a two-body V, with a gamma_3 closure."""

    def __init__(self, one_terms, two_terms, h, V, closure):
        self.M = h.shape[0]
        self.h = h
        self.V = V
        self.closure = closure
        self.one = compile_terms(one_terms, "rs")
        self.two = compile_terms(two_terms, "rtsu")

    def ddt(self, g2):
        g3 = self.closure(g2)
        gam = {2: g2, 3: g3}
        rhs = evaluate(self.one, self.h, gam, self.M)
        rhs = rhs + 0.5 * evaluate(self.two, self.V, gam, self.M)
        return 1j * rhs


def closure_gamma3_zero(M):
    z = np.zeros((M,) * 6, dtype=complex)
    return lambda g2: z


def antisymmetrize(g2):
    """Enforce the two-particle antisymmetry gamma2[pvqw] = -gamma2[vpqw] etc."""
    a = g2 - np.transpose(g2, (1, 0, 2, 3))
    a = a - np.transpose(a, (0, 1, 3, 2))
    return a / 4.0


def rk4(f, y, dt):
    k1 = f(y)
    k2 = f(y + 0.5 * dt * k1)
    k3 = f(y + 0.5 * dt * k2)
    k4 = f(y + dt * k3)
    return y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
