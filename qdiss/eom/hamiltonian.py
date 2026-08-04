"""The single definition of the Hamiltonian, imported by every derivation and by FullCI.

    H = sum_rs h^r_s a+_r a_s  +  TWO_BODY_PREFACTOR * sum_rtsu V^rt_su a+_r a+_t a_u a_s

with V antisymmetric in (r,t) and in (s,u) and Hermitian under (rt) <-> (su).

The prefactor is written once, here. Every derivation and the Fock-space builder import it,
so no derivation can be brought into agreement with another by adjusting a factor: they are
each held against FullCI, which is the only arbiter. If a derivation disagrees with FullCI
the derivation is wrong, and the remedy is never to change this constant.
"""
import numpy as np

TWO_BODY_PREFACTOR = 0.25


def random_h(M, rng):
    A = rng.normal(size=(M, M)) + 1j * rng.normal(size=(M, M))
    return 0.25 * (A + A.conj().T)


def random_V(M, rng):
    B = rng.normal(size=(M,) * 4) + 1j * rng.normal(size=(M,) * 4)
    V = B - np.transpose(B, (1, 0, 2, 3))
    V = V - np.transpose(V, (0, 1, 3, 2))
    return 0.25 * (V + np.transpose(V.conj(), (2, 3, 0, 1)))


def check_symmetries(V):
    """Return the three symmetry residuals V must satisfy for the convention to hold."""
    return {
        "antisym upper": float(np.abs(V + np.transpose(V, (1, 0, 2, 3))).max()),
        "antisym lower": float(np.abs(V + np.transpose(V, (0, 1, 3, 2))).max()),
        "hermitian": float(np.abs(V - np.transpose(V.conj(), (2, 3, 0, 1))).max()),
    }
