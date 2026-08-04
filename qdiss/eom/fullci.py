"""Ground truth. No symbolic derivation of any kind.

Builds H as a matrix on the full Fock space directly from the convention in
hamiltonian.py, and computes the RDMs and their time derivatives by exact linear algebra:

    gamma_k[P;Q]     = <psi| a+_P a_Q |psi>
    gamma_k-dot[P;Q] = i <psi| [H, a+_P a_Q] |psi>

Nothing here consumes a term list, so agreement between a derived EOM and this module is
independent evidence about the derivation. Disagreement means the derivation is wrong.
"""
import itertools
import numpy as np
from .hamiltonian import TWO_BODY_PREFACTOR


def ad(p, v):
    out = np.zeros_like(v)
    bit = 1 << p
    for s in np.nonzero(v)[0]:
        if s & bit:
            continue
        out[s | bit] += ((-1.0) ** bin(s & (bit - 1)).count("1")) * v[s]
    return out


def aa(p, v):
    out = np.zeros_like(v)
    bit = 1 << p
    for s in np.nonzero(v)[0]:
        if not (s & bit):
            continue
        out[s ^ bit] += ((-1.0) ** bin(s & (bit - 1)).count("1")) * v[s]
    return out


def random_state(M, N, rng):
    dim = 1 << M
    popc = np.array([bin(s).count("1") for s in range(dim)])
    psi = np.where(popc == N, rng.normal(size=dim) + 1j * rng.normal(size=dim), 0.0)
    return psi / np.linalg.norm(psi)


def build_H(h, V, M):
    """H as a dense matrix, from the convention in hamiltonian.py and nothing else."""
    dim = 1 << M
    H = np.zeros((dim, dim), dtype=complex)
    basis = np.eye(dim, dtype=complex)
    for p in range(M):
        for q in range(M):
            if abs(h[p, q]) < 1e-14:
                continue
            for k in range(dim):
                H[:, k] += h[p, q] * ad(p, aa(q, basis[:, k]))
    for p, t, q, u in itertools.product(*([range(M)] * 4)):
        if abs(V[p, t, q, u]) < 1e-14:
            continue
        for k in range(dim):
            H[:, k] += TWO_BODY_PREFACTOR * V[p, t, q, u] * ad(p, ad(t, aa(u, aa(q, basis[:, k]))))
    return H


def _ladder(psi, M, k):
    """Rows indexed by k-tuples P, each row a_{P_k} ... a_{P_1} |psi>."""
    out = []
    for idx in itertools.product(*([range(M)] * k)):
        v = psi
        for p in idx:
            v = aa(p, v)
        out.append(v)
    return np.array(out)


def rdm(psi, M, k):
    L = _ladder(psi, M, k)
    return (L.conj() @ L.T).reshape((M,) * (2 * k))


def rdm_dot(psi, H, M, k):
    """i <[H, a+_P a_Q]> by exact linear algebra."""
    Hp = H @ psi
    L = _ladder(psi, M, k)
    Lh = _ladder(Hp, M, k)
    fwd = Lh.conj() @ L.T
    rev = L.conj() @ Lh.T
    return (1j * (fwd - rev)).reshape((M,) * (2 * k))
