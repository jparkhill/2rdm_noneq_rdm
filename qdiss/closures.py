"""Cumulant reconstructions lambda_2[gamma_1] for the MKNO 1-RDM equation of motion.

These are CUMULANT reconstructions, not RDMFT. Nothing here minimises an energy functional or
builds a coupling matrix; the RDMFT functionals named below supply only the algebraic shape of
the cumulant. The object produced is lambda_2, to be substituted into the unapproximated MKNO
rank-1 EOM from qdiss.eom.gwt, which is a function of gamma_1 and lambda_2 alone and is
validated against qdiss.eom.fullci at 1e-15 by qdiss.eom.validate.

Structure. Every reconstruction except mp2 is diagonal in the natural-orbital pair basis and
is fixed by a single symmetric coefficient c(n_p, n_v):

    lambda_2[pv,qw] = (1/2) c(n_p,n_v) ( d_pw d_vq - d_pq d_vw )

which is antisymmetric in both index pairs because c is symmetric, and reproduces the same
energy as the non-antisymmetrised form since V is antisymmetric. Writing X for the exchange
weight that replaces the Hartree-Fock n_p n_v, every such c is c = n_p n_v - X.

    hf        X = n_p n_v, so lambda_2 = 0. Substituted into the rank-1 EOM this gives exactly
              TDHF, i[(h+W)^T, gamma_1] with W[p,q] = V[p,r,q,s] gamma_1[r,s], verified at
              3.9e-16.
    muller_rdmft  X = sqrt(n_p n_v), the RDMFT Mueller / Buijse-Baerends exchange shape. NOT to
              be confused with qdiss.muller2405, the 3-RDM reconstruction of Mueller et al.
              arXiv:2405.05093 -- different authors, different object, different rank.
    power     X = (n_p n_v)^alpha, alpha ~ 0.55-0.65. Sharma et al. shape.
    bbc       Mueller patched by occupation block: c = 0 on the diagonal p = v, so the cumulant
              carries no self-pair term (Goedecker-Umrigar); and the sign of the square root is
              flipped when both orbitals are weakly occupied, removing Mueller's spurious
              enhanced exchange among virtual-like orbitals (BBC1). The content in cumulant
              language is a sign change between blocks, c < 0 for occupied and mixed pairs and
              c > 0 for virtual-virtual. This is BBC1 plus GU, NOT BBC3: it lacks BBC3's
              frontier-pair treatment and BBC2's strong-strong restoration.
    blocked   ours. lambda_2 = A^f + f^A with A = gamma_1 - gamma_1^2 and f = g(gamma_1) fixed
              by demanding the partial-trace identity exactly. In the natural-orbital basis,
              with a_i = n_i(1-n_i), T = sum_j a_j and Q = sum_j a_j/(T - 2 a_j),
                  f_i = -a_i / [ (1+Q) (T - 2 a_i) ]
              Both factors are proportional to a_i, so each vanishes as an orbital empties or
              fills. It is the only member of this family satisfying the identity by
              construction.

and one that is not occupation-diagonal, which is the whole point of it:

    mp2       the first-order amplitude,
                  T[pv,qw] = V[pv,qw] n_q n_w (1-n_p)(1-n_v) / (e_q + e_w - e_p - e_v)
              with e the diagonal Fock energies in the natural-orbital basis, antisymmetrised
              and Hermitised. Because it depends on V it has coherences in that basis, and by
              the no-go in tools/lambda2_nogo.py it is the only candidate here that can move an
              occupation at all: any lambda_2 that is a matrix function of gamma_1 commutes with
              gamma_1, is diagonal in the natural-orbital pair basis, and contributes nothing to
              the diagonal of gamma_1-dot. That is the known TD-RDMFT statement that adiabatic
              occupation-number functionals leave occupations frozen.

Partial-trace identity. Every physical lambda_2 satisfies

    sum_x lambda_2[px,qx] = -(gamma_1 - gamma_1^2) = -gamma_1 eta_1

so the contraction of the cumulant IS the particle-hole blocking product, and since
gamma_1 eta_1 >= 0 exactly when 0 <= gamma_1 <= 1, the identity is the Pauli bound written in
cumulant variables. `trace_correct` enforces it on any reconstruction: for antisymmetric W with
partial trace S,

    Corr(S) = (1/(M-2)) (S^I + I^S) - [tr S/((M-2)(M-1))] (I^I)

has partial trace exactly S, so W - Corr(S_W) + blocked satisfies the identity while keeping
every coherence of W. Only the trace channel is replaced. Note that the blocked structure is
one particular solution among an affine family; the choice is not neutral, and on H2 it removes
a spurious near-zero mode from mp2 while destroying Mueller's 0.86 pole.
"""
import numpy as np

NAMES = ("hf", "blocked", "muller_rdmft", "power", "bbc", "mp2")


def wedge(a, b):
    return np.einsum("pq,vw->pvqw", a, b) - np.einsum("pw,vq->pvqw", a, b)


def sym_wedge(a, b):
    return wedge(a, b) + wedge(b, a)


def _from_c(c, U, M):
    D = np.zeros((M,) * 4, dtype=complex)
    for p in range(M):
        for v in range(M):
            D[p, v, v, p] += 0.5 * c[p, v]
            D[p, v, p, v] -= 0.5 * c[p, v]
    return np.einsum("ap,bv,pvqw,cq,dw->abcd", U, U, D, U.conj(), U.conj(), optimize=True)


def coefficient(name, n, alpha=0.6):
    M = len(n)
    prod = np.outer(n, n)
    root = np.sqrt(np.clip(prod, 0.0, None))
    if name == "hf":
        return np.zeros((M, M))
    if name == "muller_rdmft":
        return prod - root
    if name == "power":
        return prod - np.clip(prod, 0.0, None) ** alpha
    if name == "bbc":
        X = root.copy()
        weak = np.outer(n < 0.5, n < 0.5)
        X = np.where(weak & ~np.eye(M, dtype=bool), -X, X)
        np.fill_diagonal(X, np.diag(prod))
        return prod - X
    if name == "blocked":
        a = n - n * n
        T = a.sum()
        den = T - 2.0 * a
        if np.abs(den).min() < 1e-12:
            return np.zeros((M, M))
        Q = (a / den).sum()
        f = -a / ((1.0 + Q) * den)
        return -2.0 * (np.outer(a, f) + np.outer(f, a))
    raise ValueError(name)


def _mp2(g1, h, V, M, n, U):
    F = h + np.einsum("prqs,rs->pq", V, g1)
    e = np.real(np.diag(U.conj().T @ F @ U))
    Vn = np.einsum("ap,bv,abcd,cq,dw->pvqw", U.conj(), U.conj(), V, U, U, optimize=True)
    occ = np.einsum("q,w,p,v->pvqw", n, n, 1.0 - n, 1.0 - n)
    den = (e[None, None, :, None] + e[None, None, None, :]
           - e[:, None, None, None] - e[None, :, None, None])
    den = np.where(np.abs(den) < 1e-6, np.inf, den)
    T = Vn * occ / den
    T = T - np.transpose(T, (1, 0, 2, 3))
    T = T - np.transpose(T, (0, 1, 3, 2))
    L = 0.25 * (T + np.transpose(T.conj(), (2, 3, 0, 1)))
    return np.einsum("ap,bv,pvqw,cq,dw->abcd", U, U, L, U.conj(), U.conj(), optimize=True)


def lambda2(g1, name, h=None, V=None, alpha=0.6):
    """Reconstruct the 2-cumulant from gamma_1. h and V are needed only for mp2."""
    M = g1.shape[0]
    if name == "hf":
        return np.zeros((M,) * 4, dtype=complex)
    n, U = np.linalg.eigh(0.5 * (g1 + g1.conj().T))
    n = np.clip(n, 1e-10, 1.0)
    if name == "mp2":
        if h is None or V is None:
            raise ValueError("mp2 needs h and V")
        return _mp2(g1, h, V, M, n, U)
    return _from_c(coefficient(name, n, alpha), U, M)


def trace_target(g1):
    """The value sum_x lambda_2[px,qx] must take: -(gamma_1 - gamma_1^2)."""
    return -(g1 - g1 @ g1)


def trace_residual(l2, g1):
    return float(np.abs(np.einsum("pxqx->pq", l2) - trace_target(g1)).max())


def trace_correct(l2, g1):
    """Enforce the partial-trace identity, keeping every coherence of l2."""
    M = g1.shape[0]
    I = np.eye(M, dtype=complex)
    S = np.einsum("pxqx->pq", l2)
    corr = ((1.0 / (M - 2)) * sym_wedge(S, I)
            - (np.trace(S) / ((M - 2) * (M - 1))) * wedge(I, I))
    return l2 - corr + lambda2(g1, "blocked")
