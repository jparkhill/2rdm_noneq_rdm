#!/usr/bin/env python3
"""What extensivity demands of the transfer matrix G: the ladder operators must be fragment-local.

THE DERIVATION. For two non-interacting fragments in a product state, the k-annihilated states
factorise. Splitting a k-tuple I into I_A and I_B with j = |I_A|,

    u_I = +- (A_{I_A}|Psi_A>) (x) (A_{I_B}|Psi_B>)

so range(R) = direct sum over j of V_A^{(j)} (x) V_B^{(k-j)}, where V_f^{(j)} is spanned by the
j-fold annihilations of fragment f. Since annihilators map V_f^{(j)} INTO V_f^{(j+1)} exactly, the
compressed ladder operators

    frak_a_x = (a_x restricted to range(R)) : level k -> level k+1,        G_xy = frak_a_x^dag frak_a_y

come out FACTORISED:

    x in A :  frak_a_x = alpha^A_x (x) I_B                (raises j by 1)
    y in B :  frak_a_y = Pi_A (x) alpha^B_y               (leaves j alone), Pi_A = (-1)^{N_A} parity

Same-fragment blocks of G are therefore (ladder^dag ladder on that fragment) tensor identity on the
other; cross-fragment blocks are (ladder^dag on one) tensor (ladder on the other). So:

    EXTENSIVITY  <=>  G_xy = frak_a_x^dag frak_a_y with every frak_a_x FRAGMENT-LOCAL,
                      i.e. acting as a ladder on the factor containing x and as identity (up to the
                      parity string) on the other.

This is the same statement as "connected contraction" for cumulants and "fragment-local jump
operator" for GKSL, in a third language. Note that the ladder form ALSO delivers G >= 0 for free,
since G is then a Gram matrix of the frak_a_x, and `sum_x G_xx = sum_x frak_a_x^dag frak_a_x` is the
modelled number operator, so the contraction identity is one condition on the ladders.

THE TEST. On a composite of H3 bent and H2 as an exact product state, with k = 2 so the reconstructed
object is gamma_3, compare the cumulant leak of:

    exact G            local ladders                    control, must be 0
    G' = V^dag G V     V a random unitary on range(R)    delocalises the ladders
    G' = V^dag G V     V random but BLOCK-DIAGONAL in j  keeps the grading, breaks the tensor
                                                         factorisation inside each block
    G = D (x) I_r      D the hole 1-RDM, tr D = N-k      the crudest ansatz

Every one of these satisfies `G >= 0` and `sum_x G_xx = (N-k) I_r` (a unitary conjugation preserves
both), so the ONLY thing distinguishing them is locality. `leak` is
max |lambda_3| over index tuples straddling the fragments, divided by max |lambda_3| inside a
fragment; the exact cumulant has leak zero, verified at 8.9e-16 in tools/ansatz_extensivity.py.

Symbols: gamma_k[p_1..p_k, q_1..q_k] = <a+_p1..a+_pk a_qk..a_q1>; Mat_k(T)[I,J] = k! T[I+J] over
sorted k-tuples; Mat_{k+1}[(I,x),(J,y)] = (k+1) R_c(:,I)^dag G_xy R_c(:,J) as verified in
tools/G_transfer.py. Ground truth qdiss.eom.fullci.
"""
import sys
from math import factorial
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from molecular_2rdm import build_H_sector, molecule, sector
from qdiss import ansatz as az
from qdiss.eom import fullci

H3_BENT = "H 0.00 0.00 0.0; H 1.00 0.00 0.0; H 0.35 0.90 0.0"
H2 = "H 0 0 0; H 0 0 0.9"
MA, NA = 6, 3
MB, NB = 4, 2
M, N, K = MA + MB, NA + NB, 2


def frag_state(atom, spin, Mx, Nx):
    mol, mf, Mchk, h, V, enuc = molecule(atom, "sto-3g", 0, spin)
    assert Mchk == Mx
    dets = sector(Mx, Nx)
    es, vec = np.linalg.eigh(build_H_sector(h, V, Mx, dets))
    p = np.zeros(1 << Mx, dtype=complex)
    p[dets] = vec[:, 0]
    return p


psiA = frag_state(H3_BENT, 1, MA, NA)
psiB = frag_state(H2, 0, MB, NB)
psi = np.zeros(1 << M, dtype=complex)
for dA in np.nonzero(psiA)[0]:
    for dB in np.nonzero(psiB)[0]:
        psi[dA | (dB << MA)] = psiA[dA] * psiB[dB]
print("composite H3 bent (M=6,N=3) x H2 (M=4,N=2) -> M=%d N=%d, k=%d, norm %.12f"
      % (M, N, K, float(np.linalg.norm(psi))))

TUP_K, TUP_K1 = az.tuples(M, K), az.tuples(M, K + 1)
SEC_HI, SEC_LO = sector(M, N - K), sector(M, N - K - 1)
posK = {t: i for i, t in enumerate(TUP_K)}


def annihilate(p, tup):
    out = p
    for t in tup:
        out = fullci.aa(t, out)
    return out


g1 = fullci.rdm(psi, M, 1)
l2 = az.cumulant(2, fullci.rdm(psi, M, 2), {("G", 1): g1})
l3_exact = az.cumulant(3, fullci.rdm(psi, M, 3), {("G", 1): g1, ("L", 2): l2})
SEP3 = az.separable(3, {("G", 1): g1, ("L", 2): l2})

U = np.zeros((len(SEC_HI), len(TUP_K)), dtype=complex)
for a, I in enumerate(TUP_K):
    U[:, a] = annihilate(psi, I)[SEC_HI]
R = np.sqrt(float(factorial(K))) * U
Q, S, Vh = np.linalg.svd(R, full_matrices=False)
keep = S > 1e-10 * S[0]
r = int(keep.sum())
P = Q[:, keep]
Rc = np.diag(S[keep]) @ Vh[keep, :]
print("  (N-k)-sector dim %d, r = rank(Mat_2) = %d, G is %d x %d" % (len(SEC_HI), r, M * r, M * r))

# exact compressed ladders frak_a_x = a_x P, as (N-k-1)-sector x r matrices
lad = np.zeros((M, len(SEC_LO), r), dtype=complex)
for x in range(M):
    for al in range(r):
        b = np.zeros(1 << M, dtype=complex)
        b[SEC_HI] = P[:, al]
        lad[x, :, al] = fullci.aa(x, b)[SEC_LO]


def G_from_ladders(L):
    return np.array([[L[x].conj().T @ L[y] for y in range(M)] for x in range(M)])


def mat3_from_G(Gb):
    Y = np.zeros((len(TUP_K1),) * 2, dtype=complex)
    for a, L in enumerate(TUP_K1):
        I, x = L[:K], L[K]
        for b, Lp in enumerate(TUP_K1):
            J, y = Lp[:K], Lp[K]
            Y[a, b] = (K + 1) * (Rc[:, posK[I]].conj() @ (Gb[x, y] @ Rc[:, posK[J]]))
    return Y


def leak_of(Y):
    g3 = az.coord_to_tensor(Y / factorial(K + 1), M, K + 1)
    lam = g3 - SEP3
    A = np.abs(lam)
    ins = float(A[tuple(slice(0, MA) for _ in range(2 * (K + 1)))].max())
    keepm = np.zeros((M,) * (2 * (K + 1)), dtype=bool)
    for lo, hi in ((0, MA), (MA, M)):
        keepm[tuple(slice(lo, hi) for _ in range(2 * (K + 1)))] = True
    out = float(A[~keepm].max())
    return ins, out / max(ins, 1e-30), float(np.abs(lam - l3_exact).max())


# the j-grading of range(R): eigenvalues of the A-particle-number operator restricted to range(R)
NA_op = np.zeros((len(SEC_HI),) * 2, dtype=complex)
for i, d in enumerate(SEC_HI):
    NA_op[i, i] = bin(d & ((1 << MA) - 1)).count("1")
jgrade = np.round(np.real(np.diag(P.conj().T @ NA_op @ P))).astype(int)
Pj = P.conj().T @ NA_op @ P
ev_j = np.round(np.real(np.linalg.eigvalsh(0.5 * (Pj + Pj.conj().T)))).astype(int)
print("  A-particle numbers on range(R): %s" % np.array2string(np.sort(ev_j)))

rng = np.random.default_rng(5)


def rand_unitary(n):
    A = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    Qm, _ = np.linalg.qr(A)
    return Qm


Gex = G_from_ladders(lad)
print()
print("  %-40s %11s %11s %11s %11s"
      % ("G ansatz", "|inside|", "leak", "|G_xx-(N-k)I|", "min eig G"))


def report(name, Gb):
    ins, lk, err = leak_of(mat3_from_G(Gb))
    Sg = sum(Gb[x, x] for x in range(M))
    big = np.block([[Gb[x, y] for y in range(M)] for x in range(M)])
    mn = np.linalg.eigvalsh(0.5 * (big + big.conj().T)).min()
    print("  %-40s %11.3e %11.3e %11.3e %11.3e"
          % (name, ins, lk, np.abs(Sg - (N - K) * np.eye(r)).max(), mn))


report("exact G (local ladders)", Gex)

Vr = rand_unitary(r)
report("V^dag G V, V random on range(R)", np.array([[Vr.conj().T @ Gex[x, y] @ Vr
                                                    for y in range(M)] for x in range(M)]))

# a unitary that preserves the j-grading but mixes freely inside each graded block
w, Uj = np.linalg.eigh(0.5 * (Pj + Pj.conj().T))
lab = np.round(np.real(w)).astype(int)
Vb = np.zeros((r, r), dtype=complex)
for j in sorted(set(lab)):
    idx = np.where(lab == j)[0]
    Vb[np.ix_(idx, idx)] = rand_unitary(len(idx))
Vb = Uj @ Vb @ Uj.conj().T
report("V^dag G V, V random but j-block diagonal", np.array([[Vb.conj().T @ Gex[x, y] @ Vb
                                                             for y in range(M)] for x in range(M)]))

D = np.array([[np.trace(Gex[x, y]) / r for y in range(M)] for x in range(M)])
report("G = D (x) I_r, D the hole 1-RDM", np.array([[D[x, y] * np.eye(r)
                                                    for y in range(M)] for x in range(M)]))
print()
print("  All four satisfy G >= 0 and sum_x G_xx = (N-k) I_r; only LOCALITY of the ladders differs.")
