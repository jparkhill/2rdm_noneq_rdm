#!/usr/bin/env python3
"""Is the 3.0 real, or a normalization artifact?

Single pair:  a+_p a+_v a_v a_p = n_p n_v, eigenvalues 0 or 1, so the DIAGONAL
element of gamma_2 is at most 1, and at most 2 in the factor-of-2 pair convention
used by dqg.antisym_matrix.

Eigenvalue: the eigenvector is a SUPERPOSITION |g> = sum_{p<v} c_pv |pv>, whose
operator is G+G with G = sum c_pv a_v a_p. Nothing bounds <G+G> by the single-pair
value once the c are coherent -- that is off-diagonal long-range order.

Checks, on the (8,4) geminal-power condensate:
  1. largest DIAGONAL element of D, which should respect the single-pair bound
  2. largest EIGENVALUE of D
  3. <G+G> built explicitly from the top eigenvector, compared against both,
     so the factor of 2 is pinned down rather than assumed
"""
import numpy as np

M, N = 8, 4
dim = 1 << M
PAIRS = [(p, v) for p in range(M) for v in range(p + 1, M)]
NP = len(PAIRS)


def ad(p, v):
    out = np.zeros_like(v); bit = 1 << p
    for s in np.nonzero(v)[0]:
        if s & bit: continue
        out[s | bit] += ((-1.0) ** bin(s & (bit - 1)).count("1")) * v[s]
    return out


def aa(p, v):
    out = np.zeros_like(v); bit = 1 << p
    for s in np.nonzero(v)[0]:
        if not (s & bit): continue
        out[s ^ bit] += ((-1.0) ** bin(s & (bit - 1)).count("1")) * v[s]
    return out


vac = np.zeros(dim, dtype=complex); vac[0] = 1.0
npair = M // 2


def Gdag_cond(v):
    out = np.zeros_like(v)
    for k in range(npair):
        out += ad(2 * k, ad(2 * k + 1, v))
    return out


psi = vac
for _ in range(N // 2):
    psi = Gdag_cond(psi)
psi = psi / np.linalg.norm(psi)

# gamma_2 in the pair basis, with the project's factor of 2
kets = {b: aa(w, aa(q, psi)) for b, (q, w) in enumerate(PAIRS)}
D = np.zeros((NP, NP), dtype=complex)
for a, (p, v) in enumerate(PAIRS):
    left = aa(v, aa(p, psi))
    for b in range(NP):
        D[a, b] = 2.0 * np.vdot(left, kets[b])
D = 0.5 * (D + D.conj().T)

diag = np.real(np.diag(D))
w, U = np.linalg.eigh(D)
print("tr D = %.4f   (N(N-1) = %d)" % (np.trace(D).real, N * (N - 1)))
print()
print("largest DIAGONAL element of D            : %.6f" % diag.max())
print("  same divided by the convention factor 2: %.6f   (single-pair <n_p n_v>)"
      % (diag.max() / 2))
print("largest EIGENVALUE of D                  : %.6f" % w.max())
print()

c = U[:, -1]
print("top eigenvector support: %d pairs with |c| > 0.1" % int((np.abs(c) > 0.1).sum()))


def Gapply(coef, v):
    """G = sum_{p<v} conj(c_pv) a_v a_p"""
    out = np.zeros_like(v)
    for idx, (p, q) in enumerate(PAIRS):
        if abs(coef[idx]) < 1e-14:
            continue
        out += np.conj(coef[idx]) * aa(q, aa(p, v))
    return out


Gpsi = Gapply(c, psi)
GdG = float(np.real(np.vdot(Gpsi, Gpsi)))
print("<G+G> from the top eigenvector           : %.6f" % GdG)
print("2 * <G+G>                                : %.6f" % (2 * GdG))
print()
if abs(2 * GdG - w.max()) < 1e-8:
    print("so eigenvalue = 2 <G+G>, the convention factor is confirmed, and")
    print("<G+G> = %.4f EXCEEDS the single-pair bound of 1." % GdG)
    print("A coherent superposition of pairs is not limited by n_p n_v <= 1;")
    print("that is condensation, and it is why the eigenvalue passes 2.")
else:
    print("MISMATCH: eigenvalue %.6f vs 2<G+G> %.6f -- convention not as assumed"
          % (w.max(), 2 * GdG))
