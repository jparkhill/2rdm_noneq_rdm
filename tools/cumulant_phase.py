#!/usr/bin/env python3
"""Why does even the mp2 reconstruction leave occupations frozen?

mp2 is NOT diagonal in the natural-orbital pair basis, so the earlier no-go argument does not
apply to it, yet max|ndot| is 4e-17. Hypothesis: the amplitude T = V * (real occupation and
energy factors) is IN PHASE with V, so the occupation channel, which is an imaginary part,
cancels. If so, the anti-Hermitian-in-phase combination i(T - T-dagger), which is equally
Hermitian as a pair matrix but out of phase with V, should move occupations.

That would make the no-go broader than stated: not "lambda_2 must depend on V" but "lambda_2
must have a component out of phase with V", which is a retardation effect no instantaneous
reconstruction can supply.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from physical_hamiltonian import spin_orbital_H
from qdiss import closures
from qdiss.eom import fullci, gwt
from qdiss.eom.evaluate import assemble

CUM1 = gwt.terms_cumulant(1)
mol, mf, M, h, V, enuc = spin_orbital_H("H 0 0 0; H 0 0 0.74", "6-31g", 0)
Hm = fullci.build_H(h, V, M)
popc = np.array([bin(s).count("1") for s in range(1 << M)])
sub = np.where(popc == 2)[0]
es, vec = np.linalg.eigh(Hm[np.ix_(sub, sub)])
psi0 = np.zeros(1 << M, dtype=complex)
psi0[sub] = vec[:, 0]
rr = np.random.default_rng(5)
kick = np.zeros(1 << M, dtype=complex)
kick[sub] = rr.normal(size=len(sub)) + 1j * rr.normal(size=len(sub))
psik = psi0 + 0.3 * kick / np.linalg.norm(kick)
psik /= np.linalg.norm(psik)
g1 = fullci.rdm(psik, M, 1)
l2ex = fullci.rdm(psik, M, 2) - closures.wedge(g1, g1)
n, U = np.linalg.eigh(0.5 * (g1 + g1.conj().T))
n = np.clip(n, 1e-10, 1.0)


def amplitude():
    F = h + np.einsum("prqs,rs->pq", V, g1)
    e = np.real(np.diag(U.conj().T @ F @ U))
    Vn = np.einsum("ap,bv,abcd,cq,dw->pvqw", U.conj(), U.conj(), V, U, U, optimize=True)
    occ = np.einsum("q,w,p,v->pvqw", n, n, 1 - n, 1 - n)
    den = (e[None, None, :, None] + e[None, None, None, :]
           - e[:, None, None, None] - e[None, :, None, None])
    den = np.where(np.abs(den) < 1e-6, np.inf, den)
    T = Vn * occ / den
    T = T - np.transpose(T, (1, 0, 2, 3))
    return T - np.transpose(T, (0, 1, 3, 2))


def back(L):
    return np.einsum("ap,bv,pvqw,cq,dw->abcd", U, U, L, U.conj(), U.conj(), optimize=True)


def ndot(l2):
    t = {("D", 1): np.eye(M, dtype=complex), ("G", 1): g1, ("L", 2): l2}
    d = assemble(CUM1[0], CUM1[1], h, V, t, "pq", M)
    return np.abs(np.real(np.diag(U.conj().T @ d @ U))).max()


T = amplitude()
Th = np.transpose(T.conj(), (2, 3, 0, 1))
cases = (("exact lambda_2", l2ex),
         ("mp2  0.25(T + T^H)", back(0.25 * (T + Th))),
         ("out of phase 0.25 i(T - T^H)", back(0.25j * (T - Th))),
         ("hf (zero)", np.zeros((M,) * 4, dtype=complex)))
print("max|ndot| in the natural-orbital basis, kicked complex state")
for name, l2 in cases:
    herm = np.abs(l2 - np.transpose(l2.conj(), (2, 3, 0, 1))).max()
    print("   %-30s |l2| %8.5f  hermiticity %8.1e  max|ndot| %11.3e"
          % (name, np.abs(l2).max(), herm, ndot(l2)))
print()
ip = back(0.25 * (T + Th))
op = back(0.25j * (T - Th))
num = float(np.real(np.vdot(l2ex.ravel(), op.ravel())))
print("overlap of the exact cumulant with the out-of-phase component : %.5f" % num)
print("overlap of the exact cumulant with the in-phase component     : %.5f"
      % float(np.real(np.vdot(l2ex.ravel(), ip.ravel()))))

Ipair = np.eye(M, dtype=complex)


def tracefree(W):
    S = np.einsum("pxqx->pq", W)
    return (W - (1.0 / (M - 2)) * closures.sym_wedge(S, Ipair)
            + (np.trace(S) / ((M - 2) * (M - 1))) * closures.wedge(Ipair, Ipair))


print()
print("raw V perturbations, which are pair-Hermitian and therefore in phase")
Vas = V - np.transpose(V, (1, 0, 2, 3))
for lab, W in (("trace-free V", tracefree(Vas)),
               ("trace-free i*V", tracefree(1j * Vas)),
               ("V, no trace correction", Vas)):
    hz = np.abs(W - np.transpose(W.conj(), (2, 3, 0, 1))).max()
    print("   %-24s |W| %8.4f  hermiticity %8.1e  max|ndot| %11.3e"
          % (lab, np.abs(W).max(), hz, ndot(0.05 * W)))
print()
print("if a real pair-Hermitian V perturbation gives ndot ~ 0 here, the earlier random-H")
print("result at 4.935 came from a COMPLEX V, where V is Hermitian but not real symmetric,")
print("so it carries an out-of-phase part that a real molecular V cannot.")
