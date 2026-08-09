#!/usr/bin/env python3
"""Does the EXACT 1-RDM collision term admit a fermionic gain-loss (Lindblad) form?

THE QUESTION. The MKNO rank-1 equation of motion is exact in (gamma_1, lambda_2):

    gamma_1-dot = i[f, gamma_1] + C[lambda_2],    f[p,q] = h[p,q] + sum_rs V[p,r,q,s] gamma_1[r,s]

where f is the generalised Fock operator and C is the collision term, defined here as exactly
those terms of gwt.terms_cumulant(1) that carry a ("L",2) factor. The commutator is unitary and
cannot move the spectrum of gamma_1, so every violation of 0 <= gamma_1 <= 1 must come from C.

A fermionic Lindblad master equation with loss channels (jump operators a_p) and gain channels
(jump operators a_p^dagger) produces a 1-RDM equation of the form

    gamma_1-dot = i[f, gamma_1] + (1/2){Gam_plus, eta_1} - (1/2){Gam_minus, gamma_1}

with eta_1 = 1 - gamma_1 the hole density matrix, {,} the anticommutator, and Gam_plus,
Gam_minus both POSITIVE SEMIDEFINITE. Gam_plus is the gain (in-scattering) matrix and Gam_minus
the loss (out-scattering) matrix. That form makes the two-sided Pauli bound automatic: at an
empty orbital (gamma_1 |phi> = 0, so eta_1 |phi> = |phi>) the rate is <phi|Gam_plus|phi> >= 0, and
at a full orbital (gamma_1 |phi> = |phi>) it is -<phi|Gam_minus|phi> <= 0. Nothing is checked; the
sign follows from where the factors sit. This is the fermionic analogue of Lindblad form, and it
is what we would impose on a CLOSURE so that positivity is structural rather than monitored.

NOTE on an earlier error of mine: the sandwich form eta_1 K gamma_1 + gamma_1 K^dag eta_1 is NOT
the right target. It gives exactly ZERO rate at both boundaries, so empty orbitals could never
fill. The anticommutator gain-loss form above is the correct structure.

WHAT IS TESTED HERE. Two exact necessary conditions, which follow from the form by compressing
onto the extreme eigenspaces of gamma_1. Writing P0 for the projector onto the null space of
gamma_1 (exactly empty orbitals) and P1 for the projector onto its unit eigenspace (exactly full):

    P0 C P0 = P0 Gam_plus P0  >= 0        (the loss term is annihilated by gamma_1 |phi> = 0)
    P1 C P1 = -P1 Gam_minus P1 <= 0       (the gain term is annihilated by eta_1 |phi> = 0)

These are block generalisations of the Nagumo tangency condition. They must hold for the exact
collision term evaluated on any genuine quantum state, since exact dynamics preserves the Pauli
bound; failure would mean the collision term has been misidentified. If they hold, the gain-loss
form is consistent with the exact equation and becomes a candidate structural condition to impose
on approximate closures.

The test state is an exact FCI ground state in a RESTRICTED determinant space of H4/sto-3g chosen
so that both extreme eigenspaces are nonempty and exact: spin-orbitals 0 and 1 occupied in every
determinant (so n = 1 exactly) and 6, 7 empty in every determinant (so n = 0 exactly), with the
remaining 2 electrons correlated in the active spin-orbitals 2..5. It is a genuine N-electron
state, so its gamma_1 and lambda_2 are exactly N-representable; it is not an eigenstate of the
full Hamiltonian, which is irrelevant -- the Pauli bound must be respected by the flow at every
state, not only at eigenstates.
"""
import sys
import numpy as np
sys.path.insert(0, "."); sys.path.insert(0, "tools")
import stationary_lagrangian as S
from qdiss import closures
from qdiss.eom import fullci, gwt

M, N, h, V = S.M, S.N, S.h, S.V
EYE = np.eye(M, dtype=complex)
L2KEY = ("L", 2)

one, two = gwt.terms_cumulant(1)
COLL = ([], [])
REST = ([], [])
for which, lst in ((0, one), (1, two)):
    for s, f in lst:
        keys = [(n, r) for n, r, u, l in f]
        subs = ["".join(u) + "".join(l) for n, r, u, l in f]
        (COLL if L2KEY in keys else REST)[which].append((s, subs, keys))

print("rank-1 MKNO EOM, terms carrying a lambda_2 factor (the collision term C):")
print("   %d one-body + %d two-body carry lambda_2;  %d + %d do not"
      % (len(COLL[0]), len(COLL[1]), len(REST[0]), len(REST[1])))
for s, subs, keys in COLL[0] + COLL[1]:
    print("      sign %+d   subs %-22s keys %s" % (s, str(subs), keys))
print()


def build(terms, tab):
    tot = np.zeros((M, M), dtype=complex)
    for s, subs, keys in terms[0]:
        tot += s * np.einsum(",".join(["rs"] + subs) + "->pq", h,
                             *[tab[k] for k in keys], optimize=True)
    for s, subs, keys in terms[1]:
        tot += s * 0.25 * np.einsum(",".join(["rtsu"] + subs) + "->pq", V,
                                    *[tab[k] for k in keys], optimize=True)
    return 1j * tot


def collision(g1, l2):
    return build(COLL, {("D", 1): EYE, ("G", 1): g1, L2KEY: l2})


# ---- restricted-space exact state: n = 1 exactly on {0,1}, n = 0 exactly on {6,7} ----
FULL, EMPTY = (0, 1), (6, 7)
dets = [d for d in range(1 << M)
        if bin(d).count("1") == N
        and all(d >> b & 1 for b in FULL) and not any(d >> b & 1 for b in EMPTY)]
Hs = np.zeros((len(dets), len(dets)), dtype=complex)
idx = {d: i for i, d in enumerate(dets)}
hnz = [(p, q) for p in range(M) for q in range(M) if abs(h[p, q]) > 1e-14]
vnz = [tuple(x) for x in np.argwhere(np.abs(V) > 1e-14)]
for j, d in enumerate(dets):
    basis = np.zeros(1 << M, dtype=complex); basis[d] = 1.0
    col = np.zeros(1 << M, dtype=complex)
    for p, q in hnz:
        col += h[p, q] * fullci.ad(p, fullci.aa(q, basis))
    for p, t, q, u in vnz:
        col += 0.25 * V[p, t, q, u] * fullci.ad(p, fullci.ad(t, fullci.aa(u, fullci.aa(q, basis))))
    for d2, i in idx.items():
        if col[d2] != 0:
            Hs[i, j] = col[d2]
Hs = 0.5 * (Hs + Hs.conj().T)
es, vec = np.linalg.eigh(Hs)
psi = np.zeros(1 << M, dtype=complex); psi[dets] = vec[:, 0]
g1 = fullci.rdm(psi, M, 1)
g2 = fullci.rdm(psi, M, 2)
l2 = g2 - closures.wedge(g1, g1)
occ = np.linalg.eigvalsh(0.5 * (g1 + g1.conj().T))
print("test state: %d determinants;  occupations %s"
      % (len(dets), np.array2string(np.sort(occ)[::-1], precision=6, suppress_small=True)))
print("   |lambda_2| = %.5f   tr gamma_1 = %.10f" % (np.abs(l2).max(), np.trace(g1).real))

C = collision(g1, l2)
print("   |C| = %.6f   hermiticity |C - C^dag| = %.2e" % (np.abs(C).max(),
                                                          np.abs(C - C.conj().T).max()))
# the commutator part must not move the extreme occupations; confirm before attributing to C
f = h + np.einsum("prqs,rs->pq", V, g1)
comm = 1j * (f.conj().T @ g1 - g1 @ f.conj().T)
print("   |commutator part| = %.6f" % np.abs(comm).max())
print()

n, U = np.linalg.eigh(0.5 * (g1 + g1.conj().T))
Cn = U.conj().T @ C @ U
i0 = np.where(np.abs(n) < 1e-10)[0]
i1 = np.where(np.abs(n - 1.0) < 1e-10)[0]
print("exact necessary conditions from the gain-loss form:")
for tag, ix, want in (("P0 C P0 (empty block) must be >= 0", i0, +1),
                      ("P1 C P1 (full block)  must be <= 0", i1, -1)):
    B = Cn[np.ix_(ix, ix)]
    ev = np.linalg.eigvalsh(0.5 * (B + B.conj().T))
    ok = (ev.min() > -1e-10) if want > 0 else (ev.max() < 1e-10)
    print("   %-38s dim %d  eigenvalues %s   %s"
          % (tag, len(ix), np.array2string(ev, precision=3e-1 and 8, suppress_small=False),
             "PASS" if ok else "FAIL"))
print()
print("diagonal rates at the extreme orbitals (these are dn_p/dt, the commutator contributes 0):")
for p in list(i1) + list(i0):
    print("   orbital %d  n = %.10f   dn/dt from C = %+.6e" % (p, n[p], np.real(Cn[p, p])))
