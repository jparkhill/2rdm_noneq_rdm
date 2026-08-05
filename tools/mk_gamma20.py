#!/usr/bin/env python3
"""MK_gamma20: propagate the rank-1 closure gamma_2 = 0 and check 1-RDM positivity.

Derivation: qdiss.eom.wick / qdiss.eom.gwt term lists (identical), validated against
qdiss.eom.fullci at 1e-15. State: gamma_1 alone, with lambda_2 := -(gamma_1 wedge gamma_1)
so the implied gamma_2 is exactly zero -- the rank-1 analogue of MK_gamma30.

Compared against TDHF (lambda_2 = 0) and exact propagation. Reported: the full natural
occupation spectrum, its drift from t=0, whether it leaves [0,1], and accuracy against FullCI.
"""
import sys
import numpy as np
sys.path.insert(0, ".")
from qdiss.eom import fullci, wick
from qdiss.eom.hamiltonian import random_V, random_h

M, N = 6, 3
EYE = np.eye(M, dtype=complex)


def wedge(a, b):
    return np.einsum("pq,vw->pvqw", a, b) - np.einsum("pw,vq->pvqw", a, b)


def rk4(f, y, dt):
    k1 = f(y); k2 = f(y + dt / 2 * k1); k3 = f(y + dt / 2 * k2); k4 = f(y + dt * k3)
    return y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


rng = np.random.default_rng(7)
h, V = random_h(M, rng), random_V(M, rng)
psi0 = fullci.random_state(M, N, rng)
H = fullci.build_H(h, V, M)
w, U = np.linalg.eigh(H)
g1_0 = fullci.rdm(psi0, M, 1)
occ0 = np.linalg.eigvalsh(0.5 * (g1_0 + g1_0.conj().T))


def d1(g1, mode):
    g2 = np.zeros((M,) * 4, dtype=complex) if mode == "g20" else wedge(g1, g1)
    return wick.gamma_dot(1, h, V, {("D", 1): EYE, ("G", 1): g1, ("G", 2): g2}, M)


print("initial natural occupations: %s" % np.array2string(occ0[::-1], precision=5))
print()
dt, nsteps, every = 0.002, 500, 100
st = {m: g1_0.copy() for m in ("g20", "hf")}
print("%6s | %-30s | %-30s | %9s"
      % ("t", "MK_gamma20 (gamma_2=0)", "TDHF (lambda_2=0)", "exact"))
print("%6s | %9s %9s %9s | %9s %9s %9s | %9s"
      % ("", "n_min", "spec drift", "err vs FCI", "n_min", "spec drift", "err vs FCI", "n_min"))
for step in range(nsteps + 1):
    if step % every == 0:
        t = step * dt
        psit = U @ (np.exp(-1j * w * t) * (U.conj().T @ psi0))
        e1 = fullci.rdm(psit, M, 1)
        row = []
        for m in ("g20", "hf"):
            g = st[m]
            o = np.linalg.eigvalsh(0.5 * (g + g.conj().T))
            row.append((o.min(), np.abs(np.sort(o) - np.sort(occ0)).max(),
                        np.abs(g - e1).max()))
        print("%6.2f | %9.5f %9.2e %9.4f | %9.5f %9.2e %9.4f | %9.5f"
              % (t, row[0][0], row[0][1], row[0][2], row[1][0], row[1][1], row[1][2],
                 np.linalg.eigvalsh(0.5 * (e1 + e1.conj().T)).min()))
    if step < nsteps:
        for m in ("g20", "hf"):
            st[m] = rk4(lambda y, mm=m: d1(y, mm), st[m], dt)
print()
print("spec drift is the largest change in the sorted occupation spectrum since t=0.")
print("Both closures are commutator flows, hence isospectral, hence 0 <= n <= 1 is")
print("preserved for free. Neither can violate Pauli; neither can change an occupation.")
