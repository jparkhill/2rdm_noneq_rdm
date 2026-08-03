#!/usr/bin/env python3
"""Which closure is exact for N electrons: gamma_3 = 0 or lambda_3 = 0?

For N = 2, gamma_3 vanishes identically -- three annihilators on a two-electron
state give zero. But lambda_3 = gamma_3 - (antisymmetrized products of gamma_1
and lambda_2), so lambda_3 need not vanish. If it does not, then the closure
that is EXACT at N = 2 is gamma_3 = 0, which in cumulant language is a specific
NONZERO reconstruction functional for lambda_3 -- one of the free knobs in
docs/02.

Measured, not argued.
"""
import itertools
import sys
import numpy as np
sys.path.insert(0, "tools")
from gno import Fock

M = 6
rng = np.random.default_rng(7)
dim = 1 << M
popc = np.array([bin(s).count("1") for s in range(dim)])

print("M = %d spin orbitals; random correlated N-electron reference\n" % M)
print("%3s  %14s  %14s  %14s" % ("N", "max|gamma_3|", "max|lambda_3|", "max|lambda_2|"))

for N in (2, 3, 4):
    psi = np.where(popc == N, rng.normal(size=dim), 0.0)
    if not psi.any():
        continue
    F = Fock(M, psi)

    g3 = l3 = l2 = 0.0
    idx = range(M)
    for up in itertools.combinations(idx, 3):
        for lo in itertools.combinations(idx, 3):
            g3 = max(g3, abs(F.rdm(list(up), list(lo))))
            l3 = max(l3, abs(F.lam(list(up), list(lo))))
    for up in itertools.combinations(idx, 2):
        for lo in itertools.combinations(idx, 2):
            l2 = max(l2, abs(F.lam(list(up), list(lo))))

    print("%3d  %14.3e  %14.3e  %14.3e" % (N, g3, l3, l2))

print()
print("Reading: gamma_3 = 0 exactly at N = 2 as expected. If lambda_3 is NOT")
print("zero there, the exact 2-electron closure is gamma_3 = 0, i.e. a nonzero")
print("lambda_3 reconstruction, and 'lambda_3 = 0' is already an approximation")
print("even in the case that ought to be exact.")
