"""Three independent producers of an EOM, and one arbiter.

    wick     ordinary vacuum Wick, bare moments
    gwt      MKNO generalized normal ordering, bare or cumulant form
    fullci   exact Fock-space linear algebra -- the ground truth

Every derived EOM is validated against fullci. The two derivations are never used to
correct one another, and the two-body prefactor lives only in hamiltonian.py.
"""
