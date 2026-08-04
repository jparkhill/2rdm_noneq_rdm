"""Validate each derivation against FullCI independently.

Ground truth is qdiss.eom.fullci, which builds H on the Fock space and differentiates by
exact linear algebra with no symbolic input. Each derivation is compared to it separately.
The Wick-versus-GWT column is reported last and is a CONSEQUENCE, not evidence: if both
agree with FullCI they necessarily agree with each other, and if they agreed with each
other but not with FullCI they would both be wrong.

Cumulant tensors are supplied for the GWT cumulant form; lambda_3 here is the exact
remainder, so this validates the representation, not any truncation.
"""
import numpy as np
from . import fullci, gwt, wick
from .hamiltonian import check_symmetries, random_h, random_V
from qdiss.symbolic import mkno


def tensors_for(psi, M, max_rank):
    g = {k: fullci.rdm(psi, M, k) for k in range(1, max_rank + 2)}
    t = {("D", 1): np.eye(M, dtype=complex)}
    for k, v in g.items():
        t[("G", k)] = v
    g1 = g[1]
    t[("L", 2)] = g[2] - (np.einsum("pq,vw->pvqw", g1, g1)
                          - np.einsum("pw,vq->pvqw", g1, g1))
    if 3 in g:
        up, lo = ("p", "v", "x"), ("q", "w", "y")
        sep = np.zeros((M,) * 6, dtype=complex)
        for s, facs in mkno.gamma_to_cumulants(up, lo):
            if any(n == "L" and r == 3 for n, r, u, l in facs):
                continue
            sub = ",".join("".join(u) + "".join(l) for n, r, u, l in facs) + "->pvxqwy"
            sep = sep + s * np.einsum(sub, *[t[(n, r)] for n, r, u, l in facs], optimize=True)
        t[("L", 3)] = g[3] - sep
    return t


def run(cases=((4, 2), (6, 3), (8, 4)), ranks=(1, 2), seed=17):
    print("ground truth: qdiss.eom.fullci (Fock-space H, exact linear algebra)")
    print("  derivation A: %s" % wick.METHOD)
    print("  derivation B: %s" % gwt.METHOD)
    print()
    hdr = ("%-9s %5s | %12s %12s %12s | %12s"
           % ("(M,N)", "rank", "A vs FullCI", "B vs FullCI", "B_cum vs FCI", "A vs B"))
    print(hdr)
    print("-" * len(hdr))
    worst = 0.0
    for M, N in cases:
        rng = np.random.default_rng(seed)
        h, V = random_h(M, rng), random_V(M, rng)
        sym = check_symmetries(V)
        if max(sym.values()) > 1e-12:
            raise AssertionError("V violates the convention: %s" % sym)
        psi = fullci.random_state(M, N, rng)
        H = fullci.build_H(h, V, M)
        for rank in ranks:
            t = tensors_for(psi, M, rank)
            truth = fullci.rdm_dot(psi, H, M, rank)
            a = wick.gamma_dot(rank, h, V, t, M)
            b = gwt.gamma_dot(rank, h, V, t, M)
            bc = gwt.gamma_dot(rank, h, V, t, M, cumulant=True)
            ea = np.abs(a - truth).max()
            eb = np.abs(b - truth).max()
            ec = np.abs(bc - truth).max()
            worst = max(worst, ea, eb, ec)
            print("%-9s %5d | %12.3e %12.3e %12.3e | %12.3e"
                  % ("(%d,%d)" % (M, N), rank, ea, eb, ec, np.abs(a - b).max()))
    print()
    print("worst error of any derivation against FullCI: %.3e" % worst)
    print("VERDICT: %s" % ("ALL DERIVATIONS VALIDATED" if worst < 1e-10 else "FAILURE"))
    return worst


if __name__ == "__main__":
    import sys
    sys.exit(0 if run() < 1e-10 else 1)
