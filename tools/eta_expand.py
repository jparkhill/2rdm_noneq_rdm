#!/usr/bin/env python3
"""Apply the notebook's own DeltaToEta and read off the gain/loss structure.

MKNO cell 62 carries

    DeltaToEta :  delta^p_q  ->  eta^p_q + gamma^p_q

and CONTINUITY section 4 names eta the blocking factor. docs/01 section 4.2 says
the (1-n) factors are literally these hole contractions, falling out of normal
ordering rather than inserted.

The substitution is an IDENTITY, since eta = delta - gamma, so it cannot change any
flux. What it changes is visibility: after it, every term is unambiguously a gain
term (carrying gamma on the relevant index) or a loss term (carrying eta). The
failure mode in docs/01 was not missing blocking factors but a truncation that KEPT
gain terms while DROPPING their hole partners -- that is what gave Boltzmann and
2.81 electrons.

So: expand gamma_1-dot in the eta/gamma basis, and check whether the truncations we
have been using preserve the pairing or break it.
"""
import sys
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from qdiss.symbolic import mkno
from tools.eom1 import one as ONE1, two as TWO1


def delta_to_eta(terms):
    """Each delta^p_q becomes eta^p_q + gamma^p_q. Exact rewriting."""
    out = []
    for sign, facs in terms:
        pieces = [(sign, ())]
        for f in facs:
            name, rank, up, lo = f
            if name == "D" and rank == 1:
                nxt = []
                for s, acc in pieces:
                    nxt.append((s, acc + (("H", 1, up, lo),)))   # eta
                    nxt.append((s, acc + (("G", 1, up, lo),)))   # gamma
                pieces = nxt
            else:
                pieces = [(s, acc + (f,)) for s, acc in pieces]
        out.extend(pieces)
    return mkno.collect(out)


def show(label, terms):
    print("%s : %d terms" % (label, len(terms)))
    for s, f in terms:
        parts = []
        for n, r, u, l in f:
            sym = {"G": "g", "H": "eta", "D": "d", "L": "lam"}.get(n, n)
            parts.append("%s%d[%s;%s]" % (sym, r, "".join(u), "".join(l)))
        print("   %+d  %s" % (s, " ".join(parts)))
    print()


print("=== one-body piece of gamma_1-dot ===")
show("before", ONE1)
show("after DeltaToEta", delta_to_eta(ONE1))

print("=== two-body piece of gamma_1-dot (this is where lambda_2 enters) ===")
show("before", TWO1)
e2 = delta_to_eta(TWO1)
show("after DeltaToEta", e2)

gain = [(s, f) for s, f in e2 if any(n == "G" and r == 1 for n, r, u, l in f)]
loss = [(s, f) for s, f in e2 if any(n == "H" and r == 1 for n, r, u, l in f)]
print("gain terms (carry gamma^1): %d      loss terms (carry eta^1): %d"
      % (len(gain), len(loss)))
print()
print("Pairing: every original term had exactly one delta, so each splits into one")
print("gain and one loss partner. A truncation is representability-safe only if it")
print("keeps or drops the two together. Both gamma_3 = 0 and the generated lambda_2")
print("act on the gamma_2 / lambda_2 FACTOR, which is common to both partners, so")
print("neither breaks the pairing -- yet both still escaped. That means the escape")
print("is not a broken pairing, and the blocking factors alone are not sufficient.")
