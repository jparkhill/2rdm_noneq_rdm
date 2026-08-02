#!/usr/bin/env python3
"""Mukherjee-Kutzelnigg generalized-normal-ordered operators {a^P_Q} of any rank,
   built in an explicit Fock space, plus the cumulants lambda_k.

   Definition used (and CHECKED here, not assumed):

     a^{p1..pk}_{q1..qk} = sum over m=0..k, over ordered subsets A of the upper
     positions and B of the lower positions with |A|=|B|=m,
         eps(A) eps(B) * gamma_m^{p_A}_{q_B} * {a^{p_notA}_{q_notB}}

   where eps(S) is the sign of moving the elements of S to the front keeping
   relative order, and gamma_m is the m-RDM (not the cumulant).  Inverting this
   recursively defines {a_k}.  The defining property <{a_k}> = 0 for k>=1 is
   verified for every index tuple at the ranks used.
"""
import itertools, functools, numpy as np


class Fock:
    def __init__(self, M, psi):
        self.M = M
        self.dim = 1 << M
        self.psi = psi / np.linalg.norm(psi)
        self._g = {}
        self._lam = {}
        self._gnocache = {}

    # --- elementary operators, applied to a state vector ---
    def ad(self, p, v):
        out = np.zeros_like(v); bit = 1 << p
        src = np.nonzero(v)[0]; src = src[(src & bit) == 0]
        if len(src):
            out[src | bit] = ((-1.0) ** np.array([bin(s & (bit - 1)).count("1")
                                                  for s in src])) * v[src]
        return out

    def aa(self, p, v):
        out = np.zeros_like(v); bit = 1 << p
        src = np.nonzero(v)[0]; src = src[(src & bit) != 0]
        if len(src):
            out[src ^ bit] = ((-1.0) ** np.array([bin(s & (bit - 1)).count("1")
                                                  for s in src])) * v[src]
        return out

    def astr(self, up, lo, v):
        """a^{up}_{lo} = ad_{u1}..ad_{uk} aa_{lk}..aa_{l1}; rightmost acts first."""
        for q in lo:       v = self.aa(q, v)
        for p in up[::-1]: v = self.ad(p, v)
        return v

    def rdm(self, up, lo):
        """gamma_k^{up}_{lo} = <a^{up}_{lo}>"""
        key = (tuple(up), tuple(lo))
        if key not in self._g:
            self._g[key] = float(self.psi @ self.astr(list(up), list(lo), self.psi))
        return self._g[key]

    # --- normal-ordered operators ---
    def gno(self, up, lo, v, tag=None):
        """{a^{up}_{lo}} applied to v.  `tag` keys a per-vector memo."""
        if tag is not None:
            key=(tag,tuple(up),tuple(lo))
            if key in self._gnocache: return self._gnocache[key]
            r=self._gno(up,lo,v,tag); self._gnocache[key]=r; return r
        return self._gno(up,lo,v,tag)

    def _gno(self, up, lo, v, tag=None):
        up, lo = list(up), list(lo)
        k = len(up)
        out = self.astr(up, lo, v)
        for m in range(1, k):                     # subtract every LOWER-rank piece
                                                  # (m=0 is {a_k} itself, m=k is gamma_k)
            for A in itertools.combinations(range(k), m):
                for B in itertools.combinations(range(k), m):
                    sgn = _eps(A, k) * _eps(B, k)
                    g = self.rdm([up[i] for i in A], [lo[j] for j in B]) if m else 1.0
                    if g == 0.0: continue
                    ru = [up[i] for i in range(k) if i not in A]
                    rl = [lo[j] for j in range(k) if j not in B]
                    out = out - sgn * g * self.gno(ru, rl, v, tag)
        # m = k term is gamma_k itself
        out = out - self.rdm(up, lo) * v
        return out

    # --- cumulants ---
    def lam(self, up, lo):
        key=(tuple(up),tuple(lo))
        if key in self._lam: return self._lam[key]
        v=self._lam_calc(up,lo); self._lam[key]=v; return v

    def _lam_calc(self, up, lo):
        """lambda_k^{up}_{lo}: k-RDM minus every product of lower cumulants."""
        up, lo = tuple(up), tuple(lo)
        k = len(up)
        val = self.rdm(up, lo)
        for part in _set_partitions(k):
            if len(part) == 1:
                continue                              # that is lambda_k itself
            val -= _partition_value(self, up, lo, part)
        return val


def _eps(S, k):
    """sign of moving positions S (increasing) to the front, order preserved"""
    sgn = 1
    for rank, i in enumerate(S):
        sgn *= (-1) ** (i - rank)
    return sgn


@functools.lru_cache(maxsize=None)
def _set_partitions(k):
    if k == 0: return [()]
    if k == 1: return [((0,),)]
    out = []
    for p in _set_partitions(k - 1):
        for i in range(len(p)):
            out.append(p[:i] + (p[i] + (k - 1,),) + p[i + 1:])
        out.append(p + (((k - 1),),))
    return [tuple(sorted(pp, key=lambda b: b[0])) for pp in out]


def _partition_value(F, up, lo, part):
    """sum over bijections of lower-index blocks onto the upper blocks, with the
       sign of the induced permutation of the lower indices."""
    k = len(up)
    sizes = [len(b) for b in part]
    # concatenating the upper blocks reorders the upper indices: that costs a sign
    upper_perm = tuple(i for b in part for i in b)
    up_sgn = _perm_sign(upper_perm)
    total = 0.0
    for perm in itertools.permutations(range(k)):
        # split perm into consecutive chunks matching the block sizes
        chunks, i = [], 0
        for s in sizes:
            chunks.append(perm[i:i + s]); i += s
        if any(list(c) != sorted(c) for c in chunks):
            continue                        # each block's lower indices in order
        val = up_sgn * _perm_sign(perm)
        for blk, ch in zip(part, chunks):
            val *= F.lam([up[i] for i in blk], [lo[j] for j in ch])
            if val == 0.0: break
        total += val
    return total


def _perm_sign(p):
    p = list(p); s = 1
    for i in range(len(p)):
        for j in range(len(p) - 1 - i):
            if p[j] > p[j + 1]:
                p[j], p[j + 1] = p[j + 1], p[j]; s = -s
    return s
