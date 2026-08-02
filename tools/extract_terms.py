"""Recover derived tensor-term lists from MKNO.nb cells as clean Python data."""
import sys, re, json
sys.path.insert(0,'/root/workspace/tools')
from nbparse import load, cells
from mmaparse import load_cell_struct

SYM = {'\\[Gamma]':'gamma', '\\[Eta]':'eta', '\\[Delta]':'delta', '\\[Lambda]':'lambda',
       '\\[CapitalDelta]':'Delta', '\\[Epsilon]':'eps'}
PRETTY = {'gamma':'γ','eta':'η','delta':'δ','lambda':'λ','ta':'ã','a':'a','VTs':'V','Vs':'V','Vt':'V'}

def toks(x, out):
    if isinstance(x, str):
        out.append(x)
    elif isinstance(x, list):
        for y in x: toks(y, out)
    elif isinstance(x, tuple):
        out.append(x[0])
        for y in x[1]: toks(y, out)
    else:
        out.append(str(x))
    return out

def canon(x):
    t = [s for s in toks(x, []) if s not in (',',)]
    return t

def extract(struct):
    """Return list of (sign, [(tensorname, [up], [lo]), ...])."""
    # normalise: strip quotes, map symbols
    def norm(s):
        s = s.strip()
        if s.startswith('"') and s.endswith('"'): s = s[1:-1]
        return SYM.get(s, s)
    flat = [norm(s) for s in canon(struct)]
    terms, cur, sign = [], [], None
    i, n = 0, len(flat)
    while i < n:
        s = flat[i]
        if s == 'Sign':
            if sign is not None:
                terms.append((sign, cur)); cur = []
            j = i+1
            neg = 1
            while j < n and flat[j] in ('-','+'):
                if flat[j] == '-': neg = -neg
                j += 1
            sign = neg * int(flat[j]) if flat[j].lstrip('-').isdigit() else neg
            i = j+1; continue
        if s in ('gamma','eta','delta','lambda','ta','a','VTs','Vs','Vt','V'):
            # next token is rank, then rank*2 index labels
            try:
                rank = int(flat[i+1])
            except (ValueError, IndexError):
                i += 1; continue
            up = flat[i+2:i+2+rank]
            lo = flat[i+2+rank:i+2+2*rank]
            if len(lo) == rank and all(re.fullmatch(r'[a-z]', z) for z in up+lo):
                cur.append((s, up, lo))
                i += 2 + 2*rank; continue
        i += 1
    if sign is not None:
        terms.append((sign, cur))
    return terms

def fmt(t):
    sg, ts = t
    return ('+ ' if sg>0 else '- ') + ' '.join(
        f"{PRETTY.get(n,n)}^{{{''.join(u)}}}_{{{''.join(l)}}}" for n,u,l in ts)

if __name__ == '__main__':
    s = load('/root/workspace/refs/MukherjeeKutzelniggNO/MKNO.nb')
    allc = list(cells(s))
    for idx in [int(x) for x in sys.argv[1:]]:
        a,b,st = allc[idx]
        T = extract(load_cell_struct(s[a+1:b]))
        from collections import Counter
        c = Counter(tuple(sorted(n for n,_,_ in ts)) for _,ts in T)
        print(f"\n### cell {idx} [{st}] -> {len(T)} terms")
        for k,v in c.most_common(12):
            print("   ", '·'.join(PRETTY.get(x,x) for x in k), "×", v)
        for t in T[:8]:
            print("   ", fmt(t))
