"""Parse MKNO.nb's internal derived-term representation into Python.

A term is  {Sign->s, {{tensor,rank}, {upper...}, {lower...}}, ... }
rendered by nbparse as e.g.
  {Sign-1}, {{g,1},{"x"},{"u"}}, {{d,1},{"v"},{"q"}} ...
This recovers a list of (sign, [(name, [up], [lo]), ...]).
"""
import re, sys
sys.path.insert(0,'tools')
from nbparse import load, cells, match_bracket, split_args, strip_str

def brace_match(s, i):
    assert s[i] == '{'
    d = 0; n = len(s)
    while i < n:
        c = s[i]
        if c == '"':
            i += 1
            while i < n and s[i] != '"':
                i += 2 if s[i] == '\\' else 1
        elif c == '{': d += 1
        elif c == '}':
            d -= 1
            if d == 0: return i
        i += 1
    return -1

def clean(x):
    x = x.strip()
    x = re.sub(r'RowBox\[|\]|"|\\<|\\>|\s', '', x)
    return x

def parse_terms(raw):
    """raw = the rendered text of a term-list output cell."""
    terms = []
    # each term begins with {Sign  (rendered as {\<Sign\>N  or "Sign"->N)
    for m in re.finditer(r'\{\s*(?:\\<)?Sign(?:\\>)?\s*(-?1)', raw):
        # find enclosing brace of this term
        st = raw.rfind('{', 0, m.start()+1)
        en = brace_match(raw, st)
        if en < 0: continue
        body = raw[st+1:en]
        sign = int(m.group(1))
        tensors = []
        # tensor factors look like {{name,rank},{up..},{lo..}}
        for tm in re.finditer(r'\{\s*\{\s*([A-Za-z]+)\s*,\s*(\d+)\s*\}\s*,', body):
            ts = tm.start()
            te = brace_match(body, ts)
            if te < 0: continue
            args = split_args(body[ts+1:te])
            if len(args) < 3: continue
            name = tm.group(1)
            up = [clean(y) for y in split_args(args[1].strip().strip('{}'))]
            lo = [clean(y) for y in split_args(args[2].strip().strip('{}'))]
            tensors.append((name, [u for u in up if u], [l for l in lo if l]))
        terms.append((sign, tensors))
    return terms

def fmt(term):
    sign, ts = term
    s = '+' if sign > 0 else '-'
    parts = []
    for name, up, lo in ts:
        nm = {'g':'γ','d':'δ','eta':'η','lam':'λ','ta':'ã','a':'a'}.get(name, name)
        parts.append(f"{nm}^{{{''.join(up)}}}_{{{''.join(lo)}}}")
    return s + ' ' + ' '.join(parts)

if __name__ == '__main__':
    src, cellidx = sys.argv[1], int(sys.argv[2])
    s = load(src)
    from nbparse import render
    allc = list(cells(s))
    a,b,st = allc[cellidx]
    raw = render(s[a+1:b])
    T = parse_terms(raw)
    print(f"cell {cellidx} [{st}]  ->  {len(T)} terms")
    from collections import Counter
    c = Counter()
    for sg,ts in T:
        c[tuple(sorted(n for n,_,_ in ts))] += 1
    print("term signatures (tensor multiset -> count):")
    for k,v in c.most_common():
        print("   ", k, v)
    for t in T[:12]:
        print("   ", fmt(t))
