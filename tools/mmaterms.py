"""Extract derived tensor terms from MKNO.nb cells rendered in *typeset* form.

Cells such as 430 render each term as a product of
    TemplateBox[{name, {sub...}, {sup...}}, "Subsuperscript"]
factors, optionally preceded by a rational coefficient (FractionBox) and
separated by "+"/"-" tokens.  This is the reliable representation: the
list-form outputs (417/426/448) are the ones where the notebook itself
prints "Unable to Parse Expression back into List...".

Subsuperscript[x, sub, sup]  ->  we record (name, sub_indices, sup_indices).
"""
import re
from fractions import Fraction

TB = re.compile(r'TemplateBox\[\{(.*?)\},\s*\n?\s*"Subsuperscript"\]', re.S)
IDX = re.compile(r'\\\\?"([a-zA-Z])\\\\?"')
SYM = {'[Gamma]': 'gamma', '[Eta]': 'eta', '[Delta]': 'delta',
       '[Lambda]': 'lambda'}


def _name(tok):
    t = tok.strip().replace('\\\\', '').replace('\\"', '').replace('"', '')
    t = t.replace('\\<', '').replace('\\>', '').lstrip('\\').strip()
    return SYM.get(t, t)


def _tb_parse(inner):
    """inner = 'name,RowBox[{...sub...}],RowBox[{...sup...}]' -> (name, sub, sup)"""
    # split at top level commas
    depth, parts, cur = 0, [], ''
    for ch in inner:
        if ch in '[{(':
            depth += 1
        elif ch in ']})':
            depth -= 1
        if ch == ',' and depth == 0:
            parts.append(cur); cur = ''
        else:
            cur += ch
    parts.append(cur)
    if len(parts) != 3:
        return None
    name = _name(parts[0])
    def idxs(p):
        return IDX.findall(p)
    return (name, idxs(parts[1]), idxs(parts[2]))


def extract_typeset(raw):
    """Return [(Fraction coeff, [(name, sub, sup), ...]), ...]."""
    # tokenise the cell into an ordered stream of signs, fractions, tensors
    toks = []
    for m in re.finditer(
            r'TemplateBox\[\{.*?\},\s*\n?\s*"Subsuperscript"\]'
            r'|FractionBox\[\s*"(\d+)"\s*,\s*"(\d+)"\s*\]'
            r'|"([-+])"', raw, re.S):
        s = m.group(0)
        if s.startswith('TemplateBox'):
            inner = TB.match(s)
            t = _tb_parse(inner.group(1)) if inner else None
            if t:
                toks.append(('T', t))
        elif s.startswith('FractionBox'):
            toks.append(('C', Fraction(int(m.group(1)), int(m.group(2)))))
        else:
            toks.append(('S', m.group(3)))

    terms, coeff, sign, fac = [], Fraction(1), 1, []
    def flush():
        if fac:
            terms.append((sign * coeff, list(fac)))
    for kind, v in toks:
        if kind == 'S':
            flush()
            coeff, fac = Fraction(1), []
            sign = -1 if v == '-' else 1
        elif kind == 'C':
            coeff = v
        else:
            fac.append(v)
    flush()
    return terms


def fmt(t):
    PR = {'gamma': 'γ', 'eta': 'η', 'delta': 'δ', 'lambda': 'λ'}
    c, fs = t
    body = ' '.join(f"{PR.get(n,n)}^{{{''.join(sup)}}}_{{{''.join(sub)}}}"
                    for n, sub, sup in fs)
    return f"{'+' if c > 0 else '-'} {abs(c)} {body}"


if __name__ == '__main__':
    import sys, collections
    sys.path.insert(0, '/root/workspace/tools')
    from nbparse import load, cells
    s = load('/root/workspace/refs/MukherjeeKutzelniggNO/MKNO.nb')
    allc = list(cells(s))
    for idx in [int(x) for x in sys.argv[1:]]:
        a, b, st = allc[idx]
        T = extract_typeset(s[a + 1:b])
        print(f"\n### cell {idx} [{st}] -> {len(T)} terms")
        prof = collections.Counter(len(f) for _, f in T)
        nm = collections.Counter(n for _, f in T for n, _, _ in f)
        print("   factors/term:", dict(sorted(prof.items())), " names:", dict(nm))
        for t in T[:6]:
            print("   ", fmt(t))
