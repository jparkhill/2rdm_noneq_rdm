"""Index-based (non-recursive) parser for Mathematica .nb boxed format.

The 4.9MB file loads fine; the earlier crash was a recursive tree build.
Everything here is iterative and index-based.
"""
import re, sys
from collections import Counter

def load(path):
    return open(path, 'r', encoding='utf-8', errors='replace').read()

def match_bracket(s, i):
    """s[i]=='[' -> index of matching ']', string-aware. Iterative."""
    depth = 0
    n = len(s)
    while i < n:
        ch = s[i]
        if ch == '"':
            i += 1
            while i < n:
                if s[i] == '\\':
                    i += 2
                    continue
                if s[i] == '"':
                    break
                i += 1
        elif ch == '[':
            depth += 1
        elif ch == ']':
            depth -= 1
            if depth == 0:
                return i
        i += 1
    return -1

def cells(s):
    """Yield (start, end, style) for every Cell[...] whose 2nd arg is a style string."""
    for m in re.finditer(r'\bCell\[', s):
        a = m.end() - 1
        b = match_bracket(s, a)
        if b < 0:
            continue
        body = s[a + 1:b]
        # style = last top-level quoted token, cheap heuristic on the tail
        sm = re.search(r'"(Input|Output|Text|Print|Section|Subsection|Subsubsection|Title|Subtitle|Message|MessageName|Item)"', body[-300:])
        yield a, b, (sm.group(1) if sm else None)

# ---------- rendering ----------
def strip_str(t):
    t = t.strip()
    if t.startswith('"') and t.endswith('"'):
        t = t[1:-1]
    return t.replace('\\"', '"')

GREEK = {'\\[Gamma]': 'g', '\\[Eta]': 'eta', '\\[Delta]': 'd', '\\[Lambda]': 'lam',
         '\\[CapitalDelta]': 'Delta', '\\[Epsilon]': 'eps', '\\[Tau]': 'tau',
         '\\[ImaginaryI]': 'I', '\\[Omega]': 'w', '\\[Rule]': '->', '\\[Equal]': '==',
         '\\[Element]': 'in', '\\[Infinity]': 'oo', '\\[Integral]': 'Int',
         '\\[DifferentialD]': 'd', '\\[ExponentialE]': 'E', '\\[Pi]': 'Pi',
         '\\[IndentingNewLine]': '\n', '\\[LeftDoubleBracket]': '[[',
         '\\[RightDoubleBracket]': ']]', '\\[NoBreak]': '', '\\[VeryThinSpace]': '',
         '\\[InvisibleSpace]': '', '\\[CenterEllipsis]': '...',
         '\\[LeftSkeleton]': '<<', '\\[RightSkeleton]': '>>'}

def split_args(s):
    """Split a comma-separated Mathematica arg list at depth 0."""
    out, depth, cur, i, n = [], 0, [], 0, len(s)
    while i < n:
        ch = s[i]
        if ch == '"':
            j = i + 1
            while j < n:
                if s[j] == '\\':
                    j += 2; continue
                if s[j] == '"':
                    break
                j += 1
            cur.append(s[i:j + 1]); i = j + 1; continue
        if ch in '[{(':
            depth += 1
        elif ch in ']})':
            depth -= 1
        if ch == ',' and depth == 0:
            out.append(''.join(cur)); cur = []
        else:
            cur.append(ch)
        i += 1
    if cur:
        out.append(''.join(cur))
    return [x.strip() for x in out]

def render(s):
    """Iteratively collapse box structures into readable math text."""
    prev = None
    t = s
    for _ in range(60):
        if t == prev:
            break
        prev = t
        # Subsuperscript template -> base^{sup}_{sub}
        def sub_tmpl(m):
            inner = m.group(1)
            args = split_args(inner)
            if len(args) != 3:
                return m.group(0)
            base = strip_str(args[0])
            base = GREEK.get(base, base)
            def idx(x):
                x = x.strip()
                if x.startswith('{') and x.endswith('}'):
                    x = x[1:-1]
                return ''.join(strip_str(y) for y in split_args(x))
            return f"{base}^{{{idx(args[2])}}}_{{{idx(args[1])}}}"
        t = re.sub(r'TemplateBox\[\{((?:[^\[\]]|\[[^\[\]]*\])*)\},\s*\n?\s*"Subsuperscript"\]',
                   sub_tmpl, t)
        # RowBox[{a,b,c}] -> concatenation
        def rb(m):
            inner = m.group(1)
            return ''.join(strip_str(x) for x in split_args(inner))
        t = re.sub(r'RowBox\[\{((?:[^\[\]{}]|\{[^{}]*\})*)\}\]', rb, t)
        t = re.sub(r'SuperscriptBox\[([^\[\]]*),([^\[\]]*)\]',
                   lambda m: f"{strip_str(m.group(1))}^{{{strip_str(m.group(2))}}}", t)
        t = re.sub(r'SubscriptBox\[([^\[\]]*),([^\[\]]*)\]',
                   lambda m: f"{strip_str(m.group(1))}_{{{strip_str(m.group(2))}}}", t)
        t = re.sub(r'FractionBox\[([^\[\]]*),([^\[\]]*)\]',
                   lambda m: f"({strip_str(m.group(1))})/({strip_str(m.group(2))})", t)
        t = re.sub(r'StyleBox\[([^\[\]]*?),\s*"[A-Za-z]+"\]', lambda m: strip_str(m.group(1)), t)
        t = re.sub(r'TagBox\[([^\[\]]*?),[^\[\]]*\]', lambda m: strip_str(m.group(1)), t)
    for k, v in GREEK.items():
        t = t.replace(k, v)
    t = re.sub(r'\s+', ' ', t)
    return t.strip()

if __name__ == '__main__':
    s = load(sys.argv[1])
    c = Counter()
    tot = 0
    for a, b, st in cells(s):
        tot += 1
        c[st] += 1
    print("Cell[ occurrences:", tot)
    for k, v in c.most_common():
        print(f"  {k}: {v}")
