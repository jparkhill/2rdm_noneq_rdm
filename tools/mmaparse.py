"""Minimal, iterative Mathematica box-expression parser.

parse_expr(text) -> nested Python structure:
   Head(name, args) | list | str | int | float | symbol
RowBox[{...}] is collapsed by concatenating its children's *source text*,
which is exactly how the FrontEnd renders it, then re-parsing.
"""
import re, sys

TOK = re.compile(r'''
    (?P<str>"(?:[^"\\]|\\.)*")
  | (?P<num>-?\d+\.?\d*(?:`[\d.]*)?(?:\*\^-?\d+)?)
  | (?P<sym>\\\[[A-Za-z]+\]|[A-Za-z$][A-Za-z0-9$]*)
  | (?P<op>->|\[|\]|\{|\}|\(|\)|,|;|\*|\+|/|\^|&)
  | (?P<ws>\s+)
''', re.X)

def tokenize(s):
    i, n = 0, len(s)
    out = []
    while i < n:
        m = TOK.match(s, i)
        if not m:
            i += 1
            continue
        i = m.end()
        k = m.lastgroup
        if k == 'ws':
            continue
        out.append((k, m.group()))
    return out

def unq(t):
    t = t[1:-1]
    t = t.replace('\\<','').replace('\\>','')
    t = re.sub(r'\\\n\s*','',t)
    t = t.replace('\\"','"').replace('\\\\','\\')
    return t

class P:
    def __init__(self, toks):
        self.t = toks; self.i = 0
    def peek(self):
        return self.t[self.i] if self.i < len(self.t) else (None, None)
    def next(self):
        v = self.t[self.i]; self.i += 1; return v
    def parse(self):
        k, v = self.peek()
        if k is None: return None
        if k == 'str':
            self.next(); return unq(v)
        if k == 'num':
            self.next()
            try: return int(v)
            except ValueError:
                try: return float(v.split('`')[0])
                except ValueError: return v
        if k == 'op' and v == '{':
            self.next(); items = []
            while True:
                k2, v2 = self.peek()
                if k2 is None: break
                if v2 == '}': self.next(); break
                if v2 == ',': self.next(); continue
                items.append(self.parse())
            return items
        if k == 'sym':
            self.next(); name = v
            k2, v2 = self.peek()
            if v2 == '[':
                self.next(); args = []
                while True:
                    k3, v3 = self.peek()
                    if k3 is None: break
                    if v3 == ']': self.next(); break
                    if v3 == ',': self.next(); continue
                    args.append(self.parse())
                return (name, args)
            return name
        self.next()
        return v

def parse(s):
    return P(tokenize(s)).parse()

def flatten_rowbox(e):
    """Recursively drop RowBox/StyleBox/TagBox wrappers, keeping content."""
    if isinstance(e, tuple):
        h, a = e
        if h in ('RowBox','StyleBox','TagBox','InterpretationBox','BoxData','FormBox'):
            inner = [flatten_rowbox(x) for x in a]
            if len(inner) == 1: return inner[0]
            return inner
        return (h, [flatten_rowbox(x) for x in a])
    if isinstance(e, list):
        return [flatten_rowbox(x) for x in e]
    return e

def clean_struct(e):
    """After flatten, RowBox lists contain literal '{' ',' '}' strings; rebuild real lists."""
    if isinstance(e, list):
        parts = [clean_struct(x) for x in e]
        # if this list is a rendered brace group: ['{', a, ',', b, '}']
        flat = []
        for p in parts:
            if isinstance(p, list) and any(isinstance(q,str) and q in '{}' for q in p):
                flat.append(p)
            else:
                flat.append(p)
        if parts and parts[0] == '{' and parts[-1] == '}':
            items, cur = [], []
            for p in parts[1:-1]:
                if p == ',':
                    if len(cur) == 1: items.append(cur[0])
                    elif cur: items.append(cur)
                    cur = []
                else:
                    cur.append(p)
            if len(cur) == 1: items.append(cur[0])
            elif cur: items.append(cur)
            return items
        return parts
    if isinstance(e, tuple):
        return (e[0], [clean_struct(x) for x in e[1]])
    return e

def load_cell_struct(raw):
    return clean_struct(flatten_rowbox(parse(raw)))
