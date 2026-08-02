"""Extract the MKNO-stamped term comments from the EE2 C++ sources.

Lines in RhoDotSec / RhoDotSecFull carry their provenance, e.g.

  //61{{Sign, -1}, {{Vs, 2}, {p, u}, {p, r}}, {{Vt, 2}, {p, s}, {p, t}},
       {{\[Gamma], 1}, {p}, {p}}, {{\[Gamma], 1}, {u}, {t}},
       {{\[Eta], 1}, {p}, {q}}, {{\[Eta], 1}, {s}, {r}}}

which is the derived secular term list in machine-readable form -- the exact
output of the notebook, transcribed by hand into code.
"""
import re, sys, collections

SYM = {'\\[Gamma]': 'gamma', '\\[Eta]': 'eta', '\\[Delta]': 'delta',
       'Vs': 'Vs', 'Vt': 'Vt'}
STAMP = re.compile(r'//\s*(\d+)\s*(\{\{Sign.*\}\})\s*$')
FACTOR = re.compile(r'\{\{\s*([^,{}]+?)\s*,\s*(\d)\s*\}((?:\s*,\s*\{[^{}]*\})+)\}')
GRP = re.compile(r'\{([^{}]*)\}')


def parse_stamp(txt):
    m = re.search(r'\{\s*"?Sign"?\s*,\s*(-?\s*1)\s*\}', txt)
    if not m:
        return None
    sign = -1 if '-' in m.group(1) else 1
    facs = []
    for fm in FACTOR.finditer(txt):
        name = SYM.get(fm.group(1).strip(), fm.group(1).strip())
        rank = int(fm.group(2))
        groups = [g.strip() for g in GRP.findall(fm.group(3))]
        if len(groups) != 2:
            continue
        up = [x.strip() for x in groups[0].split(',') if x.strip()]
        lo = [x.strip() for x in groups[1].split(',') if x.strip()]
        if len(up) == rank and len(lo) == rank:
            facs.append((name, tuple(up), tuple(lo)))
    return sign, tuple(facs)


def extract_file(path):
    out = []
    for ln, line in enumerate(open(path, errors='replace'), 1):
        m = STAMP.search(line.rstrip())
        if not m:
            continue
        t = parse_stamp(m.group(2))
        if t:
            out.append((ln, int(m.group(1)), t[0], t[1]))
    return out


def fmt(sign, facs):
    PR = {'gamma': 'γ', 'eta': 'η', 'delta': 'δ'}
    return ('+ ' if sign > 0 else '- ') + ' '.join(
        f"{PR.get(n,n)}^{{{''.join(u)}}}_{{{''.join(l)}}}" for n, u, l in facs)


if __name__ == '__main__':
    for path in sys.argv[1:]:
        T = extract_file(path)
        print(f"\n### {path}: {len(T)} stamped terms")
        byg = collections.Counter(g for _, g, _, _ in T)
        print("   groups:", dict(sorted(byg.items())))
        sigs = collections.Counter(tuple(n for n, _, _ in f) for _, _, _, f in T)
        for k, v in sigs.most_common(8):
            print("   ", '·'.join(k), "x", v)
        for ln, g, s, f in T[:5]:
            print(f"   L{ln} [{g}] {fmt(s,f)}")
