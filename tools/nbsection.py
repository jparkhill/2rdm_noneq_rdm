import sys, re
sys.path.insert(0, 'tools')
from nbparse import load, cells, render, strip_str

s = load(sys.argv[1])
want = sys.argv[2].lower()
maxcells = int(sys.argv[3]) if len(sys.argv) > 3 else 40

allc = list(cells(s))
# find heading cells
heads = []
for i, (a, b, st) in enumerate(allc):
    if st in ('Section', 'Subsection', 'Subsubsection', 'Title'):
        body = s[a+1:b]
        m = re.match(r'\s*"((?:[^"\\]|\\.)*)"', body)
        title = strip_str(m.group(1)) if m else '?'
        title = title.replace('\\\n', '').replace('\\', '')
        heads.append((i, a, st, title.strip()))

hit = [h for h in heads if want in h[3].lower()]
if not hit:
    print("no heading matched. available:")
    for _, _, st, t in heads:
        print(f"  [{st}] {t}")
    sys.exit(1)

for (i, a, st, title) in hit:
    # end = next heading of same or higher level, by cell index
    order = {'Title':0,'Section':1,'Subsection':2,'Subsubsection':3}
    nxt = len(allc)
    for (j, aj, stj, tj) in heads:
        if j > i and order.get(stj,9) <= order.get(st,9):
            nxt = j; break
    print("="*90)
    print(f"[{st}] {title}   (cells {i}..{nxt})")
    print("="*90)
    n = 0
    for k in range(i+1, nxt):
        aa, bb, sk = allc[k]
        if sk not in ('Input','Output','Text','Print'):
            continue
        body = s[aa+1:bb]
        if sk == 'Text':
            m = re.match(r'\s*"((?:[^"\\]|\\.)*)"', body, re.S)
            txt = strip_str(m.group(1)) if m else body[:200]
            print(f"\n--- TEXT ---\n{txt.replace(chr(92)+chr(10),'')[:900]}")
            continue
        r = render(body)
        r = re.sub(r',\s*(CellChangeTimes|CellLabel|CellID|FontSize|Background|ExpressionUUID)->.*$', '', r, flags=re.S)
        if len(r) > 4000:
            r = r[:4000] + f"  ...[truncated, {len(r)} chars]"
        print(f"\n--- {sk} (cell {k}) ---\n{r}")
        n += 1
        if n >= maxcells:
            print("\n...[cell limit reached]")
            break
