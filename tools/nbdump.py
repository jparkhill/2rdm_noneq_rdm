"""Dump the whole notebook to readable markdown: headings + rendered Input/Output/Text."""
import re, sys
sys.path.insert(0,'tools')
from nbparse import load, cells, render, strip_str

def title_of(body):
    m = re.match(r'\s*"((?:[^"\\]|\\.)*)"', body, re.S)
    if m:
        t = m.group(1).replace('\\<','').replace('\\>','').replace('\\\n','')
        t = re.sub(r'\\+','',t)
        return re.sub(r'\s+',' ',t).strip()
    r = render(body[:3000])
    return re.sub(r'\s+',' ',r)[:120]

src, dst = sys.argv[1], sys.argv[2]
s = load(src)
lvl = {'Title':1,'Section':2,'Subsection':3,'Subsubsection':4}
out = []
for i,(a,b,st) in enumerate(cells(s)):
    body = s[a+1:b]
    if st in lvl:
        out.append(f"\n\n{'#'*lvl[st]} [{i}] {title_of(body)}\n")
    elif st == 'Text':
        m = re.match(r'\s*"((?:[^"\\]|\\.)*)"', body, re.S)
        t = strip_str(m.group(1)) if m else ''
        t = t.replace('\\\n','').replace('\\','')
        if t.strip():
            out.append(f"> {re.sub(chr(10),chr(10)+'> ',t.strip())}\n")
    elif st in ('Input','Output'):
        if 'RasterBox' in body or 'GraphicsBox' in body:
            out.append(f"`[{st} {i}: graphics omitted]`\n")
            continue
        r = render(body)
        r = re.sub(r',\s*(CellChangeTimes|CellLabel|CellID|ExpressionUUID|FontSize|Background|CellChangeTimes)->.*$','',r,flags=re.S)
        r = r.replace('BoxData','').strip()
        if len(r) > 6000:
            r = r[:6000] + f' ...[TRUNC {len(r)}]'
        if r:
            out.append(f"**{st} [{i}]**\n```\n{r}\n```\n")
open(dst,'w').write(''.join(out))
print("wrote", dst, len(''.join(out)), "chars")
