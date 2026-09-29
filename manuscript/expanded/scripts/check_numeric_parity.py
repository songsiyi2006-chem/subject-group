"""Compare numeric tokens in corresponding bilingual tables without unit pooling."""
import json,re,hashlib
from pathlib import Path
BASE=Path(__file__).resolve().parents[1]
def tables(text):
    return [m.group().splitlines()[2:] for m in re.finditer(r'(?m)^\|.*(?:\n\|.*)+',text)]
def tokens(row):
    # Normalize common typographic signs and thousands grouping; preserve order.
    row=re.sub(r'(?<=\d),(?=\d{3}(?:\D|$))','',row)
    row=row.replace('−','-').replace('–','-').replace('—','-')
    return re.findall(r'(?<![A-Za-z\d])[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?',row)
texts={x:(BASE/f'manuscript_{x}.md').read_text('utf-8') for x in ['english','chinese']}
en=tables(texts['english']);cn=tables(texts['chinese']);differences=[]
for i,(et,ct) in enumerate(zip(en,cn),1):
    if len(et)!=len(ct):differences.append({'table':i,'type':'row_count','en':len(et),'cn':len(ct)})
    for j,(er,cr) in enumerate(zip(et,ct),1):
        a,b=tokens(er),tokens(cr)
        if a!=b:differences.append({'table':i,'row':j,'english_tokens':a,'chinese_tokens':b})
result={'scope':'Literal numeric-token comparison of corresponding table rows; does not prove semantic translation equivalence.','table_count':{'english':len(en),'chinese':len(cn)},'source_sha256':{k:hashlib.sha256(v.encode()).hexdigest() for k,v in texts.items()},'passed':len(en)==len(cn) and not differences,'differences':differences}
(BASE/'results/numeric_parity.json').write_bytes((json.dumps(result,ensure_ascii=False,indent=2)+'\n').encode())
print(json.dumps(result,ensure_ascii=False,indent=2))
