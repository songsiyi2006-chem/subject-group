"""Build an offline explorer from actual saved electronic energies."""
from pathlib import Path
import csv,hashlib,json

BASE=Path(__file__).resolve().parents[1]


def main():
    source=BASE/'results/extensions/correlation/curve.csv'
    template=BASE/'scripts/quantum_explorer_template.html'
    rows=list(csv.DictReader(source.open(encoding='utf-8-sig',newline='')))
    data=[]
    for row in rows:
        if row['status'] not in ('success','converged','ok','passed'):continue
        value=row.get('S2','')
        data.append({'basis':row['basis'],'method':row['method'].upper(),'R_A':float(row['R_A']),
                     'energy_Hartree':float(row['energy_Hartree']),'S2':float(value) if value and value.lower() not in ('nan','none') else None})
    if len(data)!=150:raise ValueError('Expected the 150 converged main curve records before building.')
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    text=template.read_text(encoding='utf-8').replace('__DATA_SHA__',sha(source)).replace('__DATA_JSON__',json.dumps(data,ensure_ascii=False,allow_nan=False).replace('</','<\\/'))
    output=BASE/'reports/quantum_explorer.html'
    output.write_text(text,encoding='utf-8',newline='\n')
    record={'rows':len(data),'quantum_calls_in_browser':0,'network_dependencies':[],
            'scope':'Offline display of saved total energies and spin diagnostics; no additional scientific calculation.',
            'inputs_sha256':{p.relative_to(BASE).as_posix():sha(p) for p in (source,template,Path(__file__))},
            'output_sha256':sha(output)}
    (BASE/'results/extensions/explorer_manifest.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'rows':len(data),'file':output.relative_to(BASE).as_posix()}))


if __name__=='__main__':main()
