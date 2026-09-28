"""Execute the full supplied program; preserve unsupported source claims as source data."""
from pathlib import Path
import hashlib,json,os,platform,subprocess,sys,time
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1]
def main():
    script=ROOT/'source/run_deployment_and_figures.py';folder=ROOT/'results/original';folder.mkdir(parents=True,exist_ok=True)
    env=os.environ.copy();env['MPLBACKEND']='Agg';env['PYTHONIOENCODING']='utf-8'
    for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:env[k]='1'
    start=time.perf_counter();p=subprocess.run([sys.executable,str(script)],cwd=folder,env=env,capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=180)
    for name,txt in [('stdout',p.stdout),('stderr',p.stderr)]:
        for local in sorted({str(ROOT.resolve()),sys.executable,sys.prefix,sys.base_prefix},key=len,reverse=True):txt=txt.replace(local,'<local-runtime-or-project>')
        (folder/f'{name}.log').write_text(txt,encoding='utf-8',newline='\n')
    record=dict(returncode=p.returncode,script_sha256=hashlib.sha256(script.read_bytes()).hexdigest(),timestamp_utc=datetime.now(timezone.utc).isoformat(),elapsed_s=time.perf_counter()-start,python=platform.python_version(),threads=1,result_written=(folder/'deployment_toolkit_results.json').is_file(),raw_instrument_data_supplied=False,DFT_executed=False,grant_submitted=False,meaning='Execution success only; original concordance, wet-lab, DFT and equipment claims are not endorsed')
    (folder/'execution.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps(record));return p.returncode
if __name__=='__main__':sys.exit(main())
