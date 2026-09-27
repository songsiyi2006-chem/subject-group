"""Execute the exact extracted user script; do not execute the DFT inputs."""
from pathlib import Path
import hashlib,json,os,platform,subprocess,sys,time
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1]
def main():
    script=ROOT/'source/run_closed_loop_platform.py'; folder=ROOT/'results/original'
    folder.mkdir(parents=True,exist_ok=True); start=time.perf_counter()
    env=os.environ.copy()
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:env[key]='1'
    p=subprocess.run([sys.executable,str(script)],cwd=folder,env=env,capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=180)
    for name,txt in [('stdout',p.stdout),('stderr',p.stderr)]:
        for local in sorted({str(ROOT.resolve()),sys.executable,sys.prefix,sys.base_prefix},key=len,reverse=True):txt=txt.replace(local,'<local-runtime-or-project>')
        (folder/f'{name}.log').write_text(txt,encoding='utf-8',newline='\n')
    record=dict(returncode=p.returncode,script_sha256=hashlib.sha256(script.read_bytes()).hexdigest(),elapsed_s=time.perf_counter()-start,timestamp_utc=datetime.now(timezone.utc).isoformat(),threads=1,timeout_s=180,python=platform.python_version(),result_written=(folder/'closed_loop_results.json').is_file(),dft_executed=False,wet_lab_executed=False)
    (folder/'execution.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(record));return p.returncode
if __name__=='__main__':sys.exit(main())
