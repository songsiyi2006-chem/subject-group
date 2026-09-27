"""Execute the unchanged supplied script and preserve its result and diagnostics."""
import os
for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[name]='1'
import hashlib,json,platform,subprocess,sys,time
from datetime import datetime,timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def main():
    script=ROOT/'source/run_research_engine_original.py';folder=ROOT/'results/original'
    folder.mkdir(parents=True,exist_ok=True);start=time.perf_counter()
    result=subprocess.run([sys.executable,str(script)],cwd=folder,capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=300)
    for label,txt in [('stdout',result.stdout),('stderr',result.stderr)]:
        for local in sorted({str(ROOT.resolve()),sys.executable,sys.prefix,sys.base_prefix},key=len,reverse=True):txt=txt.replace(local,'<local-runtime-or-project>')
        (folder/f'{label}.log').write_text(txt,encoding='utf-8')
    record=dict(returncode=result.returncode,script_sha256=hashlib.sha256(script.read_bytes()).hexdigest(),elapsed_s=time.perf_counter()-start,timestamp_utc=datetime.now(timezone.utc).isoformat(),threads=1,timeout_s=300,python=platform.python_version(),result_written=(folder/'research_grade_results.json').exists())
    (folder/'execution.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(record));return result.returncode
if __name__=='__main__':sys.exit(main())
