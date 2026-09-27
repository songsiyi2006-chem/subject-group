"""Explicitly update hashes AFTER reviewing intentional changes and staging files.

The manifest and mutable validation/test execution records are excluded.
"""
from pathlib import Path
import hashlib,subprocess
ROOT=Path(__file__).resolve().parents[1]
EXCLUDE={'provenance/SHA256SUMS.txt','provenance/release_validation.json','provenance/test_results.txt'}
def main():
    names=subprocess.check_output(['git','ls-files','-z'],cwd=ROOT).decode('utf-8').split('\0')
    records=[]
    for name in sorted(n for n in names if n and n not in EXCLUDE):
        p=ROOT/name
        if not p.is_file():raise SystemExit('Missing tracked file: '+name)
        records.append(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+name)
    (ROOT/'provenance/SHA256SUMS.txt').write_text('\n'.join(records)+'\n',encoding='utf-8',newline='\n')
    print(f'Hashed {len(records)} reviewed files.')
if __name__=='__main__':main()
