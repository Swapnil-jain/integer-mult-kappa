"""Portable exact finite checks; Python 3.11+, standard library only.

Generated entirely by OpenAI Codex; no independent human mathematical review.
"""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]


def main():
    if sys.flags.optimize:raise RuntimeError('Do not use -O: producer assertions are required')
    manifest=json.loads((HERE/'manifest.json').read_text())
    for name,digest in manifest['files'].items():
        path=HERE/name
        if hashlib.sha256(path.read_bytes()).hexdigest()!=digest:
            raise RuntimeError('Certificate/source integrity mismatch: '+name)
    for name,digest in manifest['upstream_files'].items():
        # Text was normalized to LF when the pinned sources were retrieved.
        content=(REPO/name).read_text(encoding='utf-8').replace('\r\n','\n').rstrip()+'\n'
        if hashlib.sha256(content.encode()).hexdigest()!=digest:
            raise RuntimeError('Upstream dependency changed; review required: '+name)
    print('Pinned inputs and source hashes PASS',flush=True)
    commands=[
        [str(HERE/'verify.py')],
        [str(HERE/'test_certificates.py')],
        [str(HERE/'reproduce_candidate.py')],
        [str(HERE/'audit_community.py'),'--verify-saved'],
        ['-m','unittest','discover','-s',str(REPO/'tests'),'-p','test_round6.py','-v'],
    ]
    for command in commands:subprocess.run([sys.executable,*command],cwd=REPO,check=True)
    print('PASS: finite arithmetic and circuit audits; inherited multiplication arguments remain assumptions.')


if __name__=='__main__':main()
