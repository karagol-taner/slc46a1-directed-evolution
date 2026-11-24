#!/usr/bin/env bash
set -euo pipefail

mkdir -p seq
cd seq

curl -sS -L --fail --max-time 90 "https://rest.uniprot.org/uniprotkb/Q96NT5.fasta" -o Q96NT5.fasta
curl -sS -L --fail --max-time 90 "https://rest.uniprot.org/uniprotkb/F1NJ67.fasta" -o F1NJ67.fasta
curl -sS -L --fail --max-time 90 "https://www.rcsb.org/fasta/entry/7BC6" -o 7BC6.fasta
curl -sS -L --fail --max-time 90 "https://files.rcsb.org/download/7BC6.cif" -o 7BC6.cif
curl -sS -L --fail --max-time 90 "https://files.rcsb.org/download/7BC7.cif" -o 7BC7.cif

python3 - <<'PY'
import hashlib, os, sys
expect = {
    'Q96NT5.fasta': 459,
    '7BC6.cif': None,
    '7BC7.cif': None,
}
seq = ''.join(l.strip() for l in open('Q96NT5.fasta') if not l.startswith('>'))
if len(seq) != 459:
    sys.exit(f'FATAL Q96NT5 length {len(seq)}, expected 459')
ent = {}
for rec in open('7BC6.fasta').read().split('>'):
    if rec.strip():
        h, *r = rec.strip().split('\n')
        ent[h.split('|')[0]] = ''.join(r)
if len(ent['7BC6_1']) != 480:
    sys.exit(f'FATAL 7BC6 chain A entity length {len(ent["7BC6_1"])}, expected 480')
for f in ('7BC6.cif', '7BC7.cif'):
    if os.path.getsize(f) < 100000:
        sys.exit(f'FATAL {f} is {os.path.getsize(f)} bytes, truncated download')
    print(f, 'sha256', hashlib.sha256(open(f,'rb').read()).hexdigest()[:16])
print('sources fetched and checked')
PY
