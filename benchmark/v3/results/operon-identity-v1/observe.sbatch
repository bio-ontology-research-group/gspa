#!/bin/bash
#SBATCH --job-name=gspa-native-gene-ids
#SBATCH --partition=debug
#SBATCH --nodelist=node004
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --no-requeue
#SBATCH --output=/data/hohndor/gspa-v3/operon-identity-v1/slurm/%j.log
set -euo pipefail
root=/data/hohndor/gspa-v3/operon-identity-v1
out=$root/observed/$SLURM_JOB_ID
work=/scratch/$SLURM_JOB_ID/native-genes
mkdir -p "$out" "$work/libs"
export TMPDIR=/scratch/$SLURM_JOB_ID/tmp
mkdir -p "$TMPDIR"
finish() { code=$?; trap - EXIT; cp "$0" "$out/executed.sbatch"; printf '%s\n' "$code" > "$out/exit-code.txt"; exit "$code"; }
trap finish EXIT
cd "$root"
sha256sum -c observe-packages.sha256 > "$out/package-check.log"
/usr/bin/python3 - "$root/wheels" "$work/libs" "$out" <<'PY'
import hashlib,json,sys,zipfile
from pathlib import Path
wheels,libs,out=map(Path,sys.argv[1:])
for p in wheels.glob('*.whl'):
 with zipfile.ZipFile(p) as z:z.extractall(libs)
native=Path('/data/hohndor/gspa-v3/antismash-observed-v1/GCA_003598085.1.json')
data=json.loads(native.read_text())
with (out/'genome.fna').open('w') as f:
 for r in data['records']:f.write('>'+(r.get('original_id') or r['id'])+'\n'+r['seq']['data']+'\n')
(out/'input.json').write_text(json.dumps(dict(native_json_sha256=hashlib.sha256(native.read_bytes()).hexdigest(),record_count=len(data['records']),fasta_sha256=hashlib.sha256((out/'genome.fna').read_bytes()).hexdigest(),python=sys.version),indent=2)+'\n')
PY
export PYTHONPATH="$work/libs"
/usr/bin/python3 -m pyrodigal --version > "$out/version.txt"
for mode in single meta; do
 /usr/bin/python3 -m pyrodigal -i "$out/genome.fna" -o "$out/$mode.gff" -a "$out/$mode.faa" -f gff -p "$mode" -j 2 > "$out/$mode.log" 2>&1
 /usr/bin/python3 "$root/before-predict_operons.py" --gff "$out/$mode.gff" --out-dir "$out/$mode-operons-before" > "$out/$mode-operons-before.log"
done
/usr/bin/python3 - "$out" <<'PY'
from pathlib import Path
import json,sys,hashlib
out=Path(sys.argv[1]);r={}
for mode in ('single','meta'):
 ids={l[1:].split()[0] for l in (out/(mode+'.faa')).read_text().splitlines() if l.startswith('>')}
 rows=[l.split('\t') for l in (out/(mode+'.gff')).read_text().splitlines() if l and not l.startswith('#')]
 native={dict(p.split('=',1) for p in row[8].split(';') if '=' in p)['ID'] for row in rows if row[2]=='CDS'}
 ops={p for l in (out/(mode+'-operons-before/operons_for_integrate.tsv')).read_text().splitlines() for p in l.split('\t')}
 r[mode]=dict(proteins=len(ids),cds=len(native),shared_native_gff_protein_ids=len(ids&native),operon_members=len(ops),operon_members_matching_proteins=len(ops&ids))
 assert r[mode]['proteins']>0 and r[mode]['operon_members']>0
(out/'identity-observation.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps(r,indent=2))
PY
