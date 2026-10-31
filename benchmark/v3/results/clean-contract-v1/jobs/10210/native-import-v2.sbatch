#!/bin/bash
#SBATCH --job-name=gspa-clean-import
#SBATCH --partition=debug
#SBATCH --cpus-per-task=3
#SBATCH --nodelist=node002
#SBATCH --mem=16G
#SBATCH --time=00:40:00
#SBATCH --output=/data/hohndor/gspa-v3/clean-contract-v1/slurm/%j.log
set -euo pipefail
package=/data/hohndor/gspa-v3/empty-input-v1
current=/data/hohndor/gspa-v3/clean-contract-v1
work=/scratch/$SLURM_JOB_ID/gspa-clean-import
jobout=$current/jobs/$SLURM_JOB_ID
mkdir -p "$work/fs" "$jobout"
export TMPDIR=/scratch/$SLURM_JOB_ID/tmp
mkdir -p "$TMPDIR"
finish() {
 code=$?
 trap - EXIT
 if [ -d "$work/fs/home/leechuck/tmp/dbcan-clean-disabled-route" ]; then rsync -a "$work/fs/home/leechuck/tmp/dbcan-clean-disabled-route/" "$jobout/dbcan-route/"; fi
 if [ -d "$work/fs/home/leechuck/tmp/clean-route" ]; then rsync -a "$work/fs/home/leechuck/tmp/clean-route/" "$jobout/route/"; fi
 if [ -f "$work/fs/home/leechuck/Public/software/gspa/gspa-cli/build/libs/gspa-3.0.0.jar" ]; then cp "$work/fs/home/leechuck/Public/software/gspa/gspa-cli/build/libs/gspa-3.0.0.jar" "$jobout/gspa-3.0.0.jar"; fi
 printf '%s\n' "$code" > "$jobout/exit-code.txt"
 for module in gspa-core gspa-predictors gspa-cli; do
  if [ -d "$work/fs/home/leechuck/Public/software/gspa/$module/build/test-results" ]; then
   mkdir -p "$jobout/final-tests/$module"
   rsync -a "$work/fs/home/leechuck/Public/software/gspa/$module/build/test-results/" "$jobout/final-tests/$module/" || true
  fi
 done
 exit "$code"
}
trap finish EXIT
hostname > "$jobout/hostname.txt"
cp "$current/native-import-v2.sbatch" "$jobout/"
sha256sum "$package/before.tar.zst" "$package/python-runtime.tar.zst" "$current/source-v7.tar.gz" > "$jobout/packages.sha256"
zstd -dc "$package/before.tar.zst" | tar -xf - -C "$work/fs"
zstd -dc "$package/python-runtime.tar.zst" | tar -xf - -C "$work/fs"
tar -xf "$current/source-v7.tar.gz" -C "$work/fs"
repo=$work/fs/home/leechuck/Public/software/gspa
mkdir -p "$work/fs/home/leechuck/tmp"
container=(singularity exec --no-home --no-mount tmp --pwd /home/leechuck/Public/software/gspa
 --bind /usr,/lib,/lib64,/bin,/etc --bind /data:/data
 --bind "$work/fs/home/leechuck:/home/leechuck"
 --bind "$work/fs/usr/lib/jvm/java-21-openjdk-amd64:/usr/lib/jvm/java-21-openjdk-amd64"
 --env TMPDIR=/home/leechuck/tmp --env JAVA_TOOL_OPTIONS=-Djava.io.tmpdir=/home/leechuck/tmp
 --env JAVA_HOME=/usr/lib/jvm/java-21-openjdk-amd64
 --env PATH=/home/leechuck/tmp/gspa-v3-testenv/bin:/usr/lib/jvm/java-21-openjdk-amd64/bin:/usr/bin:/bin
 --env GSPA_TEST_PYTHON=/home/leechuck/tmp/gspa-v3-testenv/bin/python
 --env GRADLE_USER_HOME=/home/leechuck/.gradle
 /data/hohndor/gspa-v3/observed-gap-context-modes-v1-relocation/container-root)
mkdir -p "$repo/gspa-cli/build/libs"
cp "$current/jobs/10203/gspa-3.0.0.jar" "$repo/gspa-cli/build/libs/gspa-3.0.0.jar"
sha256sum "$current/jobs/10203/gspa-3.0.0.jar" "$repo/gspa-cli/build/libs/gspa-3.0.0.jar" > "$jobout/validated-jar.sha256"
cp "$current/audit_native.py" "$jobout/"
python3 "$jobout/audit_native.py" "$current/jobs/10201" --output "$jobout/native-audit.json" > "$jobout/native-audit.log" 2>&1
for dataset in new price; do
 "${container[@]}" java -jar "$current/jobs/10203/gspa-3.0.0.jar" clean-claims \
  --tsv "$current/jobs/10201/predictions/$dataset/$dataset.clean.tsv" \
  --receipt "$current/jobs/10201/predictions/$dataset/$dataset.clean.receipt.json" \
  --proteins "$current/jobs/10201/inputs/$dataset.fasta" --score-kind gmm_confidence \
  --output "$jobout/$dataset.claims.jsonl" > "$jobout/$dataset.import.log" 2>&1
done
python3 - "$current/jobs/10201" "$jobout" <<'CHECK' > "$jobout/import-parity.json"
import json,sys,hashlib,math
from pathlib import Path
native,out=map(Path,sys.argv[1:]);report={}
for dataset in ['new','price']:
 tsv=native/'predictions'/dataset/(dataset+'.clean.tsv')
 expected={}
 for line in tsv.read_text().splitlines()[1:]:
  pid,term,score,atype,metadata=line.split('\t');expected[pid,term]=(float(score),json.loads(metadata))
 claims=[json.loads(line) for line in (out/(dataset+'.claims.jsonl')).read_text().splitlines()]
 assert len(claims)==len(expected)
 seen=set()
 for row in claims:
  key=row['protein_id'],row['function_id'];assert key not in seen;seen.add(key)
  score,metadata=expected[key]
  assert row['function_type']=='EC' and row['source']=='clean' and row['raw_score']==score
  assert row['metadata']==metadata
  probability=row['calibrated_prob'];assert math.isfinite(probability) and 0<=probability<=1
 report[dataset]={'claims':len(claims),'proteins':len({r['protein_id'] for r in claims}),'all_scores_and_metadata_equal':True,'calibrated_prob_note':'Existing default ClaimExtractor heuristic; not native GMM output or learned calibration.','derived_probability_range':[min(r['calibrated_prob'] for r in claims),max(r['calibrated_prob'] for r in claims)]}
print(json.dumps({'status':'verified','datasets':report},indent=2))
CHECK
