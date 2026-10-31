#!/bin/bash
#SBATCH --job-name=gspa-antismash-contract
#SBATCH --partition=debug
#SBATCH --cpus-per-task=3
#SBATCH --nodelist=node004
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#SBATCH --output=/data/hohndor/gspa-v3/antismash-observed-v1/slurm/%j.log
set -euo pipefail
base=/data/hohndor/gspa-v3/empty-input-v1
root=/data/hohndor/gspa-v3/antismash-observed-v1
work=/scratch/$SLURM_JOB_ID/antismash-contract
jobout=$root/jobs/$SLURM_JOB_ID
mkdir -p "$work/fs" "$jobout"
export TMPDIR=/scratch/$SLURM_JOB_ID/tmp
mkdir -p "$TMPDIR"
finish() {
 code=$?
 trap - EXIT
 set +e
 repo=$work/fs/home/leechuck/Public/software/gspa
 if [ -f "$repo/gspa-cli/build/libs/gspa-3.0.0.jar" ]; then
  cp "$repo/gspa-cli/build/libs/gspa-3.0.0.jar" "$jobout/gspa-3.0.0.jar" || code=90
  sha256sum "$jobout/gspa-3.0.0.jar" > "$jobout/jar.sha256"
 fi
 if [ -d "$work/fs/home/leechuck/tmp/antismash-cli/results" ]; then rsync -a "$work/fs/home/leechuck/tmp/antismash-cli/results/" "$jobout/cli/" || code=92; fi
 for module in gspa-core gspa-predictors gspa-cli; do
  if [ -d "$repo/$module/build/test-results" ]; then
   mkdir -p "$jobout/final-tests/$module"
   rsync -a "$repo/$module/build/test-results/" "$jobout/final-tests/$module/" || code=91
  fi
 done
 
 printf '%s\n' "$code" > "$jobout/exit-code.txt"
 date --iso-8601=seconds > "$jobout/finished-at.txt"
 exit "$code"
}
trap finish EXIT
cd "$root"
sha256sum -c build-packages-v2.sha256 > "$jobout/package-check.log"
hostname > "$jobout/hostname.txt"
cp "$0" "$jobout/executed.sbatch"
zstd -dc "$base/before.tar.zst" | tar -xf - -C "$work/fs"
zstd -dc "$base/python-runtime.tar.zst" | tar -xf - -C "$work/fs"
tar -xzf "$root/build-source-v2.tar.gz" -C "$work/fs"
repo=$work/fs/home/leechuck/Public/software/gspa
/usr/bin/python3 - "$repo" "$root/build-source-files-v2.json" <<'PY' > "$jobout/source-check.log"
import hashlib,json,sys
from pathlib import Path
root=Path(sys.argv[1]);manifest=json.loads(Path(sys.argv[2]).read_text())
for name,digest in manifest.items():
 with (root/name).open('rb') as stream:actual=hashlib.file_digest(stream,'sha256').hexdigest()
 assert actual==digest,name
print('verified',len(manifest),'source files')
PY
mkdir -p "$work/fs/home/leechuck/tmp"
container=(singularity exec --no-home --no-mount tmp --pwd /home/leechuck/Public/software/gspa
 --bind /usr,/lib,/lib64,/bin,/etc
 --bind "$work/fs/home/leechuck:/home/leechuck"
 --bind "$work/fs/usr/lib/jvm/java-21-openjdk-amd64:/usr/lib/jvm/java-21-openjdk-amd64"
 --env XDG_CACHE_HOME=/home/leechuck/tmp/cache --env PYTHONPATH=/home/leechuck/Public/software/gspa/scripts/tests:/home/leechuck/Public/software/gspa/benchmark/v3 --env TMPDIR=/home/leechuck/tmp --env JAVA_TOOL_OPTIONS=-Djava.io.tmpdir=/home/leechuck/tmp
 --env JAVA_HOME=/usr/lib/jvm/java-21-openjdk-amd64
 --env PATH=/home/leechuck/tmp/gspa-v3-testenv/bin:/usr/lib/jvm/java-21-openjdk-amd64/bin:/usr/bin:/bin
 --env GSPA_TEST_PYTHON=/home/leechuck/tmp/gspa-v3-testenv/bin/python
 --env GRADLE_USER_HOME=/home/leechuck/.gradle
 /data/hohndor/gspa-v3/observed-gap-context-modes-v1-relocation/container-root)
builder=$repo/gspa-predictors/src/main/groovy/gspa/predictor/specialized/AntiSmashPredictor.groovy
cp "$builder" "$work/after-builder.groovy"
cp "$root/before-AntiSmashPredictor.groovy" "$builder"
set +e
"${container[@]}" ./gradlew :gspa-predictors:test --tests gspa.predictor.specialized.AntiSmashObservedSpec --offline --no-daemon --max-workers=3 > "$jobout/before-focused.log" 2>&1
before=$?
set -e
printf '%s\n' "$before" > "$jobout/before-exit-code.txt"
mkdir -p "$jobout/before-tests"
rsync -a "$repo/gspa-predictors/build/test-results/" "$jobout/before-tests/"
test "$before" -ne 0
cp "$work/after-builder.groovy" "$builder"
"${container[@]}" ./gradlew :gspa-predictors:test --tests gspa.predictor.specialized.AntiSmashObservedSpec --tests gspa.predictor.AdapterRunContractSpec --tests gspa.predictor.AllPredictorParsingSpec :gspa-cli:shadowJar --offline --no-daemon --max-workers=3 > "$jobout/focused.log" 2>&1
mkdir -p "$jobout/focused-tests"
rsync -a "$repo/gspa-predictors/build/test-results/" "$jobout/focused-tests/"
"${container[@]}" ./gradlew :gspa-core:test :gspa-predictors:test :gspa-cli:test --continue --offline --no-daemon --max-workers=3 > "$jobout/full.log" 2>&1
mkdir -p "$work/fs/home/leechuck/tmp/antismash-cli"
cp /data/hohndor/gspa-v3/extracellular-builder-v1/jobs/10422/gspa-3.0.0.jar "$work/fs/home/leechuck/tmp/antismash-cli/before.jar"
cp "$root/GCA_003598085.1.json" "$work/fs/home/leechuck/tmp/antismash-cli/native.json"
"${container[@]}" python benchmark/v3/check_antismash_route.py --jar gspa-cli/build/libs/gspa-3.0.0.jar --before-jar /home/leechuck/tmp/antismash-cli/before.jar --native /home/leechuck/tmp/antismash-cli/native.json --java /usr/lib/jvm/java-21-openjdk-amd64/bin/java --out /home/leechuck/tmp/antismash-cli/results > "$jobout/cli.log" 2>&1
