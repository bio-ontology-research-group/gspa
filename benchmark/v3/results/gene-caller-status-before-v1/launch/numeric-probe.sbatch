#!/bin/bash
#SBATCH --job-name=gspa-numeric-reader
#SBATCH --partition=debug
#SBATCH --cpus-per-task=3
#SBATCH --nodelist=node004
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#SBATCH --output=/data/hohndor/gspa-v3/gene-caller-status-v1/slurm/%j.log
set -euo pipefail
base=/data/hohndor/gspa-v3/empty-input-v1
root=/data/hohndor/gspa-v3/gene-caller-status-v1
work=/scratch/$SLURM_JOB_ID/antismash-nextflow
jobout=$root/jobs/$SLURM_JOB_ID
mkdir -p "$work/fs" "$jobout"
export TMPDIR=/scratch/$SLURM_JOB_ID/tmp
mkdir -p "$TMPDIR"
finish() {
 code=$?
 trap - EXIT
 set +e
 repo=$work/fs/home/leechuck/Public/software/gspa

 if [ -d "$work/fs/home/leechuck/tmp/gene-probe/results" ]; then rsync -a --exclude=libraries "$work/fs/home/leechuck/tmp/gene-probe/results/" "$jobout/probe/" || code=92; fi
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
sha256sum -c numeric-probe-packages.sha256 > "$jobout/package-check.log"
hostname > "$jobout/hostname.txt"
cp "$0" "$jobout/executed.sbatch"
zstd -dc "$base/before.tar.zst" | tar -xf - -C "$work/fs"
zstd -dc "$base/python-runtime.tar.zst" | tar -xf - -C "$work/fs"
tar -xzf "/data/hohndor/gspa-v3/antismash-nextflow-v1/source-v2.tar.gz" -C "$work/fs"
repo=$work/fs/home/leechuck/Public/software/gspa
/usr/bin/python3 - "$repo" "/data/hohndor/gspa-v3/antismash-nextflow-v1/source-files-v2.json" <<'PY' > "$jobout/source-check.log"
import hashlib,json,sys
from pathlib import Path
root=Path(sys.argv[1]);manifest=json.loads(Path(sys.argv[2]).read_text())
for name,digest in manifest.items():
 with (root/name).open('rb') as stream:actual=hashlib.file_digest(stream,'sha256').hexdigest()
 assert actual==digest,name
print('verified',len(manifest),'source files')
PY
for module in gspa-core gspa-predictors gspa-cli; do
 rm -rf "$repo/$module/build/test-results"
done
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


probe=$work/fs/home/leechuck/tmp/gene-probe
mkdir -p "$probe"
cp "$root/numeric.groovy" "$probe/"
cp /data/hohndor/gspa-v3/antismash-nextflow-v1/jobs/10481/gspa-3.0.0.jar "$probe/before.jar"
test "$(sha256sum "$probe/before.jar" | cut -d' ' -f1)" = 069263f7b876222815928a637d02106b9c3b31a756edd1cf692432acc2e52303
cp -r "$root/jobs/10510/observed" "$probe/observed"
"${container[@]}" java -cp /home/leechuck/tmp/gene-probe/before.jar groovy.ui.GroovyMain /home/leechuck/tmp/gene-probe/numeric.groovy /home/leechuck/tmp/gene-probe/observed /home/leechuck/tmp/gene-probe/results > "$jobout/numeric.log" 2>&1
