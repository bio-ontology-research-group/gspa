#!/bin/bash
#SBATCH --job-name=gspa-dbcan-route
#SBATCH --partition=debug
#SBATCH --cpus-per-task=3
#SBATCH --nodelist=node002
#SBATCH --mem=16G
#SBATCH --time=00:40:00
#SBATCH --output=/data/hohndor/gspa-v3/dbcan-route-v1/slurm/%j.log
set -euo pipefail
package=/data/hohndor/gspa-v3/empty-input-v1
current=/data/hohndor/gspa-v3/dbcan-route-v1
work=/scratch/$SLURM_JOB_ID/gspa-dbcan-route
jobout=$current/jobs/$SLURM_JOB_ID
mkdir -p "$work/fs" "$jobout"
export TMPDIR=/scratch/$SLURM_JOB_ID/tmp
mkdir -p "$TMPDIR"
finish() {
 code=$?
 trap - EXIT
 if [ -f "$work/fs/home/leechuck/Public/software/gspa/gspa-cli/build/libs/gspa-3.0.0.jar" ]; then cp "$work/fs/home/leechuck/Public/software/gspa/gspa-cli/build/libs/gspa-3.0.0.jar" "$jobout/gspa-3.0.0.jar"; fi
 printf '%s\n' "$code" > "$jobout/exit-code.txt"
 for module in gspa-core gspa-predictors gspa-cli; do
  if [ -d "$work/fs/home/leechuck/Public/software/gspa/$module/build/test-results" ]; then
   mkdir -p "$jobout/final-tests/$module"
   rsync -a "$work/fs/home/leechuck/Public/software/gspa/$module/build/test-results/" "$jobout/final-tests/$module/" || true
  fi
 done
 if [ -d "$work/fs/home/leechuck/tmp/gspa-dbcan-route-v1/replay" ]; then
  rsync -a --exclude=nf-work --exclude=nf-missing-work --exclude=.nextflow "$work/fs/home/leechuck/tmp/gspa-dbcan-route-v1/replay/" "$jobout/replay/" || true
 fi
 exit "$code"
}
trap finish EXIT
hostname > "$jobout/hostname.txt"
cp "$current/validate-v5.sbatch" "$jobout/"
sha256sum "$package/before.tar.zst" "$package/python-runtime.tar.zst" "$current/source-v4.tar" > "$jobout/packages.sha256"
zstd -dc "$package/before.tar.zst" | tar -xf - -C "$work/fs"
zstd -dc "$package/python-runtime.tar.zst" | tar -xf - -C "$work/fs"
tar -xf "$current/source-v4.tar" -C "$work/fs"
repo=$work/fs/home/leechuck/Public/software/gspa
mkdir -p "$work/fs/home/leechuck/tmp"
container=(singularity exec --no-home --no-mount tmp --pwd /home/leechuck/Public/software/gspa
 --bind /usr,/lib,/lib64,/bin,/etc
 --bind "$work/fs/home/leechuck:/home/leechuck"
 --bind "$work/fs/usr/lib/jvm/java-21-openjdk-amd64:/usr/lib/jvm/java-21-openjdk-amd64"
 --env TMPDIR=/home/leechuck/tmp --env JAVA_TOOL_OPTIONS=-Djava.io.tmpdir=/home/leechuck/tmp
 --env JAVA_HOME=/usr/lib/jvm/java-21-openjdk-amd64
 --env PATH=/home/leechuck/tmp/gspa-v3-testenv/bin:/usr/lib/jvm/java-21-openjdk-amd64/bin:/usr/bin:/bin
 --env GSPA_TEST_PYTHON=/home/leechuck/tmp/gspa-v3-testenv/bin/python
 --env GRADLE_USER_HOME=/home/leechuck/.gradle
 /data/hohndor/gspa-v3/observed-gap-context-modes-v1-relocation/container-root)
# Save the current source overlay and reproduce absence of route support.
python3 - "$work/fs" <<'PYBEFORE'
from pathlib import Path
import shutil,sys
fs=Path(sys.argv[1]);base=fs/'home/leechuck/tmp/gspa-dbcan-route-v1/before'
repo=fs/'home/leechuck/Public/software/gspa'
for p in base.rglob('*'):
 if p.is_file():shutil.copyfile(p,repo/p.relative_to(base))
PYBEFORE
set +e
"${container[@]}" ./gradlew :gspa-predictors:test --tests gspa.predictor.specialized.DbCanRouteSpec :gspa-cli:test --tests gspa.cli.DbCanClaimsSpec --continue --offline --no-daemon --max-workers=3 > "$jobout/before.log" 2>&1
before=$?
set -e
printf '%s\n' "$before" > "$jobout/before-exit-code.txt"
for module in gspa-predictors gspa-cli; do
 mkdir -p "$jobout/before-tests/$module"
 if [ -d "$repo/$module/build/test-results" ]; then rsync -a "$repo/$module/build/test-results/" "$jobout/before-tests/$module/"; fi
done
tar -xf "$current/source-v4.tar" -C "$work/fs"
"${container[@]}" ./gradlew :gspa-predictors:test --tests gspa.predictor.specialized.DbCanRouteSpec --tests gspa.predictor.specialized.DbCanContractSpec :gspa-cli:test --tests gspa.cli.DbCanClaimsSpec --offline --no-daemon --max-workers=3 > "$jobout/focused.log" 2>&1
for module in gspa-predictors gspa-cli; do
 mkdir -p "$jobout/focused-tests/$module"
 rsync -a "$repo/$module/build/test-results/" "$jobout/focused-tests/$module/"
done
"${container[@]}" ./gradlew :gspa-predictors:test :gspa-cli:test :gspa-cli:shadowJar --offline --no-daemon --max-workers=3 > "$jobout/full.log" 2>&1
"${container[@]}" python /home/leechuck/Public/software/gspa/benchmark/v3/check_dbcan_route.py --jar /home/leechuck/Public/software/gspa/gspa-cli/build/libs/gspa-3.0.0.jar --nextflow /home/leechuck/bin/nextflow --out /home/leechuck/tmp/gspa-dbcan-route-v1/replay > "$jobout/routes.log" 2>&1
