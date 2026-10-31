#!/bin/bash
#SBATCH --job-name=gspa-dbcan-contract
#SBATCH --partition=debug
#SBATCH --cpus-per-task=3
#SBATCH --nodelist=node002
#SBATCH --mem=16G
#SBATCH --time=00:40:00
#SBATCH --output=/data/hohndor/gspa-v3/dbcan-contract-v1/slurm/%j.log
set -euo pipefail
package=/data/hohndor/gspa-v3/empty-input-v1
current=/data/hohndor/gspa-v3/dbcan-contract-v1
work=/scratch/$SLURM_JOB_ID/gspa-dbcan
jobout=$current/jobs/$SLURM_JOB_ID
mkdir -p "$work/fs" "$jobout"
export TMPDIR=/scratch/$SLURM_JOB_ID/tmp
mkdir -p "$TMPDIR"
finish() {
 code=$?
 trap - EXIT
 printf '%s\n' "$code" > "$jobout/exit-code.txt"
 if [ -d "$work/fs/home/leechuck/Public/software/gspa/gspa-predictors/build/test-results" ]; then
  rsync -a "$work/fs/home/leechuck/Public/software/gspa/gspa-predictors/build/test-results/" "$jobout/final-tests/" || true
 fi
 exit "$code"
}
trap finish EXIT
hostname > "$jobout/hostname.txt"
cp "$current/validate-v2.sbatch" "$jobout/"
sha256sum "$package/before.tar.zst" "$package/python-runtime.tar.zst" "$current/source-v2.tar" > "$jobout/packages.sha256"
zstd -dc "$package/before.tar.zst" | tar -xf - -C "$work/fs"
zstd -dc "$package/python-runtime.tar.zst" | tar -xf - -C "$work/fs"
tar -xf "$current/source-v2.tar" -C "$work/fs"
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
adapter=$repo/gspa-predictors/src/main/groovy/gspa/predictor/specialized/DbCanPredictor.groovy
cp "$adapter" "$jobout/DbCanPredictor.after.groovy"
cp "$work/fs/home/leechuck/tmp/gspa-dbcan-contract-v1/DbCanPredictor.before.groovy" "$adapter"
set +e
"${container[@]}" ./gradlew :gspa-predictors:test --tests gspa.predictor.specialized.DbCanContractSpec --offline --no-daemon --max-workers=3 > "$jobout/before.log" 2>&1
before=$?
set -e
printf '%s\n' "$before" > "$jobout/before-exit-code.txt"
rsync -a "$repo/gspa-predictors/build/test-results/" "$jobout/before-tests/"
cp "$jobout/DbCanPredictor.after.groovy" "$adapter"
"${container[@]}" ./gradlew :gspa-predictors:test --tests gspa.predictor.specialized.DbCanContractSpec --offline --no-daemon --max-workers=3 > "$jobout/focused.log" 2>&1
rsync -a "$repo/gspa-predictors/build/test-results/" "$jobout/focused-tests/"
"${container[@]}" ./gradlew :gspa-predictors:test --offline --no-daemon --max-workers=3 > "$jobout/predictors.log" 2>&1
