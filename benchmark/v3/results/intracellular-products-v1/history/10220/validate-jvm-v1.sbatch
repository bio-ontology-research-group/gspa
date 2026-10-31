#!/bin/bash
#SBATCH --job-name=gspa-intracellular-products
#SBATCH --partition=debug
#SBATCH --cpus-per-task=3
#SBATCH --nodelist=node002
#SBATCH --mem=16G
#SBATCH --time=00:40:00
#SBATCH --output=/data/hohndor/gspa-v3/intracellular-products-v1/slurm/%j.log
set -euo pipefail
package=/data/hohndor/gspa-v3/empty-input-v1
current=/data/hohndor/gspa-v3/intracellular-products-v1
work=/scratch/$SLURM_JOB_ID/gspa-intracellular-products
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
 exit "$code"
}
trap finish EXIT
hostname > "$jobout/hostname.txt"
cp "$current/validate-jvm-v1.sbatch" "$jobout/"
sha256sum "$package/before.tar.zst" "$package/python-runtime.tar.zst" "$current/source-v1.tar.gz" > "$jobout/packages.sha256"
zstd -dc "$package/before.tar.zst" | tar -xf - -C "$work/fs"
zstd -dc "$package/python-runtime.tar.zst" | tar -xf - -C "$work/fs"
tar -xf "$current/source-v1.tar.gz" -C "$work/fs"
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
"${container[@]}" python -m unittest discover -s benchmark/neural/tests -p test_clean_contract.py -v > "$jobout/python-tests.log" 2>&1
"${container[@]}" ./gradlew :gspa-core:test --tests gspa.integration.workflow.MaskedAssignmentBuilderSpec --tests gspa.integration.workflow.ExperimentBiologySpec --offline --no-daemon --max-workers=3 > "$jobout/focused.log" 2>&1
mkdir -p "$jobout/focused-tests"
rsync -a "$repo/gspa-core/build/test-results/" "$jobout/focused-tests/"
"${container[@]}" ./gradlew :gspa-core:test :gspa-predictors:test :gspa-cli:test :gspa-cli:shadowJar --offline --no-daemon --max-workers=3 > "$jobout/full.log" 2>&1
