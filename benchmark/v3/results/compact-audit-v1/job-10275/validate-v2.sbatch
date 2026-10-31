#!/bin/bash
#SBATCH --job-name=gspa-compact-audit
#SBATCH --partition=debug
#SBATCH --cpus-per-task=3
#SBATCH --nodelist=node004
#SBATCH --mem=16G
#SBATCH --time=00:40:00
#SBATCH --output=/data/hohndor/gspa-v3/compact-audit-v1/slurm/%j.log
set -euo pipefail
package=/data/hohndor/gspa-v3/empty-input-v1
current=/data/hohndor/gspa-v3/compact-audit-v1
work=/scratch/$SLURM_JOB_ID/gspa-compact-audit
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
cp "$current/validate-v2.sbatch" "$jobout/"
sha256sum "$package/before.tar.zst" "$package/python-runtime.tar.zst" "$current/source-v2.tar.gz" > "$jobout/packages.sha256"
zstd -dc "$package/before.tar.zst" | tar -xf - -C "$work/fs"
zstd -dc "$package/python-runtime.tar.zst" | tar -xf - -C "$work/fs"
tar -xf "$current/source-v2.tar.gz" -C "$work/fs"
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
"${container[@]}" ./gradlew :gspa-core:test --tests gspa.integration.AuditChangesSpec --tests gspa.integration.IntegrationRunAuditSpec --offline --no-daemon --max-workers=3 > "$jobout/focused.log" 2>&1
mkdir -p "$jobout/focused-tests"
rsync -a "$repo/gspa-core/build/test-results/" "$jobout/focused-tests/"
"${container[@]}" ./gradlew :gspa-core:test :gspa-predictors:test :gspa-cli:test :gspa-cli:shadowJar --offline --no-daemon --max-workers=3 > "$jobout/full.log" 2>&1
oldjar=/data/hohndor/gspa-v3/price-production-v1/gspa-3.0.0.jar
newjar=$repo/gspa-cli/build/libs/gspa-3.0.0.jar
sha256sum "$oldjar" "$newjar" "$current/memory.groovy" "$current/parity.py" > "$jobout/comparison-inputs.sha256"
for version in before after; do
 if [ "$version" = before ]; then selectedjar=$oldjar; else selectedjar=$newjar; fi
 /usr/bin/time -v -o "$jobout/$version-memory.resources.txt" java -Xmx3g -Djava.io.tmpdir="$TMPDIR" -cp "$selectedjar" groovy.ui.GroovyMain "$current/memory.groovy" 1000000 > "$jobout/$version-memory.json" 2> "$jobout/$version-memory.err"
done
/data/hohndor/gspa-v3/price-production-v1/venv/bin/python "$current/parity.py" "$repo" "$newjar" "$jobout/parity" "$TMPDIR" > "$jobout/parity.log" 2>&1
