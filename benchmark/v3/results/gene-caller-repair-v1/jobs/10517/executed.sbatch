#!/bin/bash
#SBATCH --job-name=gspa-gene-caller-repair
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
 if [ -f "$repo/gspa-cli/build/libs/gspa-3.0.0.jar" ]; then
  cp "$repo/gspa-cli/build/libs/gspa-3.0.0.jar" "$jobout/gspa-3.0.0.jar" || code=90
  sha256sum "$jobout/gspa-3.0.0.jar" > "$jobout/jar.sha256"
 fi
 if [ -d "$work/fs/home/leechuck/tmp/antismash-cli/results" ]; then rsync -a --exclude=nf-work --exclude=nf-missing-work --exclude=.nextflow "$work/fs/home/leechuck/tmp/antismash-cli/results/" "$jobout/cli/" || code=92; fi
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
sha256sum -c repair-packages-v3.sha256 > "$jobout/package-check.log"
hostname > "$jobout/hostname.txt"
cp "$0" "$jobout/executed.sbatch"
zstd -dc "$base/before.tar.zst" | tar -xf - -C "$work/fs"
zstd -dc "$base/python-runtime.tar.zst" | tar -xf - -C "$work/fs"
tar -xzf "$root/repair-source-v3.tar.gz" -C "$work/fs"
repo=$work/fs/home/leechuck/Public/software/gspa
/usr/bin/python3 - "$repo" "$root/repair-source-v3.json" <<'PY' > "$jobout/source-check.log"
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

# Reproduce against original production files; new test helpers do not change them.
tar -xzf "$root/before-source.tar.gz" -C "$work/fs"
mv "$repo/gspa-predictors/src/test/groovy/gspa/predictor/GeneCallerPipelineSpec.groovy" "$work/after-only-pipeline-spec.groovy"
set +e
"${container[@]}" ./gradlew :gspa-predictors:test --tests gspa.predictor.genecalling.GeneCallerContractSpec :gspa-cli:test --tests gspa.cli.AntiSmashClaimsSpec --continue --offline --no-daemon --max-workers=3 > "$jobout/before.log" 2>&1
before_code=$?
set -e
printf '%s\n' "$before_code" > "$jobout/before-command-exit.txt"
for module in gspa-predictors gspa-cli; do
 mkdir -p "$jobout/before-tests/$module"
 if [ -d "$repo/$module/build/test-results" ]; then rsync -a "$repo/$module/build/test-results/" "$jobout/before-tests/$module/"; fi
done
/usr/bin/python3 - "$jobout" <<'PYBEFORE' > "$jobout/before-gate.log"
from pathlib import Path
import sys,xml.etree.ElementTree as ET
p=Path(sys.argv[1]);assert int((p/'before-command-exit.txt').read_text())!=0
suites=[ET.parse(f).getroot() for f in (p/'before-tests').rglob('TEST-*.xml')]
assert suites,'No runtime test results: compilation failure is not a reproduced regression'
failures=[case.attrib['name'] for s in suites for case in s.findall('testcase') if case.find('failure') is not None or case.find('error') is not None]
assert any('both callers join reversed proteins' in n for n in failures),failures
assert any('malformed or incomplete output' in n for n in failures),failures
assert any('native numeric contigs' in n for n in failures),failures
print('Observed intended before failures:',len(failures));print('\n'.join(failures))
PYBEFORE
tar -xzf "$root/repair-source-v3.tar.gz" -C "$work/fs"
"${container[@]}" ./gradlew :gspa-predictors:test --tests gspa.predictor.genecalling.GeneCallerContractSpec --tests gspa.predictor.GeneCallerPipelineSpec --tests gspa.predictor.PipelineEmptyInputSpec :gspa-cli:test --tests gspa.cli.AntiSmashClaimsSpec :gspa-cli:shadowJar --offline --no-daemon --max-workers=3 > "$jobout/focused.log" 2>&1
for module in gspa-predictors gspa-cli; do
 mkdir -p "$jobout/focused-tests/$module"
 rsync -a "$repo/$module/build/test-results/" "$jobout/focused-tests/$module/"
done
"${container[@]}" ./gradlew :gspa-core:test :gspa-predictors:test :gspa-cli:test --continue --offline --no-daemon --max-workers=3 > "$jobout/full.log" 2>&1
