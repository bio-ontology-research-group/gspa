#!/bin/bash
#SBATCH --job-name=gspa-pyrodigal-operons
#SBATCH --partition=debug
#SBATCH --cpus-per-task=3
#SBATCH --nodelist=node004
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#SBATCH --output=/data/hohndor/gspa-v3/operon-identity-v1/slurm/%j.log
set -euo pipefail
base=/data/hohndor/gspa-v3/empty-input-v1
root=/data/hohndor/gspa-v3/operon-identity-v1
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
 if [ -d "$work/fs/home/leechuck/tmp/antismash-cli/results" ]; then rsync -a --exclude=single-work --exclude=meta-work --exclude=libraries --exclude=.nextflow "$work/fs/home/leechuck/tmp/antismash-cli/results/" "$jobout/cli/" || code=92; fi
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
sha256sum -c route-packages-v1.sha256 > "$jobout/package-check.log"
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
/usr/bin/python3 - "/data/hohndor/gspa-v3/antismash-nextflow-v1/jobs/10481" <<'PYGATE' > "$jobout/parent-gate.log"
import hashlib,sys,xml.etree.ElementTree as ET
from pathlib import Path
p=Path(sys.argv[1]); suites=[ET.parse(f).getroot() for f in (p/'final-tests').rglob('TEST-*.xml')]
assert sum(int(s.attrib['tests']) for s in suites)==1521
assert sum(int(s.attrib['skipped']) for s in suites)==6
assert not any(int(s.attrib['failures']) or int(s.attrib['errors']) for s in suites)
assert 'BUILD SUCCESSFUL' in (p/'full.log').read_text()
assert hashlib.sha256((p/'gspa-3.0.0.jar').read_bytes()).hexdigest()=='069263f7b876222815928a637d02106b9c3b31a756edd1cf692432acc2e52303'
print('Parent10481:1521 cases,6 skips,zero failures/errors; exact JAR verified. Parent replay launcher failure preserved.')
PYGATE
mkdir -p "$repo/gspa-cli/build/libs"
cp "/data/hohndor/gspa-v3/antismash-nextflow-v1/jobs/10481/gspa-3.0.0.jar" "$repo/gspa-cli/build/libs/gspa-3.0.0.jar"
tar -xzf "$root/route-overlay-v1.tar.gz" -C "$work/fs"
/usr/bin/python3 - "$repo" "$root/route-source-v1.json" <<'PYSOURCE' > "$jobout/overlay-source-check.log"
import hashlib,json,sys
from pathlib import Path
root=Path(sys.argv[1]);manifest=json.loads(Path(sys.argv[2]).read_text())
for name,digest in manifest.items():assert hashlib.sha256((root/name).read_bytes()).hexdigest()==digest,name
print('All route-overlay sources verified')
PYSOURCE
mkdir -p "$work/fs/home/leechuck/tmp/native-operon/observed" "$work/fs/home/leechuck/tmp/native-operon/wheels"
cp "$root/observed/10484/genome.fna" "$root/observed/10484/input.json" "$root/observed/10484/identity-observation.json" "$root/observed/10484/single.gff" "$root/observed/10484/single.faa" "$root/observed/10484/meta.gff" "$root/observed/10484/meta.faa" "$work/fs/home/leechuck/tmp/native-operon/observed/"
cp -r "$root/observed/10484/single-operons-before" "$root/observed/10484/meta-operons-before" "$work/fs/home/leechuck/tmp/native-operon/observed/"
cp "$root"/wheels/*.whl "$work/fs/home/leechuck/tmp/native-operon/wheels/"
chmod u+x "$work/fs/home/leechuck/bin/nextflow"
"${container[@]}" python benchmark/v3/check_pyrodigal_operon_route.py --jar gspa-cli/build/libs/gspa-3.0.0.jar --nextflow /home/leechuck/bin/nextflow --observed /home/leechuck/tmp/native-operon/observed --wheels /home/leechuck/tmp/native-operon/wheels --out /home/leechuck/tmp/antismash-cli/results > "$jobout/routes.log" 2>&1
