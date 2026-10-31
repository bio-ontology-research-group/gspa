#!/bin/bash
#SBATCH --job-name=gspa-antismash-cli
#SBATCH --partition=debug
#SBATCH --cpus-per-task=3
#SBATCH --nodelist=node004
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --no-requeue
#SBATCH --output=/data/hohndor/gspa-v3/antismash-observed-v1/slurm/%j.log
set -euo pipefail
base=/data/hohndor/gspa-v3/empty-input-v1
root=/data/hohndor/gspa-v3/antismash-observed-v1
work=/scratch/$SLURM_JOB_ID/antismash-cli
jobout=$root/jobs/$SLURM_JOB_ID
mkdir -p "$work/fs" "$jobout"
export TMPDIR=/scratch/$SLURM_JOB_ID/tmp
mkdir -p "$TMPDIR"
finish() {
 code=$?
 trap - EXIT
 set +e
 if [ -d "$work/fs/home/leechuck/tmp/antismash-cli/results" ]; then rsync -a "$work/fs/home/leechuck/tmp/antismash-cli/results/" "$jobout/cli/" || code=92; fi
 printf '%s\n' "$code" > "$jobout/exit-code.txt"
 date --iso-8601=seconds > "$jobout/finished-at.txt"
 exit "$code"
}
trap finish EXIT
cd "$root"
sha256sum -c cli-packages-v1.sha256 > "$jobout/package-check.log"
cp "$0" "$jobout/executed.sbatch"
python3 - <<'CHECK' > "$jobout/test-prerequisites.json"
from pathlib import Path
import xml.etree.ElementTree as E
import json
root=Path('/data/hohndor/gspa-v3/antismash-observed-v1/jobs/10439/final-tests')
sums={k:0 for k in ('tests','failures','errors','skipped')}
for f in root.rglob('TEST-*.xml'):
 x=E.parse(f).getroot()
 for k in sums:sums[k]+=int(x.get(k,'0'))
assert sums==dict(tests=1502,failures=0,errors=0,skipped=6),sums
print(json.dumps(sums,indent=2))
CHECK
zstd -dc "$base/before.tar.zst" | tar -xf - -C "$work/fs"
zstd -dc "$base/python-runtime.tar.zst" | tar -xf - -C "$work/fs"
mkdir -p "$work/fs/home/leechuck/tmp/antismash-cli"
cp "$root/check_antismash_route.py" "$root/check_native_memberships.py" "$work/fs/home/leechuck/tmp/antismash-cli/"
cp "$root/jobs/10439/gspa-3.0.0.jar" "$work/fs/home/leechuck/tmp/antismash-cli/after.jar"
cp /data/hohndor/gspa-v3/extracellular-builder-v1/jobs/10422/gspa-3.0.0.jar "$work/fs/home/leechuck/tmp/antismash-cli/before.jar"
cp "$root/GCA_003598085.1.json" "$work/fs/home/leechuck/tmp/antismash-cli/native.json"
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

"${container[@]}" python /home/leechuck/tmp/antismash-cli/check_antismash_route.py --jar /home/leechuck/tmp/antismash-cli/after.jar --before-jar /home/leechuck/tmp/antismash-cli/before.jar --native /home/leechuck/tmp/antismash-cli/native.json --java /usr/lib/jvm/java-21-openjdk-amd64/bin/java --out /home/leechuck/tmp/antismash-cli/results > "$jobout/cli.log" 2>&1
"${container[@]}" python /home/leechuck/tmp/antismash-cli/check_native_memberships.py --native /home/leechuck/tmp/antismash-cli/native.json --cli /home/leechuck/tmp/antismash-cli/results > "$jobout/independent-memberships.json"
