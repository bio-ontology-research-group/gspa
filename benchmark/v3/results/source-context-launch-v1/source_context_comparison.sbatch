#!/bin/bash
#SBATCH --job-name=gspa-source-context
#SBATCH --partition=debug
#SBATCH --nodelist=node004
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=2-00:00:00
#SBATCH --no-requeue
# Submit with --output=IMMUTABLE_BUNDLE/slurm/%j.log and arguments BUNDLE MODE.
set -euo pipefail
root=${1:?immutable bundle directory required}
mode=${2:?preflight or full required}
case "$mode" in preflight|full) ;; *) exit 2 ;; esac
if [ "$mode" = full ] && [ "$#" -ne 4 ]; then
 echo 'full mode requires explicit JAR and verified probe receipt paths' >&2
 exit 2
fi
jar=${3:-/data/hohndor/gspa-v3/price-production-v3/gspa-3.0.0.jar}
probe=${4:-/data/hohndor/gspa-v3/multisource-comparison-probe-v3/jobs/10389/results/verification.json}
base=/data/hohndor/gspa-v3/observed-gap-context-modes-v1-relocation
supplement=/data/hohndor/gspa-v3/context-resume-v1
work=/scratch/$SLURM_JOB_ID/source-context
jobout=$root/jobs/$SLURM_JOB_ID
logical_source=/home/leechuck/tmp/source-context-current
logical_results=/home/leechuck/tmp/source-context-results
physical_results=$work/fs$logical_results
mkdir -p "$work/fs" "$jobout"
export TMPDIR=/scratch/$SLURM_JOB_ID/tmp
mkdir -p "$TMPDIR"
finish() {
 code=$?
 trap - EXIT
 set +e
 if [ -d "$physical_results" ]; then
  rsync -a "$physical_results/" "$jobout/results/" || code=90
 fi
 printf '%s\n' "$code" > "$jobout/exit-code.txt"
 date --iso-8601=seconds > "$jobout/finished-at.txt"
 exit "$code"
}
trap finish EXIT
cd "$root"
sha256sum -c packages.sha256 > "$jobout/package-check.log"
hostname > "$jobout/hostname.txt"
date --iso-8601=seconds > "$jobout/started-at.txt"
df -h "$work" "$jobout" > "$jobout/storage-before.txt"
cp "$root/packages.sha256" "$root/source.sha256" "$jobout/"
cp "$0" "$jobout/executed.sbatch"
zstd -dc "$base/package.tar.zst" | tar -xf - -C "$work/fs"
tar -xf "$supplement/supplement.tar" -C "$work/fs"
tar -xf "$root/parent-supplement.tar" -C "$work/fs"
mkdir -p "$work/fs$logical_source" "$work/fs/home/leechuck/tmp"
tar -xzf "$root/source.tar.gz" -C "$work/fs$logical_source"
cd "$work/fs$logical_source"
sha256sum -c "$root/source.sha256" > "$jobout/source-check.log"
container=(singularity exec --no-home --no-mount tmp --pwd "$logical_source"
 --bind /usr,/lib,/lib64,/bin,/etc
 --bind /data:/data:ro --bind "$jobout:$jobout"
 --bind "$work/fs/home/leechuck:/home/leechuck"
 --bind "$work/fs/usr/lib/jvm/java-21-openjdk-amd64:/usr/lib/jvm/java-21-openjdk-amd64"
 --env TMPDIR=/home/leechuck/tmp
 --env JAVA_TOOL_OPTIONS=-Djava.io.tmpdir=/home/leechuck/tmp
 --env PATH=/usr/lib/jvm/java-21-openjdk-amd64/bin:/usr/bin:/bin
 --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1
 --env PYTHONPATH="$logical_source/vendor" --env IJSON_BACKEND=python
 "$base/container-root")
python=/home/leechuck/tmp/gspa-v3-testenv/bin/python
"${container[@]}" java -version > "$jobout/java-version.txt" 2>&1
"${container[@]}" rsync --version > "$jobout/rsync-version.txt" 2>&1
"${container[@]}" "$python" -c 'import sys,numpy,ijson; print(sys.version); print(numpy.__version__,ijson.__version__,ijson.backend)' > "$jobout/python-version.txt" 2>&1
"${container[@]}" "$python" - <<'PY' > "$jobout/tests.log" 2>&1
import unittest
suite = unittest.TestSuite()
for pattern, expected in (('test_source_arm_context*.py',37), ('test_verified_batch.py',21),
                         ('test_streaming_ledger.py',11), ('test_ec_reference_panel.py',5),
                         ('test_observed_gap_context_analysis.py',4), ('test_source_context_remote.py',5),
                         ('test_cohort_sequence_inputs.py',22)):
    tests = unittest.defaultTestLoader.discover('benchmark/v3/tests', pattern=pattern)
    assert tests.countTestCases() == expected, (pattern, tests.countTestCases(), expected)
    suite.addTests(tests)
result = unittest.TextTestRunner(verbosity=2).run(suite)
assert result.wasSuccessful() and not result.skipped
PY
/usr/bin/time -v -o "$jobout/resources.txt" "${container[@]}" "$python" "$logical_source/benchmark/v3/run_source_context_remote.py" "$mode" \
 --out "$logical_results" --logs "$jobout/execution" \
 --bindings "$root/packages.sha256" \
 --jar "$jar" \
 --clean-root /data/hohndor/gspa-v3/clean-cohort-v3/collection \
 --clean-audit /data/hohndor/gspa-v3/clean-cohort-v3/jobs/10368/audit.json \
 --probe "$probe" \
 > "$jobout/execution.log" 2>&1
