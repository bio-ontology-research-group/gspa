#!/bin/bash
#SBATCH --job-name=gspa-source-probe
#SBATCH --partition=debug
#SBATCH --nodelist=node002
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --no-requeue
#SBATCH --output=/data/hohndor/gspa-v3/multisource-comparison-probe-v4/slurm/%j.log
set -euo pipefail
root=/data/hohndor/gspa-v3/multisource-comparison-probe-v4
base=/data/hohndor/gspa-v3/observed-gap-context-modes-v1-relocation
work=/scratch/$SLURM_JOB_ID/source-probe
jobout=$root/jobs/$SLURM_JOB_ID
mkdir -p "$work/fs" "$work/source" "$jobout"
export TMPDIR=/scratch/$SLURM_JOB_ID/tmp
export JAVA_TOOL_OPTIONS="-Djava.io.tmpdir=$TMPDIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export IJSON_BACKEND=python
export PYTHONPATH=/data/hohndor/gspa-v3/price-production-v1/venv/lib/python3.12/site-packages
mkdir -p "$TMPDIR"
finish() {
 code=$?
 trap - EXIT
 set +e
 if [ -d "$work/results" ]; then rsync -a "$work/results/" "$jobout/results/" || code=90; fi
 printf '%s\n' "$code" > "$jobout/exit-code.txt"
 date --iso-8601=seconds > "$jobout/finished-at.txt"
 exit "$code"
}
trap finish EXIT
cd "$root"
sha256sum -c packages.sha256 > "$jobout/package-check.log"
hostname > "$jobout/hostname.txt"
zstd -dc "$base/package.tar.zst" | tar -xf - -C "$work/fs"
tar -xzf "$root/source.tar.gz" -C "$work/source"
cd "$work/source"
sha256sum -c "$root/source.sha256" > "$jobout/source-check.log"
python=/data/hohndor/gspa-v3/clean-cohort-v2/env/bin/python
"$python" -m unittest discover -s benchmark/v3/tests -p 'test_source_arm_context*.py' -v > "$jobout/tests.log" 2>&1
"$python" -m unittest discover -s benchmark/v3/tests -p test_cohort_sequence_inputs.py -v > "$jobout/source-input-tests.log" 2>&1
java -version > "$jobout/java-version.txt" 2>&1
"$python" -c 'import sys,numpy,scipy,ijson; print(sys.version); print(numpy.__version__,scipy.__version__,ijson.__version__,ijson.backend)' > "$jobout/python-version.txt"
options=()
if [ "${1:-}" = clean ]; then options+=(--clean); fi
/usr/bin/time -v -o "$jobout/resources.txt" "$python" "$root/probe.py" --source "$work/source" --fs "$work/fs" \
 --out "$work/results" --jar /data/hohndor/gspa-v3/ec-prefix-identity-v3/jobs/10402/gspa-3.0.0.jar "${options[@]}" > "$jobout/probe.log" 2>&1
