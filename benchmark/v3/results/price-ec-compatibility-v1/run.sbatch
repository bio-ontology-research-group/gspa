#!/bin/bash
#SBATCH --job-name=gspa-ec-price
#SBATCH --partition=debug
#SBATCH --nodelist=node004
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#SBATCH --output=/data/hohndor/gspa-v3/price-ec-compatibility-v1/slurm/%j.log
set -euo pipefail
root=/data/hohndor/gspa-v3/price-ec-compatibility-v1
jobout=$root/jobs/$SLURM_JOB_ID
work=/scratch/$SLURM_JOB_ID/price-ec
mkdir -p "$jobout" "$work/results" "$work/tmp"
export TMPDIR=$work/tmp
export IJSON_BACKEND=python
export PYTHONPATH=/data/hohndor/gspa-v3/price-production-v1/venv/lib/python3.12/site-packages
unset JAVA_TOOL_OPTIONS JDK_JAVA_OPTIONS _JAVA_OPTIONS
finish() {
 code=$?
 trap - EXIT
 set +e
 rsync -a "$work/results/" "$jobout/results/" || code=90
 printf '%s\n' "$code" > "$jobout/exit-code.txt"
 date --iso-8601=seconds > "$jobout/finished-at.txt"
 exit "$code"
}
trap finish EXIT
cd "$root"
sha256sum -c packages.sha256 > "$jobout/package-check.log"
cp "$0" "$jobout/executed.sbatch"
/usr/bin/time -v -o "$jobout/resources.txt" python3 "$root/check.py" --out "$work/results" > "$jobout/run.log" 2>&1
