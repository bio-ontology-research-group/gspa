#!/bin/bash
#SBATCH --job-name=gspa-price-full-reference
#SBATCH --partition=debug
#SBATCH --nodelist=node004
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=3-00:00:00
#SBATCH --no-requeue
#SBATCH --output=/data/hohndor/gspa-v3/price-reference-full-v1/slurm/%j.log
set -euo pipefail
root=/data/hohndor/gspa-v3/price-reference-full-v1
jobout=$root/jobs/$SLURM_JOB_ID
export TMPDIR=/scratch/$SLURM_JOB_ID
export PYTHONDONTWRITEBYTECODE=1
mkdir -p "$TMPDIR" "$jobout"
finish() {
 code=$?
 trap - EXIT
 printf '%s\n' "$code" > "$jobout/exit-code.txt"
 date --iso-8601=seconds > "$jobout/finished-at.txt"
 exit "$code"
}
trap finish EXIT
cd "$root"
sha256sum -c package.sha256 > "$jobout/package-check.log"
cp "$0" "$jobout/executed.sbatch"
python3 prerequisites.py "$jobout/prerequisites.json"
cd "$root/source"
sha256sum -c "$root/source.sha256" > "$jobout/source-check.log"
python3 -m unittest discover -s benchmark/v3/tests -p test_price_reference_metrics.py -v > "$jobout/tests.log" 2>&1
/usr/bin/time -v -o "$jobout/resources.txt" python3 benchmark/v3/evaluate_price_reference.py \
 --collection /data/hohndor/gspa-v3/price-production-v3/collection \
 --manifest benchmark/manifests/suite6_price149_manifest.tsv \
 --ec-audit "$root/ec-audit.json" --predictors topec mapred enzbert dpz graphec clean \
 --out "$jobout/evaluation" > "$jobout/evaluation.log" 2>&1
python3 "$root/check_retained.py" --root "$jobout/evaluation" \
 --audit "$root/ec-audit.json" --manifest benchmark/manifests/suite6_price149_manifest.tsv \
 > "$jobout/independent-arithmetic.json"
