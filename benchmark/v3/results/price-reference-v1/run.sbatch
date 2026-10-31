#!/bin/bash
#SBATCH --job-name=gspa-price-reference
#SBATCH --partition=debug
#SBATCH --nodelist=node004
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/data/hohndor/gspa-v3/price-reference-v1/slurm/%j.log
set -euo pipefail
root=/data/hohndor/gspa-v3/price-reference-v1
export TMPDIR=/scratch/$SLURM_JOB_ID
mkdir -p "$TMPDIR" "$root/jobs/$SLURM_JOB_ID"
cd "$root/source"
sha256sum -c "$root/source.sha256" > "$root/jobs/$SLURM_JOB_ID/source-check.log"
python3 -m unittest discover -s benchmark/v3/tests -p test_price_reference_metrics.py -v > "$root/jobs/$SLURM_JOB_ID/tests.log" 2>&1
set +e
/usr/bin/time -v -o "$root/jobs/$SLURM_JOB_ID/resources.txt" python3 benchmark/v3/evaluate_price_reference.py \
  --collection /data/hohndor/gspa-v3/price-production-v1/collection \
  --manifest benchmark/manifests/suite6_price149_manifest.tsv \
  --ec-audit "$root/ec-audit.json" --predictors topec mapred \
  --out "$root/jobs/$SLURM_JOB_ID/evaluation"
code=$?
printf '%s\n' "$code" > "$root/jobs/$SLURM_JOB_ID/exit-code.txt"
exit "$code"
