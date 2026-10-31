#!/bin/bash
#SBATCH --job-name=gspa-clean-cohort-audit
#SBATCH --partition=debug
#SBATCH --nodelist=node003
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:40:00
#SBATCH --output=/data/hohndor/gspa-v3/clean-cohort-audit-v1/slurm/%j.log
set -euo pipefail
root=/data/hohndor/gspa-v3/clean-cohort-audit-v1
acquisition=/data/hohndor/gspa-v3/clean-cohort-v2
export TMPDIR=/scratch/$SLURM_JOB_ID
out=$root/jobs/$SLURM_JOB_ID
mkdir -p "$TMPDIR" "$out"
finish() { code=$?; trap - EXIT; printf '%s\n' "$code" > "$out/exit-code.txt"; exit "$code"; }
trap finish EXIT
cd "$root/source"
sha256sum -c "$root/source.sha256" > "$out/source-check.log"
python3 -m unittest discover -s benchmark/v3/tests -p 'test_*clean_cohort.py' -v > "$out/cohort-tests.log" 2>&1
python3 -m unittest discover -s benchmark/neural/tests -p test_clean_contract.py -v > "$out/adapter-tests.log" 2>&1
sha256sum -c "$acquisition/plan.sha256" > "$out/plan-check.log"
read -r plan_sha plan_file < "$acquisition/plan.sha256"
/usr/bin/time -v -o "$out/resources.txt" python3 benchmark/v3/audit_clean_cohort.py \
 --root "$acquisition/collection" --expected-plan "$plan_sha" --out "$out/audit.json"
