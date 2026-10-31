#!/bin/bash
#SBATCH --job-name=gspa-clean-cohort
#SBATCH --partition=debug
#SBATCH --nodelist=node004
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=08:00:00
#SBATCH --output=/data/hohndor/gspa-v3/clean-cohort-v2/slurm/%A_%a.log
set -euo pipefail
root=/data/hohndor/gspa-v3/clean-cohort-v2
export TMPDIR=/scratch/$SLURM_JOB_ID
mkdir -p "$TMPDIR" "$root/job-status/$SLURM_JOB_ID"
finish() { code=$?; trap - EXIT; printf '%s\n' "$code" > "$root/job-status/$SLURM_JOB_ID/exit-code.txt"; exit "$code"; }
trap finish EXIT
cd "$root/source"
sha256sum -c "$root/source.sha256" > "$root/job-status/$SLURM_JOB_ID/source-check.log"
sha256sum -c "$root/plan.sha256" > "$root/job-status/$SLURM_JOB_ID/plan-check.log"
read -r plan_sha plan_file < "$root/plan.sha256"
nvidia-smi > "$root/job-status/$SLURM_JOB_ID/nvidia-smi.txt"
hostname > "$root/job-status/$SLURM_JOB_ID/hostname.txt"
/usr/bin/time -v -o "$root/job-status/$SLURM_JOB_ID/resources.txt" python3 benchmark/v3/observe_clean_cohort.py run \
 --out "$root/collection" --expected-plan "$plan_sha" --index "$SLURM_ARRAY_TASK_ID" --scratch "$TMPDIR"
