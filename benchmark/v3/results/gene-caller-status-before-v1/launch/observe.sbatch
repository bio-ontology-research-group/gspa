#!/bin/bash
#SBATCH --job-name=gspa-gene-caller-outcomes
#SBATCH --partition=debug
#SBATCH --nodelist=node004
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --no-requeue
#SBATCH --output=/data/hohndor/gspa-v3/gene-caller-status-v1/slurm/%j.log
set -euo pipefail
root=/data/hohndor/gspa-v3/gene-caller-status-v1
jobout=$root/jobs/$SLURM_JOB_ID
mkdir -p "$jobout"
export TMPDIR=/scratch/$SLURM_JOB_ID/tmp
mkdir -p "$TMPDIR"
finish() { code=$?; trap - EXIT; printf '%s\n' "$code" > "$jobout/exit-code.txt"; date -Is > "$jobout/finished-at.txt"; exit "$code"; }
trap finish EXIT
cd "$root"
sha256sum -c package.sha256 > "$jobout/package-check.log"
cp "$0" "$jobout/executed.sbatch"
/usr/bin/python3 observe.py --out "$jobout/observed" --work /scratch/$SLURM_JOB_ID/observe > "$jobout/observe.log" 2>&1
