#!/bin/bash
#SBATCH --job-name=gspa-context-audit
#SBATCH --partition=debug
#SBATCH --nodelist=node002
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/data/hohndor/gspa-v3/context-comparison-audit-v1/slurm/%j.log
set -euo pipefail
current=/data/hohndor/gspa-v3/context-comparison-audit-v1
package=/data/hohndor/gspa-v3/observed-gap-context-modes-v1-relocation
resume=/data/hohndor/gspa-v3/context-resume-v1
work=/scratch/$SLURM_JOB_ID/context-audit
jobout=$current/jobs/$SLURM_JOB_ID
mkdir -p "$work/fs" "$jobout"
export TMPDIR=/scratch/$SLURM_JOB_ID/tmp
mkdir -p "$TMPDIR"
finish() { code=$?; trap - EXIT; printf '%s\n' "$code" > "$jobout/exit-code.txt"; date --iso-8601=seconds > "$jobout/finished-at.txt"; exit "$code"; }
trap finish EXIT
hostname > "$jobout/hostname.txt"
python3 --version > "$jobout/python-version.txt"
date --iso-8601=seconds > "$jobout/started-at.txt"
cp "$current/audit-v1.sbatch" "$current/audit-v1.py" "$jobout/"
sha256sum "$package/package.tar.zst" "$resume/supplement.tar" "$resume/supplement.json" "$jobout/audit-v1.py" > "$jobout/inputs.sha256"
zstd -dc "$package/package.tar.zst" | tar -xf - -C "$work/fs"
tar -xf "$resume/supplement.tar" -C "$work/fs"
python3 "$jobout/audit-v1.py" --package "$work/fs" --base "$package/jobs/10163/results" --final "$resume/jobs/10176/results" --supplement "$resume/supplement.json" --out "$jobout/audit" --workers 4 > "$jobout/audit.log" 2>&1
