#!/bin/bash
#SBATCH --job-name=gspa-probe-independent
#SBATCH --partition=debug
#SBATCH --nodelist=node004
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#SBATCH --output=/data/hohndor/gspa-v3/source-probe-audit-v2/slurm/%j.log
set -euo pipefail
root=/data/hohndor/gspa-v3/source-probe-audit-v2
out=$root/jobs/$SLURM_JOB_ID
mkdir -p "$out"
export TMPDIR=/scratch/$SLURM_JOB_ID
finish() { code=$?; trap - EXIT; printf '%s\n' "$code" > "$out/exit-code.txt"; exit "$code"; }
trap finish EXIT
cd "$root"
sha256sum -c source.sha256 > "$out/source-check.log"
cp "$0" "$out/executed.sbatch"
/usr/bin/time -v -o "$out/resources.txt" python3 "$root/audit.py" > "$out/audit.log" 2>&1
