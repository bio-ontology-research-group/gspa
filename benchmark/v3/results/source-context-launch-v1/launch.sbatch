#!/bin/bash
#SBATCH --job-name=gspa-full-source-context
#SBATCH --partition=debug
#SBATCH --nodelist=node004
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=2-00:00:00
#SBATCH --no-requeue
#SBATCH --output=/data/hohndor/gspa-v3/source-context-full-v1/slurm/%j.log
set -euo pipefail
root=/data/hohndor/gspa-v3/source-context-full-v1
jobout=$root/jobs/$SLURM_JOB_ID
mkdir -p "$jobout"
export TMPDIR=/scratch/$SLURM_JOB_ID
# Preserve a failed sealing attempt as well as failures in the full runner.
finish() {
 code=$?
 trap - EXIT
 if [ ! -f "$jobout/exit-code.txt" ]; then printf '%s\n' "$code" > "$jobout/exit-code.txt"; fi
 exit "$code"
}
trap finish EXIT
cp "$0" "$jobout/executed-launch.sbatch"
python3 "$root/seal.py" "$root" > "$jobout/seal.log" 2>&1
bash "$root/source_context_comparison.sbatch" "$root" full \
 /data/hohndor/gspa-v3/ec-prefix-identity-v3/jobs/10402/gspa-3.0.0.jar \
 /data/hohndor/gspa-v3/multisource-comparison-probe-v4/jobs/10403/results/verification.json
