#!/bin/bash
#SBATCH --job-name=gspa-planner-feedback
#SBATCH --partition=debug
#SBATCH --cpus-per-task=2
#SBATCH --nodelist=node002
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/data/hohndor/gspa-v3/planner-feedback-v1/slurm/%j.log
set -euo pipefail
package=/data/hohndor/gspa-v3/empty-input-v1
current=/data/hohndor/gspa-v3/planner-feedback-v1
work=/scratch/$SLURM_JOB_ID/gspa-planner-feedback
jobout=$current/jobs/$SLURM_JOB_ID
mkdir -p "$work/fs" "$jobout"
export TMPDIR=/scratch/$SLURM_JOB_ID/tmp
mkdir -p "$TMPDIR"
finish() { code=$?; trap - EXIT; if [ -d "$work/fs/home/leechuck/tmp/planner-diagnostic" ]; then rsync -a "$work/fs/home/leechuck/tmp/planner-diagnostic/" "$jobout/diagnostic/"; fi; printf '%s\n' "$code" > "$jobout/exit-code.txt"; exit "$code"; }
trap finish EXIT
hostname > "$jobout/hostname.txt"
cp "$current/validate-v4.sbatch" "$jobout/"
sha256sum "$package/before.tar.zst" "$package/python-runtime.tar.zst" "$current/source-v3.tar.gz" > "$jobout/packages.sha256"
zstd -dc "$package/before.tar.zst" | tar -xf - -C "$work/fs"
zstd -dc "$package/python-runtime.tar.zst" | tar -xf - -C "$work/fs"
tar -xf "$current/source-v3.tar.gz" -C "$work/fs"
mkdir -p "$work/fs/home/leechuck/tmp"
container=(singularity exec --no-home --no-mount tmp --pwd /home/leechuck/Public/software/gspa
 --bind "$HOME/.config/borg-llm/api-key:/run/gspa-llm-key:ro" --bind /usr,/lib,/lib64,/bin,/etc --bind "$work/fs/home/leechuck:/home/leechuck"
 --env GSPA_LLM_URL=http://10.72.186.139:8000/v1 --env PYTHONPATH=/home/leechuck/Public/software/gspa/scripts/tests --env TMPDIR=/home/leechuck/tmp --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1
 --env PATH=/home/leechuck/tmp/gspa-v3-testenv/bin:/usr/bin:/bin
 /data/hohndor/gspa-v3/observed-gap-context-modes-v1-relocation/container-root)

"${container[@]}" python benchmark/v3/check_planner_feedback.py --before /home/leechuck/tmp/planner-before.py --out /home/leechuck/tmp/planner-diagnostic --live-qwen --key-file /run/gspa-llm-key > "$jobout/diagnostic.log" 2>&1
