#!/bin/bash
#SBATCH --job-name=gspa-sulfur-control
#SBATCH --partition=debug
#SBATCH --cpus-per-task=2
#SBATCH --nodelist=node002
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/data/hohndor/gspa-v3/sulfonate-diagnostic-v1/slurm/%j.log
set -euo pipefail
package=/data/hohndor/gspa-v3/empty-input-v1
current=/data/hohndor/gspa-v3/sulfonate-diagnostic-v1
work=/scratch/$SLURM_JOB_ID/gspa-sulfonate
jobout=$current/jobs/$SLURM_JOB_ID
mkdir -p "$work/fs" "$jobout"
export TMPDIR=/scratch/$SLURM_JOB_ID/tmp
mkdir -p "$TMPDIR"
finish() {
 code=$?
 trap - EXIT
 printf '%s\n' "$code" > "$jobout/exit-code.txt"
 if [ -d "$work/fs/home/leechuck/tmp/gspa-sulfonate-biology-v1" ]; then
  rsync -a "$work/fs/home/leechuck/tmp/gspa-sulfonate-biology-v1/" "$jobout/results/" || true
 fi
 exit "$code"
}
trap finish EXIT
hostname > "$jobout/hostname.txt"
cp "$current/controls-v1.sbatch" "$jobout/"
sha256sum "$package/before.tar.zst" "$package/python-runtime.tar.zst" "$current/source-v1.tar" > "$jobout/packages.sha256"
zstd -dc "$package/before.tar.zst" | tar -xf - -C "$work/fs"
zstd -dc "$package/python-runtime.tar.zst" | tar -xf - -C "$work/fs"
tar -xf "$current/source-v1.tar" -C "$work/fs"
tar -xf "$current/controls-source-v1.tar" -C "$work/fs"
cp -a "$current/jobs/10190/results/plan" "$work/fs/home/leechuck/tmp/gspa-sulfonate-biology-v1/plan"
container=(singularity exec --no-home --no-mount tmp --pwd /home/leechuck/Public/software/gspa
 --bind /usr,/lib,/lib64,/bin,/etc
 --bind "$work/fs/home/leechuck:/home/leechuck"
 --env TMPDIR=/home/leechuck/tmp --env OPENBLAS_NUM_THREADS=1 --env OMP_NUM_THREADS=1
 --env PATH=/home/leechuck/tmp/gspa-v3-testenv/bin:/usr/bin:/bin
 /data/hohndor/gspa-v3/observed-gap-context-modes-v1-relocation/container-root)
base=/home/leechuck/tmp/gspa-sulfonate-biology-v1
"${container[@]}" python benchmark/v3/check_sulfonate_medium.py --plan "$base/plan" --model "$base/iML1515.xml" --out "$base/medium-controls.json" > "$jobout/controls.log" 2>&1
