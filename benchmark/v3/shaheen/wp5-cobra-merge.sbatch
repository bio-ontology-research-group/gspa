#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=gspa-wp5-cobra-merge
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#SBATCH --exclude=nid00022
#SBATCH --output=/scratch/project/k10018/gspa-v3/out/wp5-recovery-v1/slurm/%j.log
# Merge wp5-cobra-shards.sbatch reports (SHARDS, default 256) with merge_wp5_cobra.py from the code slice.
set -euo pipefail
module load singularity
fs=$(ls -d /scratch/project/k10018/gspa-v3/shared/arm-fs-*/ | head -1); fs=${fs%/}
O=/scratch/project/k10018/gspa-v3/out/wp5-recovery-v1 L=/data/hohndor/gspa-v3/wp5-recovery-v1-out
t=/scratch/project/k10018/gspa-v3/work/$SLURM_JOB_ID; mkdir -p "$t"; trap 'rm -rf "$t"' EXIT
sha256sum /project/k10018/code/gspa-v3/benchmark/v3/merge_wp5_cobra.py
singularity exec --no-home --no-mount tmp --bind /project/k10018/code/gspa-v3:/opt/gspa:ro --bind "$O:$L" \
  --bind "$fs/home/leechuck:/home/leechuck:ro" --bind "$t:/home/leechuck/tmp/arm-tmp" \
  --env XDG_CACHE_HOME=/home/leechuck/tmp/arm-tmp/.cache --env PYTHONPATH=/opt/gspa/benchmark/v3 \
  /project/k10018/sw/images/debian-trixie /home/leechuck/tmp/gspa-v3-testenv/bin/python /opt/gspa/benchmark/v3/merge_wp5_cobra.py \
  --runs "$L/runs" --shard-dir "$L/oracle/cobra-${ARM_STAGE:-deterministic}.shards" --shards "${SHARDS:-256}" --arm-stages "${ARM_STAGE:-deterministic}" \
  --out "$L/oracle/cobra-${ARM_STAGE:-deterministic}.merged.json"
