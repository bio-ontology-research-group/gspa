#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=gspa-wp5-evaluate
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --no-requeue
#SBATCH --exclude=nid00022
#SBATCH --output=/scratch/project/k10018/gspa-v3/out/wp5-recovery-v1/slurm/%j.log
# WP5 evaluation phases (evaluate_wp5_recovery.py from the frozen source.tar.gz) on the Shaheen bundle/runs.
# PHASES: "join calibrate" (freezes evaluation/policy.json; review its SHA-256 first), then
# PHASES="evaluate report" POLICY_SHA256=<reviewed hash>.
set -euo pipefail
: "${PHASES:?}"
module load singularity
fs=$(ls -d /scratch/project/k10018/gspa-v3/shared/arm-fs-*/ | head -1); fs=${fs%/}
input_root=/project/k10018/gspa-v3/wp5-recovery-v1
O=/scratch/project/k10018/gspa-v3/out/wp5-recovery-v1 L=/data/hohndor/gspa-v3/wp5-recovery-v1-out
t=/scratch/project/k10018/gspa-v3/work/$SLURM_JOB_ID; jobout=$O/jobs/$SLURM_JOB_ID
mkdir -p "$t/wp5-source" "$jobout"; trap 'rm -rf "$t"' EXIT
cp "$0" "$jobout/executed.sbatch"; echo "$PHASES ${POLICY_SHA256:-}" > "$jobout/phases.txt"
tar -xzf "$input_root/source.tar.gz" -C "$t/wp5-source"
( cd "$t/wp5-source" && sha256sum -c "$input_root/source.sha256" > "$jobout/source-check.log" )
src=/home/leechuck/tmp/arm-tmp/wp5-source/benchmark/v3
c=(singularity exec --no-home --no-mount tmp --pwd "$L" --bind "$O:$L"
  --bind "$fs/home/leechuck:/home/leechuck:ro" --bind "$t:/home/leechuck/tmp/arm-tmp"
  --env XDG_CACHE_HOME=/home/leechuck/tmp/arm-tmp/.cache --env PYTHONPATH=$src --env PYTHONDONTWRITEBYTECODE=1
  /project/k10018/sw/images/debian-trixie /home/leechuck/tmp/gspa-v3-testenv/bin/python)
for phase in $PHASES; do
  case $phase in
    join|calibrate) "${c[@]}" "$src/evaluate_wp5_recovery.py" "$phase" --bundle bundle --runs runs --out evaluation ;;
    evaluate) "${c[@]}" "$src/evaluate_wp5_recovery.py" evaluate --bundle bundle --runs runs --out evaluation --policy-sha256 "${POLICY_SHA256:?}" ;;
    report) "${c[@]}" "$src/evaluate_wp5_recovery.py" report --out evaluation ;;
    *) echo "unknown phase $phase" >&2; exit 2 ;;
  esac > "$jobout/$phase.log" 2>&1 || { echo "phase $phase failed" >&2; tail -20 "$jobout/$phase.log" >&2; exit 1; }
  echo "phase $phase ok"
done
[ -f "$O/evaluation/policy.json" ] && sha256sum "$O/evaluation/policy.json"
