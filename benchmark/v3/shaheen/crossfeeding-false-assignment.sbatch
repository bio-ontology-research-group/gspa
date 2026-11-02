#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=gspa-v3-xfeed
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --exclusive
#SBATCH --time=01:30:00
#SBATCH --no-requeue
#SBATCH --exclude=nid00022,nid00497,nid01033
#SBATCH --output=/scratch/project/k10018/gspa-v3/out/absence-crossfeeding-v1/slurm/%x-%j.log
#
# One false-enzyme-assignment arm of the coculture sensitivity test
# (plan: tasks/v3-row8-row9-absence-crossfeeding-2026-11-01.md, Part B).
# All 7,920 frozen coculture plans are rerun through the production community engine
# with the perturbed gene-reaction rules, sharded over the node's cores, then merged.
#
#   MODE=uniform|targeted RATE=0.05 SEED=1 sbatch benchmark/v3/shaheen/crossfeeding-false-assignment.sbatch
#
set -euo pipefail

REPO_ROOT=${REPO_ROOT:-/scratch/project/k10018/gspa-v3/absence-crossfeeding-v1/code}
INPUTS=${INPUTS:-/scratch/project/k10018/gspa-v3/absence-crossfeeding-v1/inputs}
OUT_ROOT=${OUT_ROOT:-/scratch/project/k10018/gspa-v3/out/absence-crossfeeding-v1}
PYTHON=${PYTHON:-/project/k10018/sw/envs/gspa-v3/bin/python3}
MODE=${MODE:?uniform|targeted}
RATE=${RATE:?rate}
SEED=${SEED:?seed}
SHARDS=${SHARDS:-96}

ARM="$MODE-r$RATE-s$SEED"
OUT="$OUT_ROOT/$ARM"
mkdir -p "$OUT" "$OUT_ROOT/slurm"
export TMPDIR="$OUT/tmp"; mkdir -p "$TMPDIR"
cd "$REPO_ROOT"
echo "host=$(hostname) arm=$ARM shards=$SHARDS start=$(date -Is)"
$PYTHON -c "import libsbml,numpy,scipy;print('env',libsbml.getLibSBMLDottedVersion(),numpy.__version__,scipy.__version__)"
for i in $(seq 0 $((SHARDS-1))); do
  $PYTHON benchmark/v3/absence_crossfeeding/false_assignments.py run \
    --plan-dir "$INPUTS/community" --out "$OUT/s$i" \
    --ec-model "$INPUTS/models/iML1515.xml" --bs-model "$INPUTS/models/iYO844.xml" \
    --rate "$RATE" --seed "$SEED" --mode "$MODE" \
    --shard-index "$i" --shard-count "$SHARDS" > "$OUT/s$i.log" 2>&1 &
done
wait
$PYTHON benchmark/v3/absence_crossfeeding/false_assignments.py merge \
  --plan-dir "$INPUTS/community" \
  --shards $(for i in $(seq 0 $((SHARDS-1))); do echo "$OUT/s$i"; done) \
  --out "$OUT/merged"
rm -rf "$TMPDIR"
echo "ARM DONE $ARM end=$(date -Is)"
