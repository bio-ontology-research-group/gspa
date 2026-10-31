#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=gspa-followup-score
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=12:00:00
#SBATCH --no-requeue
#SBATCH --exclude=nid00022,nid00497,nid01033
#SBATCH --output=/scratch/project/k10018/gspa-v3/out/stage2b-followups-v1/slurm/score-%j.log
#
# Price scoring of sparse update records (experiments 2, 4, 5) against stage 0 and stage 2B with the frozen metric
# functions, one process per predictor, inside the same container as the Stage 2B evaluation.
set -euo pipefail
: "${REPO_ROOT:?REPO_ROOT must point at the code checkout}"
: "${UPDATES:?UPDATES directory with <index>.<arm>.json}"
: "${ARMS:?ARMS (space separated arm labels)}"
: "${LABEL:?LABEL (output folder name)}"
OUT_ROOT=${OUT_ROOT:-/scratch/project/k10018/gspa-v3/out/stage2b-followups-v1}
COLLECTION=${COLLECTION:-/scratch/project/k10018/gspa-v3/out/price-stage2b-rc1/collection}
AUDIT=${AUDIT:-/scratch/project/k10018/gspa-v3/out/price-reference-rc1/ec-audit.json}
MANIFEST=${MANIFEST:-$REPO_ROOT/benchmark/manifests/suite6_price149_manifest.tsv}
PYLIB=${PYLIB:-/data/hohndor/gspa-v3/price-production-v1/venv/lib/python3.12/site-packages}
IMAGE=${IMAGE:-/project/k10018/sw/images/debian-trixie}
OUT=$OUT_ROOT/score-$LABEL-${SLURM_JOB_ID:-local}
module load singularity
mkdir -p "$OUT" "$OUT_ROOT/slurm"
[ "$(sha256sum "$AUDIT" | cut -d' ' -f1)" = 5302494b16ecd3f76fd96d24e48e27c298206b27608be72c5c66f1d197f4f79f ]
sha256sum "$AUDIT" "$MANIFEST" "$REPO_ROOT/tasks/v3-stage2b-followups-design-2026-10-31.md" "$REPO_ROOT/benchmark/v3/followups/price_score.py" > "$OUT/inputs.sha256"
for P in ${PREDICTORS:-topec mapred enzbert dpz graphec clean}; do
  singularity exec --no-home --no-mount tmp \
    --bind /project/k10018/sw/envs/gspa-tools:/project/k10018/sw/envs/gspa-tools:ro --bind /project/k10018/sw/envs/gspa-v3:/project/k10018/sw/envs/gspa-v3:ro \
    --bind /project/k10018/gspa-v3:/data/hohndor/gspa-v3:ro --bind /project/k10018/code:/project/k10018/code:ro \
    --bind /scratch/project/k10018/gspa-v3/out:/scratch/project/k10018/gspa-v3/out \
    --bind /scratch/project/k10018/gspa-v3/shared/arm-fs-d04c9e11ee64db47/usr/lib/jvm/java-21-openjdk-amd64:/usr/lib/jvm/java-21-openjdk-amd64:ro \
    --bind /scratch/project/k10018/gspa-v3/shared/arm-fs-d04c9e11ee64db47/home/leechuck:/home/leechuck:ro --bind /usr/bin/time:/usr/bin/time:ro \
    --env PYTHONDONTWRITEBYTECODE=1 --env "PYTHONPATH=$REPO_ROOT/benchmark/v3:$PYLIB" \
    "$IMAGE" /project/k10018/sw/envs/gspa-tools/bin/python3.12 -m followups.price_score \
    --collection "$COLLECTION" --manifest "$MANIFEST" --ec-audit "$AUDIT" --updates-dir "$UPDATES" --predictors "$P" --arms $ARMS \
    --out "$OUT/report-$P" > "$OUT/score-$P.log" 2>&1 &
done
fail=0; for j in $(jobs -p); do wait "$j" || fail=$((fail+1)); done
echo "failed predictors: $fail"; echo "exit $fail" > "$OUT/exit-code.txt"; [ "$fail" -eq 0 ]
