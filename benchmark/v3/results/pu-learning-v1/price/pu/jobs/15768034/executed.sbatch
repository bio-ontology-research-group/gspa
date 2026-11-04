#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=gspa-price-contract
#SBATCH --nodes=1
#SBATCH --time=06:00:00
#SBATCH --no-requeue
#SBATCH --exclude=nid00022,nid00497,nid01033
#SBATCH --output=/scratch/project/k10018/gspa-v3/out/pu-learning-v1/price/pu/slurm/%x-%j.log
#
# Contract-bearing DIAMOND claims on the 22 Price genomes (plan tasks/v3-price-contract-learned-2026-11-02.md,
# steps I and G). REPO_ROOT must be passed explicitly (Slurm spools this script). Nothing runs on login nodes.
#
# MODE=cells      run the genome cells in CELLS (indices of $OUT_ROOT/plan.json) with ARMS (default: all nine arms of
#                 price_contract_learned/run_cell.py), packed by memory on this node.
# MODE=py         run $PY_SCRIPT (relative to REPO_ROOT) with $PY_ARGS inside the container.
# MODE=anomalies  the production anomaly extractor (roles_closure/run_price_anomalies.run, unchanged) on all 22 genomes.
set -euo pipefail
if [ -n "${SLURM_JOB_ID:-}" ] && [ -z "${REPO_ROOT:-}" ]; then
  echo "REPO_ROOT must be passed explicitly inside Slurm" >&2; exit 2
fi
REPO_ROOT=${REPO_ROOT:?}
[ -f "$REPO_ROOT/benchmark/v3/price_contract_learned/run_cell.py" ] || { echo "REPO_ROOT has no code slice: $REPO_ROOT" >&2; exit 2; }

OUT_ROOT=${OUT_ROOT:-/scratch/project/k10018/gspa-v3/out/price-contract-learned-v1}
CTX_ROOT=${CTX_ROOT:-/scratch/project/k10018/gspa-v3/out/price-context-v1}
RC1_PLAN=${RC1_PLAN:-$OUT_ROOT/plan/plan.json}
LAYOUT=${LAYOUT:-}                     # default: plan.json "layout"
REFERENCE=${REFERENCE:-/data/hohndor/gspa-v3/price-production-v1/reference}
JAR=${JAR:-/scratch/project/k10018/gspa-v3/out/routes-learned-v1/engine/gspa-3.0.0.jar}
EXPECT_JAR_SHA256=${EXPECT_JAR_SHA256:-08f8b5a132f4f1e5bf458f45b7604b1aa3b94388c96f170682bd2c27f2559c5a}
MODEL=${MODEL:-$OUT_ROOT/models/model_logistic.json}
MODEL_GBDT=${MODEL_GBDT:-$OUT_ROOT/models/model_gbdt.json}
EXPECT_MODEL_SHA256=${EXPECT_MODEL_SHA256:-c9df93d60e7b290b7562fb8f0b25f69f72ecea5edcb9b17492251dfd1cd5512c}
GSREF=${GSREF:-/data/hohndor/gspa-v3/gapsmith-reference-bacteria-v1/reference}
PATHWAY_GROUPS=${PATHWAY_GROUPS:-/data/hohndor/gspa-v3/cohort-context-modes-v2-inputs/pathway-groups.tsv}
COHORT_DIR=${COHORT_DIR:-$CTX_ROOT/price-cohort}       # cohort_orthogroups.tsv + cohort_metadata.json (22 genomes)
MANIFEST=${MANIFEST:-$REPO_ROOT/benchmark/manifests/suite6_price149_manifest.tsv}
INPUTS=${INPUTS:-$CTX_ROOT/inputs}
PYTHON=${PYTHON:-/project/k10018/sw/envs/gspa-tools/bin/python3.12}
PYLIB=${PYLIB:-/data/hohndor/gspa-v3/price-production-v1/venv/lib/python3.12/site-packages}
FSROOT=${FSROOT:-/scratch/project/k10018/gspa-v3/shared/arm-fs-d04c9e11ee64db47}
IMAGE=${IMAGE:-/project/k10018/sw/images/debian-trixie}
WORK_ROOT=${WORK_ROOT:-/scratch/project/k10018/gspa-v3/work}
MEM_BUDGET_G=${MEM_BUDGET_G:-340}
HEAP=${HEAP:-16g} ALLOC_G=${ALLOC_G:-20} TIMEOUT=${TIMEOUT:-7200}
ARMS=${ARMS:-stage0 stage1 learned learned_gbdt gaps gaps_learned universal pathway organism obs_learned rlc_learned cohort_universal cohort_pathway cohort_organism}
MODE=${MODE:-cells}
JOB=${SLURM_JOB_ID:-local$$}
logical_data=/data/hohndor/gspa-v3
jobout=$OUT_ROOT/jobs/$JOB
work=$WORK_ROOT/price-contract-$JOB
unset JAVA_TOOL_OPTIONS JDK_JAVA_OPTIONS _JAVA_OPTIONS PYTHONPATH PYTHONHOME
module load singularity
launcher=$(readlink -f "$(command -v singularity)")
[ -f "$FSROOT/.complete" ] || { echo "shared package tree missing: $FSROOT" >&2; exit 2; }
mkdir -p "$jobout" "$work/tmp" "$OUT_ROOT/slurm"
date --iso-8601=seconds > "$jobout/started-at.txt"; hostname > "$jobout/hostname.txt"
cp "$0" "$jobout/executed.sbatch"
echo "$EXPECT_JAR_SHA256  $JAR" | sha256sum -c - > "$jobout/jar-check.log"
echo "$EXPECT_MODEL_SHA256  $MODEL" | sha256sum -c - >> "$jobout/jar-check.log"
sha256sum "$MODEL_GBDT" "$RC1_PLAN" >> "$jobout/jar-check.log"

# Options (binds, env) before the image; every path the engine reads must be bound.
prefix=("$launcher" exec --no-home --no-mount tmp
  --bind /project/k10018/sw/envs/gspa-tools:/project/k10018/sw/envs/gspa-tools:ro
  --bind "/project/k10018/gspa-v3:$logical_data:ro"
  --bind "$OUT_ROOT:$OUT_ROOT" --bind "$work:$work"
  --bind "$(dirname "$RC1_PLAN"):$(dirname "$RC1_PLAN"):ro"
  --bind "$(dirname "$JAR"):$(dirname "$JAR"):ro"
  --bind "$CTX_ROOT/inputs:$CTX_ROOT/inputs:ro" --bind "$COHORT_DIR:$COHORT_DIR:ro"
  --bind "$REPO_ROOT:$REPO_ROOT:ro"
  --bind "$FSROOT/usr/lib/jvm/java-21-openjdk-amd64:/usr/lib/jvm/java-21-openjdk-amd64:ro"
  --bind "$FSROOT/home/leechuck:/home/leechuck:ro"
  --bind /usr/bin/time:/usr/bin/time:ro
  --env PATH=/usr/lib/jvm/java-21-openjdk-amd64/bin:/usr/bin:/bin
  --env LC_ALL=C.UTF-8 --env PYTHONDONTWRITEBYTECODE=1
  --env OPENBLAS_NUM_THREADS=1 --env OMP_NUM_THREADS=1 --env MKL_NUM_THREADS=1
  --env "PYTHONPATH=$REPO_ROOT/benchmark/v3:$PYLIB")
printf '%s\n' "${prefix[@]}" > "$jobout/container-prefix.txt"

if [ -z "$LAYOUT" ]; then
  LAYOUT=$("${prefix[@]}" "$IMAGE" "$PYTHON" -c 'import json,sys; print(json.load(open(sys.argv[1]))["layout"])' "$RC1_PLAN")
fi
echo "layout=$LAYOUT" >> "$jobout/started-at.txt"

if [ "$MODE" = "py" ]; then
  # One Python program of the code slice inside the container (e.g. evaluate_stage3_agent.py prepare), PY_ARGS verbatim.
  : "${PY_SCRIPT:?}"
  # shellcheck disable=SC2086
  "${prefix[@]}" --env "TMPDIR=$work/tmp" "$IMAGE" "$PYTHON" "$REPO_ROOT/$PY_SCRIPT" $PY_ARGS > "$jobout/py.log" 2>&1 || { tail -n 40 "$jobout/py.log"; exit 1; }
  tail -n 30 "$jobout/py.log"; date --iso-8601=seconds > "$jobout/finished-at.txt"; exit 0
fi

if [ "$MODE" = "anomalies" ]; then
  # Production anomaly extractor, unchanged code (run_price_anomalies.run): arms stage1, learned_logistic, learned_gbdt.
  "${prefix[@]}" --env "TMPDIR=$work/tmp" --env "JAVA_TOOL_OPTIONS=-Djava.io.tmpdir=$work/tmp" "$IMAGE" "$PYTHON" \
    "$REPO_ROOT/benchmark/v3/price_contract_learned/run_anomalies.py" \
    --out "$OUT_ROOT/anomalies" --parts "$OUT_ROOT/parts" --layouts "$LAYOUT" \
    --cohort-metadata "$COHORT_DIR/cohort_metadata.json" --jar "$JAR" --reference "$REFERENCE" \
    --models "$OUT_ROOT/models" --batch-threads "${BATCH_THREADS:-8}" --heap "${ANOMALY_HEAP:-48g}" > "$jobout/anomalies.log" 2>&1 || { tail -n 40 "$jobout/anomalies.log"; exit 1; }
  tail -n 20 "$jobout/anomalies.log"
  rm -rf "$work"; date --iso-8601=seconds > "$jobout/finished-at.txt"; exit 0
fi

: "${CELLS:?cell indices required}"
read -r -a cells <<< "${CELLS//[,+]/ }"
total=$(( ${#cells[@]} * ALLOC_G ))
echo "job $JOB on $(hostname): cells=[${cells[*]}] alloc=${total}G budget=${MEM_BUDGET_G}G arms=[$ARMS]"
mkdir -p "$OUT_ROOT/cells"

run_cell() {
  local i=$1 cw log rc=0
  cw=$work/cell-$i; log=$jobout/cell-$i; mkdir -p "$cw/tmp" "$log"
  date --iso-8601=seconds > "$log/started-at.txt"
  # shellcheck disable=SC2086
  "${prefix[@]}" --env "TMPDIR=$cw/tmp" "$IMAGE" "$PYTHON" \
    "$REPO_ROOT/benchmark/v3/price_contract_learned/run_cell.py" \
    --plan "$RC1_PLAN" --index "$i" --arms $ARMS --layout "$LAYOUT" --reference "$REFERENCE" --jar "$JAR" \
    --inputs "$INPUTS" --orthogroups "$COHORT_DIR/cohort_orthogroups.tsv" --pathway-groups "$PATHWAY_GROUPS" \
    --reaction-graph "$GSREF/seed_reactions.tsv" --reaction-aliases "$GSREF/seed_Enzyme_Class_Reactions_Aliases_unique.tsv" \
    --model "$MODEL" --model-gbdt "$MODEL_GBDT" --layout-metadata "$COHORT_DIR/layout_metadata.json" --manifest "$MANIFEST" \
    --heap "$HEAP" --timeout "$TIMEOUT" --scratch "$cw/tmp" --out "$OUT_ROOT/cells" > "$log/cell.log" 2>&1 || rc=$?
  printf '%s\n' "$rc" > "$log/exit-code.txt"; date --iso-8601=seconds > "$log/finished-at.txt"
  [ "$rc" -ne 0 ] || rm -rf "$cw"
  return 0
}

POLL=${POLL:-15}
declare -A running=(); used=0
reap() { local p; for p in "${!running[@]}"; do if ! kill -0 "$p" 2>/dev/null; then wait "$p" || true; used=$((used - running[$p])); unset "running[$p]"; fi; done; }
for i in "${cells[@]}"; do
  need=$ALLOC_G
  while [ "${#running[@]}" -gt 0 ] && [ $((used + need)) -gt "$MEM_BUDGET_G" ]; do sleep "$POLL"; reap; done
  run_cell "$i" & running[$!]=$need; used=$((used + need))
  echo "$(date --iso-8601=seconds) started cell $i (${need}G, in use ${used}G)"
done
while [ "${#running[@]}" -gt 0 ]; do sleep "$POLL"; reap; done
fails=0
for i in "${cells[@]}"; do s=$(cat "$jobout/cell-$i/exit-code.txt" 2>/dev/null || echo missing); echo "cell $i: exit $s"; [ "$s" = 0 ] || fails=$((fails + 1)); done
date --iso-8601=seconds > "$jobout/finished-at.txt"; printf '%s\n' "$fails" > "$jobout/failures.txt"
[ "$fails" -eq 0 ] && rm -rf "$work"
echo "done, failures=$fails"; [ "$fails" -eq 0 ]
