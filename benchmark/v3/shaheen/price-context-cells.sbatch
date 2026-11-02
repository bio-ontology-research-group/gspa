#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=gspa-price-context
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --no-requeue
#SBATCH --exclude=nid00022,nid00497,nid01033
#SBATCH --output=/scratch/project/k10018/gspa-v3/out/price-context-v1/slurm/%x-%j.log
#
# Row-5 reaction-context arms on the Price cells (plan tasks/v3-row1-row5-price-calibration-context-2026-11-01.md,
# Part B). REPO_ROOT must be passed explicitly (Slurm spools this script). Nothing runs on login nodes.
#
# MODE=inputs  extract the 22 Price genomes' gapsmith tables and draft models from the frozen model batch, build the
#              observation JSON and the cohort metadata with plausibleReactions ($OUT_ROOT/inputs).
# MODE=cells   run the cells in CELLS (space/comma separated rc1 cell indices) with ARMS (default: all five) on this
#              node, packed by memory class (small <= 65, medium <= 100, large > 100, as in every Price collection).
# MODE=evaluate score all completed cells with the unchanged Price metrics and write the gap report.
set -euo pipefail
if [ -n "${SLURM_JOB_ID:-}" ] && [ -z "${REPO_ROOT:-}" ]; then
  echo "REPO_ROOT must be passed explicitly inside Slurm" >&2; exit 2
fi
REPO_ROOT=${REPO_ROOT:?}
[ -f "$REPO_ROOT/benchmark/v3/price_calibration_context/run_context_cell.py" ] || { echo "REPO_ROOT has no code slice: $REPO_ROOT" >&2; exit 2; }

OUT_ROOT=${OUT_ROOT:-/scratch/project/k10018/gspa-v3/out/price-context-v1}
RC1_PLAN=${RC1_PLAN:-/scratch/project/k10018/gspa-v3/out/price-production-rc1/collection/plan.json}
LAYOUT=${LAYOUT:-}                     # default: plan.json "layout"
REFERENCE=${REFERENCE:-/data/hohndor/gspa-v3/price-production-v1/reference}
JAR=${JAR:-/scratch/project/k10018/gspa-v3/out/routes-learned-v1/engine/gspa-3.0.0.jar}
EXPECT_JAR_SHA256=${EXPECT_JAR_SHA256:-08f8b5a132f4f1e5bf458f45b7604b1aa3b94388c96f170682bd2c27f2559c5a}
MODEL=${MODEL:-/scratch/project/k10018/gspa-v3/out/routes-learned-v1/model_logistic.json}
MODELS_TAR=${MODELS_TAR:-/scratch/project/k10018/gspa-v3/benchmark-models-v1/batches/batch_000/outputs.tar}
GSREF=${GSREF:-/data/hohndor/gspa-v3/gapsmith-reference-bacteria-v1/reference}
PATHWAY_GROUPS=${PATHWAY_GROUPS:-/data/hohndor/gspa-v3/cohort-context-modes-v2-inputs/pathway-groups.tsv}
COHORT_DIR=${COHORT_DIR:-$OUT_ROOT/price-cohort}       # cohort_orthogroups.tsv + cohort_metadata.json (22 genomes)
MANIFEST=${MANIFEST:-$REPO_ROOT/benchmark/manifests/suite6_price149_manifest.tsv}
EC_AUDIT=${EC_AUDIT:-$REPO_ROOT/benchmark/v3/results/ec-identity-v1/audit.json}
PYTHON=${PYTHON:-/project/k10018/sw/envs/gspa-tools/bin/python3.12}
PYLIB=${PYLIB:-/data/hohndor/gspa-v3/price-production-v1/venv/lib/python3.12/site-packages}
FSROOT=${FSROOT:-/scratch/project/k10018/gspa-v3/shared/arm-fs-d04c9e11ee64db47}
IMAGE=${IMAGE:-/project/k10018/sw/images/debian-trixie}
WORK_ROOT=${WORK_ROOT:-/scratch/project/k10018/gspa-v3/work}
MEM_BUDGET_G=${MEM_BUDGET_G:-340}
HEAP_SMALL=${HEAP_SMALL:-24g} HEAP_MEDIUM=${HEAP_MEDIUM:-96g} HEAP_LARGE=${HEAP_LARGE:-168g}
ALLOC_SMALL_G=${ALLOC_SMALL_G:-36} ALLOC_MEDIUM_G=${ALLOC_MEDIUM_G:-112} ALLOC_LARGE_G=${ALLOC_LARGE_G:-192}
TIMEOUT_SMALL=${TIMEOUT_SMALL:-7200} TIMEOUT_MEDIUM=${TIMEOUT_MEDIUM:-14400} TIMEOUT_LARGE=${TIMEOUT_LARGE:-36000}
ARMS=${ARMS:-stage0 gaps universal pathway organism}
FEATURE_EXPORT_MAX_INDEX=${FEATURE_EXPORT_MAX_INDEX:--1}   # feature export for cells with index <= this (sparse predictors)
KEEP_OUTPUTS=${KEEP_OUTPUTS:-ledgers}
MODE=${MODE:-cells}
JOB=${SLURM_JOB_ID:-local$$}
logical_data=/data/hohndor/gspa-v3
jobout=$OUT_ROOT/jobs/$JOB
work=$WORK_ROOT/price-context-$JOB
unset JAVA_TOOL_OPTIONS JDK_JAVA_OPTIONS _JAVA_OPTIONS PYTHONPATH PYTHONHOME
module load singularity
launcher=$(readlink -f "$(command -v singularity)")
[ -f "$FSROOT/.complete" ] || { echo "shared package tree missing: $FSROOT" >&2; exit 2; }
mkdir -p "$jobout" "$work/tmp" "$OUT_ROOT/slurm"
date --iso-8601=seconds > "$jobout/started-at.txt"; hostname > "$jobout/hostname.txt"
cp "$0" "$jobout/executed.sbatch"
echo "$EXPECT_JAR_SHA256  $JAR" | sha256sum -c - > "$jobout/jar-check.log"

# Options (binds, env) before the image; every path the engine reads must be bound.
prefix=("$launcher" exec --no-home --no-mount tmp
  --bind /project/k10018/sw/envs/gspa-tools:/project/k10018/sw/envs/gspa-tools:ro
  --bind "/project/k10018/gspa-v3:$logical_data:ro"
  --bind "$OUT_ROOT:$OUT_ROOT" --bind "$work:$work"
  --bind "$(dirname "$RC1_PLAN"):$(dirname "$RC1_PLAN"):ro"
  --bind "$(dirname "$JAR"):$(dirname "$JAR"):ro" --bind "$(dirname "$MODEL"):$(dirname "$MODEL"):ro"
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

if [ "$MODE" = "inputs" ]; then
  [ ! -e "$OUT_ROOT/inputs" ] || { echo "refusing: $OUT_ROOT/inputs exists" >&2; exit 2; }
  mkdir -p "$work/models"
  genomes=$("${prefix[@]}" "$IMAGE" "$PYTHON" -c 'import json,sys; print(" ".join(sorted(json.load(open(sys.argv[1]))["genomes"])))' "$COHORT_DIR/cohort_metadata.json")
  members=()
  for g in $genomes; do members+=("$g/$g-all-Reactions.tbl" "$g/$g-all-Pathways.tbl" "$g/$g-draft.xml" "$g/input_sha256.txt"); done
  /usr/bin/time -v -o "$jobout/extract.resources.txt" tar -xf "$MODELS_TAR" -C "$work/models" "${members[@]}"
  ( cd "$work/models" && find . -type f | sort | xargs sha256sum ) > "$jobout/extracted.sha256"
  sha256sum "$MODELS_TAR" > "$jobout/models-tar.sha256"
  mkdir -p "$OUT_ROOT/inputs-models"
  cp -r "$work/models/." "$OUT_ROOT/inputs-models/"
  "${prefix[@]}" --env "TMPDIR=$work/tmp" "$IMAGE" "$PYTHON" \
    "$REPO_ROOT/benchmark/v3/price_calibration_context/build_context_inputs.py" build \
    --models-dir "$OUT_ROOT/inputs-models" --reaction-graph "$GSREF/seed_reactions.tsv" \
    --cohort-metadata "$COHORT_DIR/cohort_metadata.json" --out "$OUT_ROOT/inputs" > "$jobout/build.log" 2>&1
  tail -n 40 "$jobout/build.log"
  rm -rf "$work"; date --iso-8601=seconds > "$jobout/finished-at.txt"; exit 0
fi

if [ "$MODE" = "evaluate" ]; then
  EVAL_OUT=${EVAL_OUT:-$OUT_ROOT/evaluation-$JOB}
  "${prefix[@]}" --env "TMPDIR=$work/tmp" "$IMAGE" "$PYTHON" \
    "$REPO_ROOT/benchmark/v3/price_calibration_context/evaluate_context.py" \
    --cells "$OUT_ROOT/cells" --manifest "$MANIFEST" --ec-audit "$EC_AUDIT" \
    --rc1-retained "${RC1_RETAINED:-/scratch/project/k10018/gspa-v3/out/price-reference-rc1/jobs/15713513/evaluation/retained}" \
    ${ROW1_POLICY:+--row1-policy "$ROW1_POLICY"} ${ROW1_SPLIT:+--row1-split "$ROW1_SPLIT"} \
    --out "$EVAL_OUT" > "$jobout/evaluate.log" 2>&1
  mkdir -p "$EVAL_OUT/gaps"
  for f in "$OUT_ROOT"/cells/*/*.gaps.jsonl; do cp "$f" "$EVAL_OUT/gaps/"; done
  "${prefix[@]}" --env "TMPDIR=$work/tmp" "$IMAGE" "$PYTHON" \
    "$REPO_ROOT/benchmark/v3/price_calibration_context/build_context_inputs.py" gap-report \
    --gaps-dir "$EVAL_OUT/gaps" --manifest "$MANIFEST" --reaction-graph "$GSREF/seed_reactions.tsv" --reaction-aliases "$GSREF/seed_Enzyme_Class_Reactions_Aliases_unique.tsv" \
    --currency "$GSREF/diffusion_mets.tsv" --pathway-groups "$PATHWAY_GROUPS" \
    --cohort-metadata "$OUT_ROOT/inputs/cohort_metadata.json" --out "$EVAL_OUT/gap_report.json" > "$jobout/gap-report.log" 2>&1
  tail -n 30 "$jobout/evaluate.log" "$jobout/gap-report.log"
  rm -rf "$work"; date --iso-8601=seconds > "$jobout/finished-at.txt"; exit 0
fi

: "${CELLS:?cell indices required}"
read -r -a cells <<< "${CELLS//[,+]/ }"
class_of() { if [ "$1" -le 65 ]; then echo SMALL; elif [ "$1" -le 100 ]; then echo MEDIUM; else echo LARGE; fi; }
total=0
for i in "${cells[@]}"; do c=$(class_of "$i"); a=ALLOC_${c}_G; total=$((total + ${!a})); done
echo "job $JOB on $(hostname): cells=[${cells[*]}] alloc=${total}G budget=${MEM_BUDGET_G}G arms=[$ARMS]"
mkdir -p "$OUT_ROOT/cells"

run_cell() {
  local i=$1 c heap timeout cw log rc=0 fe=()
  c=$(class_of "$i"); heap=HEAP_$c; timeout=TIMEOUT_$c; heap=${!heap}; timeout=${!timeout}
  cw=$work/cell-$i; log=$jobout/cell-$i; mkdir -p "$cw/tmp" "$log"
  [ "$i" -le "$FEATURE_EXPORT_MAX_INDEX" ] && fe=(--feature-export)
  date --iso-8601=seconds > "$log/started-at.txt"
  # shellcheck disable=SC2086
  "${prefix[@]}" --env "TMPDIR=$cw/tmp" "$IMAGE" "$PYTHON" \
    "$REPO_ROOT/benchmark/v3/price_calibration_context/run_context_cell.py" \
    --plan "$RC1_PLAN" --index "$i" --arms $ARMS --layout "$LAYOUT" --reference "$REFERENCE" --jar "$JAR" \
    --inputs "$OUT_ROOT/inputs" --orthogroups "$COHORT_DIR/cohort_orthogroups.tsv" --pathway-groups "$PATHWAY_GROUPS" \
    --reaction-graph "$GSREF/seed_reactions.tsv" --reaction-aliases "$GSREF/seed_Enzyme_Class_Reactions_Aliases_unique.tsv" \
    --model "$MODEL" --layout-metadata "$COHORT_DIR/layout_metadata.json" --manifest "$MANIFEST" --heap "$heap" --timeout "$timeout" --scratch "$cw/tmp" \
    --out "$OUT_ROOT/cells" --keep-outputs "$KEEP_OUTPUTS" "${fe[@]}" > "$log/cell.log" 2>&1 || rc=$?
  printf '%s\n' "$rc" > "$log/exit-code.txt"; date --iso-8601=seconds > "$log/finished-at.txt"
  [ "$rc" -ne 0 ] || rm -rf "$cw"
  return 0
}

POLL=${POLL:-15}
declare -A running=(); used=0
reap() { local p; for p in "${!running[@]}"; do if ! kill -0 "$p" 2>/dev/null; then wait "$p" || true; used=$((used - running[$p])); unset "running[$p]"; fi; done; }
for i in "${cells[@]}"; do
  c=$(class_of "$i"); a=ALLOC_${c}_G; need=${!a}
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
