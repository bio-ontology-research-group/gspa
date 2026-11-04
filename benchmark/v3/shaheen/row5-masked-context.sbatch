#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=gspa-row5-masked
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --no-requeue
#SBATCH --exclude=nid00022,nid00497,nid01033
#SBATCH --output=/scratch/project/k10018/gspa-v3/out/row5-masked-context-v1/slurm/%x-%j.log
#
# Row-5 masked-locus reaction-context benchmark (plan tasks/v3-row5-context-improvement-2026-11-02.md).
# REPO_ROOT must be passed explicitly (Slurm spools this script). Nothing runs on login nodes.
#
# MODE=stage     extract the 40 genomes' gapsmith tables and draft models into $OUT_ROOT/inputs-models,
#                run registry.py loci and registry.py inputs, copy Price DIAMOND claims with sha256 check
#                against claims_report.json, run filter-claims for CLEAN/DeepProZyme/GraphEC (22 cells each)
#                into $OUT_ROOT/claims-filtered/<predictor>/<genome>/claims.jsonl, then claim-sets.
# MODE=cells     run the <locus_id>:<claim_set> tokens in CELLS with ARMS on this node, packed by memory class
#                from published claim count (<= 200k -> 24g/36G, <= 2M -> 96g/112G, else 168g/192G).
# MODE=evaluate  run evaluate.py on SCOPE (development|heldout|all) with optional CHOOSE_WINDOW=1.
set -euo pipefail
if [ -n "${SLURM_JOB_ID:-}" ] && [ -z "${REPO_ROOT:-}" ]; then
  echo "REPO_ROOT must be passed explicitly inside Slurm" >&2; exit 2
fi
REPO_ROOT=${REPO_ROOT:?}
[ -f "$REPO_ROOT/benchmark/v3/row5_masked_context/run_cell.py" ] || { echo "REPO_ROOT has no code slice: $REPO_ROOT" >&2; exit 2; }

OUT_ROOT=${OUT_ROOT:-/scratch/project/k10018/gspa-v3/out/row5-masked-context-v1}
RC1_PLAN=${RC1_PLAN:-/scratch/project/k10018/gspa-v3/out/price-production-rc1/collection/plan.json}
PRICE_LAYOUT=${PRICE_LAYOUT:-}
COHORT_LAYOUT_DIR=${COHORT_LAYOUT_DIR:-/data/hohndor/gspa-v3/training-cohort/layout}
COHORT_BATCH_MANIFEST=${COHORT_BATCH_MANIFEST:-/data/hohndor/gspa-v3/training-cohort/batch_manifest.tsv}
COHORT_LABELS=${COHORT_LABELS:-/data/hohndor/gspa-v3/training-cohort/labels/labels.tsv}
COHORT_SPLITS=${COHORT_SPLITS:-$REPO_ROOT/benchmark/v3/results/training-cohort/splits.tsv}
PRICE_SPLIT=${PRICE_SPLIT:-$REPO_ROOT/benchmark/v3/results/price-calibration-context-v1/row1/split.json}
PRICE148_MANIFEST=${PRICE148_MANIFEST:-$REPO_ROOT/benchmark/manifests/suite6_price149_manifest.tsv}
PRICE_SWISSPROT=${PRICE_SWISSPROT:-$REPO_ROOT/benchmark/v3/offgrowth/data/swissprot_experimental_ecs_price22.tsv}
COHORT_ORGANISM_METADATA=${COHORT_ORGANISM_METADATA:-/data/hohndor/gspa-v3/cohort-context-modes-v2-inputs/organism-metadata.json}
COHORT_METADATA=${COHORT_METADATA:-/data/hohndor/gspa-v3/training-cohort/cohort_inputs/cohort_metadata.json}
COHORT_ORTHOGROUPS=${COHORT_ORTHOGROUPS:-/data/hohndor/gspa-v3/training-cohort/cohort_inputs/cohort_orthogroups.tsv}
PRICE_COHORT_DIR=${PRICE_COHORT_DIR:-/scratch/project/k10018/gspa-v3/out/price-context-v1/price-cohort}
PRICE_DIAMOND_DIR=${PRICE_DIAMOND_DIR:-/scratch/project/k10018/gspa-v3/out/price-contract-learned-v1/claims}
COHORT_DIAMOND_DIR=${COHORT_DIAMOND_DIR:-$OUT_ROOT/claims-diamond-cohort}
MODELS_TAR=${MODELS_TAR:-/scratch/project/k10018/gspa-v3/benchmark-models-v1/batches/batch_000/outputs.tar}
COHORT_MODELS_TAR=${COHORT_MODELS_TAR:-}
REFERENCE=${REFERENCE:-/data/hohndor/gspa-v3/price-production-v1/reference}
JAR=${JAR:-/scratch/project/k10018/gspa-v3/out/routes-learned-v1/engine/gspa-3.0.0.jar}
EXPECT_JAR_SHA256=${EXPECT_JAR_SHA256:-08f8b5a132f4f1e5bf458f45b7604b1aa3b94388c96f170682bd2c27f2559c5a}
JAR_CTX=${JAR_CTX:-$OUT_ROOT/engine/gspa-3.1.0-dev.jar}
EXPECT_JAR_CTX_SHA256=${EXPECT_JAR_CTX_SHA256:-bdc44e5f418ce6fbf6c637137f11c167f87dc973d05861ee7f7f758ec73a6122}   # 3.1.0-dev with the EC-namespace repair (first build c51d27d9… ran the pilot only)
MODEL=${MODEL:-/scratch/project/k10018/gspa-v3/out/routes-learned-v1/model_logistic.json}
GSREF=${GSREF:-/data/hohndor/gspa-v3/gapsmith-reference-bacteria-v1/reference}
PATHWAY_GROUPS=${PATHWAY_GROUPS:-/data/hohndor/gspa-v3/cohort-context-modes-v2-inputs/pathway-groups.tsv}
PYTHON=${PYTHON:-/project/k10018/sw/envs/gspa-tools/bin/python3.12}
PYLIB=${PYLIB:-/data/hohndor/gspa-v3/price-production-v1/venv/lib/python3.12/site-packages}
FSROOT=${FSROOT:-/scratch/project/k10018/gspa-v3/shared/arm-fs-d04c9e11ee64db47}
IMAGE=${IMAGE:-/project/k10018/sw/images/debian-trixie}
WORK_ROOT=${WORK_ROOT:-/scratch/project/k10018/gspa-v3/work}
MEM_BUDGET_G=${MEM_BUDGET_G:-340}
HEAP_SMALL=${HEAP_SMALL:-24g} HEAP_MEDIUM=${HEAP_MEDIUM:-96g} HEAP_LARGE=${HEAP_LARGE:-168g}
ALLOC_SMALL_G=${ALLOC_SMALL_G:-36} ALLOC_MEDIUM_G=${ALLOC_MEDIUM_G:-112} ALLOC_LARGE_G=${ALLOC_LARGE_G:-192}
TIMEOUT_SMALL=${TIMEOUT_SMALL:-7200} TIMEOUT_MEDIUM=${TIMEOUT_MEDIUM:-14400} TIMEOUT_LARGE=${TIMEOUT_LARGE:-36000}
WINDOW_BP=${WINDOW_BP:-5000 10000 20000}
ARMS=${ARMS:-}
KEEP_OUTPUTS=${KEEP_OUTPUTS:-ledgers}
MODE=${MODE:-cells}
CELLS_DIR=${CELLS_DIR:-$OUT_ROOT/cells}      # cell outputs (the pilot wrote to $OUT_ROOT/cells; the main run may use another folder)
JOB=${SLURM_JOB_ID:-local$$}
logical_data=/data/hohndor/gspa-v3
jobout=$OUT_ROOT/jobs/$JOB
work=$WORK_ROOT/row5-masked-$JOB
unset JAVA_TOOL_OPTIONS JDK_JAVA_OPTIONS _JAVA_OPTIONS PYTHONPATH PYTHONHOME
module load singularity
launcher=$(readlink -f "$(command -v singularity)")
[ -f "$FSROOT/.complete" ] || { echo "shared package tree missing: $FSROOT" >&2; exit 2; }
mkdir -p "$jobout" "$work/tmp" "$OUT_ROOT/slurm"
date --iso-8601=seconds > "$jobout/started-at.txt"; hostname > "$jobout/hostname.txt"
cp "$0" "$jobout/executed.sbatch"
echo "$EXPECT_JAR_SHA256  $JAR" | sha256sum -c - > "$jobout/jar-check.log"
if [ -f "$JAR_CTX" ] && [ -n "$EXPECT_JAR_CTX_SHA256" ]; then
  echo "$EXPECT_JAR_CTX_SHA256  $JAR_CTX" | sha256sum -c - >> "$jobout/jar-check.log"
fi

# Options (binds, env) before the image; every path the engine reads must be bound.
prefix=("$launcher" exec --no-home --no-mount tmp
  --bind /project/k10018/sw/envs/gspa-tools:/project/k10018/sw/envs/gspa-tools:ro
  --bind "/project/k10018/gspa-v3:$logical_data:ro"
  --bind "$OUT_ROOT:$OUT_ROOT" --bind "$work:$work"
  --bind "$(dirname "$RC1_PLAN"):$(dirname "$RC1_PLAN"):ro"
  --bind "$(dirname "$JAR"):$(dirname "$JAR"):ro"
  --bind "$(dirname "$MODEL"):$(dirname "$MODEL"):ro"
  --bind "$PRICE_COHORT_DIR:$PRICE_COHORT_DIR:ro"
  --bind "$(dirname "$PRICE_DIAMOND_DIR"):$(dirname "$PRICE_DIAMOND_DIR"):ro"
  --bind "$REPO_ROOT:$REPO_ROOT:ro"
  --bind "$FSROOT/usr/lib/jvm/java-21-openjdk-amd64:/usr/lib/jvm/java-21-openjdk-amd64:ro"
  --bind "$FSROOT/home/leechuck:/home/leechuck:ro"
  --bind /usr/bin/time:/usr/bin/time:ro
  --env PATH=/usr/lib/jvm/java-21-openjdk-amd64/bin:/usr/bin:/bin
  --env LC_ALL=C.UTF-8 --env PYTHONDONTWRITEBYTECODE=1
  --env OPENBLAS_NUM_THREADS=1 --env OMP_NUM_THREADS=1 --env MKL_NUM_THREADS=1
  --env "PYTHONPATH=$REPO_ROOT/benchmark/v3:$PYLIB")
printf '%s\n' "${prefix[@]}" > "$jobout/container-prefix.txt"

if [ -z "$PRICE_LAYOUT" ]; then
  PRICE_LAYOUT=$("${prefix[@]}" "$IMAGE" "$PYTHON" -c 'import json,sys; print(json.load(open(sys.argv[1]))["layout"])' "$RC1_PLAN")
fi
echo "price_layout=$PRICE_LAYOUT" >> "$jobout/started-at.txt"

if [ "$MODE" = "stage" ]; then
  mkdir -p "$work/models" "$OUT_ROOT/inputs-models" "$OUT_ROOT/registry" "$OUT_ROOT/inputs" "$OUT_ROOT/claims-diamond-price/claims" "$OUT_ROOT/claims-filtered"
  genomes=$("${prefix[@]}" "$IMAGE" "$PYTHON" -c '
import csv, json, sys
p_split = json.load(open(sys.argv[1]))["assignment"].keys()
c_split = [r["genome_id"].strip() for r in csv.DictReader(open(sys.argv[2]), delimiter="\t") if r.get("genome_id")]
print(" ".join(sorted(set(p_split) | set(c_split))))
' "$PRICE_SPLIT" "$COHORT_SPLITS")
  members=()
  for g in $genomes; do
    members+=("$g/$g-all-Reactions.tbl" "$g/$g-all-Pathways.tbl" "$g/$g-draft.xml" "$g/input_sha256.txt")
  done
  /usr/bin/time -v -o "$jobout/extract.resources.txt" tar --ignore-failed-read -xf "$MODELS_TAR" -C "$work/models" "${members[@]}" || true
  if [ -n "$COHORT_MODELS_TAR" ] && [ -f "$COHORT_MODELS_TAR" ]; then
    /usr/bin/time -v -o "$jobout/extract-cohort.resources.txt" tar --ignore-failed-read -xf "$COHORT_MODELS_TAR" -C "$work/models" "${members[@]}" || true
  fi
  cp -r "$work/models/." "$OUT_ROOT/inputs-models/"
  ( cd "$OUT_ROOT/inputs-models" && find . -type f | sort | xargs sha256sum ) > "$jobout/extracted.sha256"
  sha256sum "$MODELS_TAR" > "$jobout/models-tar.sha256"

  # 1. Copy Price DIAMOND claims with sha256 check against claims_report.json
  "${prefix[@]}" "$IMAGE" "$PYTHON" -c '
import hashlib, json, shutil, sys
from pathlib import Path
src_dir = Path(sys.argv[1])
dst_dir = Path(sys.argv[2])
rep_path = src_dir / "claims_report.json" if (src_dir / "claims_report.json").is_file() else src_dir.parent / "claims_report.json"
rep = json.loads(rep_path.read_text())
if rep_path.resolve() != (dst_dir / "claims_report.json").resolve():
    shutil.copy2(rep_path, dst_dir / "claims_report.json")
for g, info in sorted(rep.get("genomes", {}).items()):
    c_src = src_dir / "claims" / f"{g}.claims.jsonl"
    if not c_src.is_file():
        c_src = src_dir / f"{g}.claims.jsonl"
    h = hashlib.sha256(c_src.read_bytes()).hexdigest()
    exp = info.get("claims_sha256") or info.get("sha256")
    if exp and h != exp:
        raise SystemExit(f"SHA256 mismatch for {c_src}: got {h}, expected {exp}")
    c_dst = dst_dir / "claims" / f"{g}.claims.jsonl"
    c_dst.parent.mkdir(parents=True, exist_ok=True)
    if c_src.resolve() != c_dst.resolve():
        shutil.copy2(c_src, c_dst)
' "$PRICE_DIAMOND_DIR" "$OUT_ROOT/claims-diamond-price" > "$jobout/copy-price-diamond.log" 2>&1

  # 2. Build loci registry
  "${prefix[@]}" --env "TMPDIR=$work/tmp" "$IMAGE" "$PYTHON" \
    "$REPO_ROOT/benchmark/v3/row5_masked_context/registry.py" loci \
    --models-dir "$OUT_ROOT/inputs-models" \
    --price148-manifest "$PRICE148_MANIFEST" \
    --price-swissprot "$PRICE_SWISSPROT" \
    --cohort-labels "$COHORT_LABELS" \
    --cohort-batch-manifest "$COHORT_BATCH_MANIFEST" \
    --cohort-splits "$COHORT_SPLITS" \
    --price-split "$PRICE_SPLIT" \
    --out "$OUT_ROOT/registry" > "$jobout/loci.log" 2>&1

  # 3. Build masked/unmasked/sham inputs
  "${prefix[@]}" --env "TMPDIR=$work/tmp" "$IMAGE" "$PYTHON" \
    "$REPO_ROOT/benchmark/v3/row5_masked_context/registry.py" inputs \
    --loci "$OUT_ROOT/registry/loci.json" \
    --models-dir "$OUT_ROOT/inputs-models" \
    --reaction-graph "$GSREF/seed_reactions.tsv" \
    --price-cohort-metadata "$PRICE_COHORT_DIR/cohort_metadata.json" \
    --price-layout-metadata "$PRICE_COHORT_DIR/layout_metadata.json" \
    --price-layout "$PRICE_LAYOUT" \
    --price-diamond-dir "$OUT_ROOT/claims-diamond-price" \
    --cohort-organism-metadata "$COHORT_ORGANISM_METADATA" \
    --cohort-metadata "$COHORT_METADATA" \
    --cohort-layout-dir "$COHORT_LAYOUT_DIR" \
    --cohort-diamond-dir "$COHORT_DIAMOND_DIR" \
    --scratch "$work/tmp" \
    --out "$OUT_ROOT/inputs" > "$jobout/inputs.log" 2>&1

  # 4. Run filter-claims for the 22 cells each of CLEAN, DeepProZyme, GraphEC (in parallel on this node;
  #    each filter is a two-pass stream over one published matrix, memory O(proteins + ECs)).
  "${prefix[@]}" "$IMAGE" "$PYTHON" -c '
import json, sys
from pathlib import Path
from row5_masked_context.registry import DENSE_PREDICTORS, PREDICTOR_ALIASES
plan = json.loads(Path(sys.argv[1]).read_text())
for cell in plan.get("cells", []):
    pred = PREDICTOR_ALIASES.get(cell["predictor"], cell["predictor"].lower())
    if pred not in DENSE_PREDICTORS:
        continue
    print(str(Path(cell["part"]) / "claims.jsonl"), str(Path(sys.argv[2]) / pred / cell["assembly"] / "claims.jsonl"), sep="\t")
' "$RC1_PLAN" "$OUT_ROOT/claims-filtered" > "$jobout/filter-claims.todo"
  FILTER_PARALLEL=${FILTER_PARALLEL:-22}
  export -n prefix 2>/dev/null || true
  printf '%s\0' "${prefix[@]}" > "$jobout/prefix.nul"
  while IFS=$'\t' read -r src dst; do
    printf '%s\t%s\n' "$src" "$dst"
  done < "$jobout/filter-claims.todo" | xargs -P "$FILTER_PARALLEL" -L 1 bash -c '
    src=$0; dst=$1
    mapfile -d "" pfx < '"$jobout/prefix.nul"'
    mkdir -p "$(dirname "$dst")"
    "${pfx[@]}" --env "TMPDIR='"$work/tmp"'" "'"$IMAGE"'" "'"$PYTHON"'" "'"$REPO_ROOT"'/benchmark/v3/row5_masked_context/registry.py" filter-claims --in "$src" --out "$dst" --receipt "$dst.receipt.json" > "$dst.log" 2>&1 || echo "FAILED $src" >> "'"$jobout"'/filter-claims.failed"
  ' > "$jobout/filter-claims.log" 2>&1
  [ ! -s "$jobout/filter-claims.failed" ] || { cat "$jobout/filter-claims.failed"; exit 1; }

  # 5. Assemble claim_sets.json
  "${prefix[@]}" --env "TMPDIR=$work/tmp" "$IMAGE" "$PYTHON" \
    "$REPO_ROOT/benchmark/v3/row5_masked_context/registry.py" claim-sets \
    --rc1-plan "$RC1_PLAN" \
    --diamond-price-dir "$OUT_ROOT/claims-diamond-price" \
    --diamond-cohort-dir "$COHORT_DIAMOND_DIR" \
    --price-layout "$PRICE_LAYOUT" \
    --cohort-layout-dir "$COHORT_LAYOUT_DIR" \
    --cohort-batch-manifest "$COHORT_BATCH_MANIFEST" \
    --filtered-dir "$OUT_ROOT/claims-filtered" \
    --out "$OUT_ROOT/registry/claim_sets.json" > "$jobout/claim-sets.log" 2>&1

  tail -n 20 "$jobout/loci.log" "$jobout/inputs.log" "$jobout/claim-sets.log"
  rm -rf "$work"; date --iso-8601=seconds > "$jobout/finished-at.txt"; exit 0
fi

if [ "$MODE" = "evaluate" ]; then
  SCOPE=${SCOPE:-heldout}
  CHOOSE_WINDOW=${CHOOSE_WINDOW:-0}
  PARAMETER_CHOICE=${PARAMETER_CHOICE:-$OUT_ROOT/evaluation-development/parameter_choice.json}
  EVAL_OUT=${EVAL_OUT:-$OUT_ROOT/evaluation-$SCOPE}
  cw_flags=()
  if [ "$CHOOSE_WINDOW" = "1" ]; then
    cw_flags+=(--choose-window --parameter-choice "$PARAMETER_CHOICE")
  elif [ -f "$PARAMETER_CHOICE" ] || [ "$SCOPE" = "heldout" ]; then
    cw_flags+=(--parameter-choice "$PARAMETER_CHOICE")
  fi
  # shellcheck disable=SC2086
  "${prefix[@]}" --env "TMPDIR=$work/tmp" "$IMAGE" "$PYTHON" \
    "$REPO_ROOT/benchmark/v3/row5_masked_context/evaluate.py" \
    --cells "$CELLS_DIR" \
    --loci "$OUT_ROOT/registry/loci.json" \
    --scope "$SCOPE" \
    --window-bp $WINDOW_BP \
    "${cw_flags[@]}" \
    --out "$EVAL_OUT" > "$jobout/evaluate.log" 2>&1
  tail -n 40 "$jobout/evaluate.log"
  rm -rf "$work"; date --iso-8601=seconds > "$jobout/finished-at.txt"; exit 0
fi

: "${CELLS:?cell tokens (<locus_id>:<claim_set>) required}"
read -r -a cells <<< "${CELLS//[,+]/ }"

SMALL_MAX_CLAIMS=${SMALL_MAX_CLAIMS:-500000}   # measured: 24 GB heap holds the cohort pass up to ~450k claims (row-5 B4)
MEDIUM_MAX_CLAIMS=${MEDIUM_MAX_CLAIMS:-2000000}
class_of_claims() {
  local n=$1
  if [ "$n" -le "$SMALL_MAX_CLAIMS" ]; then
    echo SMALL
  elif [ "$n" -le "$MEDIUM_MAX_CLAIMS" ]; then
    echo MEDIUM
  else
    echo LARGE
  fi
}

declare -A cell_claims=()
declare -A cell_group=()
while IFS=$'\t' read -r tok n_claims grp; do
  [ -n "$tok" ] || continue
  cell_claims["$tok"]=$n_claims
  cell_group["$tok"]=$grp
done < <("${prefix[@]}" "$IMAGE" "$PYTHON" -c '
import json, sys
loci = json.load(open(sys.argv[1]))["loci"]
csets = json.load(open(sys.argv[2]))["by_key"]
for tok in sys.argv[3:]:
    lid, cs = tok.split(":", 1)
    loc = loci[lid]
    g = loc["genome"]
    entry = csets[f"{cs}:{g}"]
    print(tok, int(entry["published_claims"]), loc["group"], sep="\t")
' "$OUT_ROOT/registry/loci.json" "$OUT_ROOT/registry/claim_sets.json" "${cells[@]}")

total=0
for tok in "${cells[@]}"; do
  c=$(class_of_claims "${cell_claims[$tok]:-50000}")
  a=ALLOC_${c}_G
  total=$((total + ${!a}))
done
echo "job $JOB on $(hostname): cells=[${cells[*]}] alloc=${total}G budget=${MEM_BUDGET_G}G arms=[${ARMS:-default}]"
mkdir -p "$CELLS_DIR"

run_one_cell() {
  local tok=$1 lid cs n_claims grp c heap timeout cw log rc=0 orthogroups layout_meta_flag=() arms_flag=()
  lid=${tok%%:*}
  cs=${tok#*:}
  n_claims=${cell_claims[$tok]:-50000}
  grp=${cell_group[$tok]:-cohort}
  c=$(class_of_claims "$n_claims")
  heap=HEAP_$c; timeout=TIMEOUT_$c; heap=${!heap}; timeout=${!timeout}
  if [ "$grp" = "price" ]; then
    orthogroups="$PRICE_COHORT_DIR/cohort_orthogroups.tsv"
    [ -f "$PRICE_COHORT_DIR/layout_metadata.json" ] && layout_meta_flag=(--layout-metadata "$PRICE_COHORT_DIR/layout_metadata.json")
  else
    orthogroups="$COHORT_ORTHOGROUPS"
  fi
  # shellcheck disable=SC2206
  [ -n "$ARMS" ] && arms_flag=(--arms $ARMS)
  cw="$work/cell-${lid}-${cs}"
  log="$jobout/cell-${lid}-${cs}"
  mkdir -p "$cw/tmp" "$log"
  date --iso-8601=seconds > "$log/started-at.txt"
  # shellcheck disable=SC2086
  "${prefix[@]}" --env "TMPDIR=$cw/tmp" "$IMAGE" "$PYTHON" \
    "$REPO_ROOT/benchmark/v3/row5_masked_context/run_cell.py" \
    --loci "$OUT_ROOT/registry/loci.json" \
    --locus-id "$lid" \
    --claim-set "$cs" \
    --claim-sets "$OUT_ROOT/registry/claim_sets.json" \
    --inputs "$OUT_ROOT/inputs" \
    "${arms_flag[@]}" \
    --jar "$JAR" \
    --jar-ctx "$JAR_CTX" \
    --reference "$REFERENCE" \
    --reaction-graph "$GSREF/seed_reactions.tsv" \
    --reaction-aliases "$GSREF/seed_Enzyme_Class_Reactions_Aliases_unique.tsv" \
    --pathway-groups "$PATHWAY_GROUPS" \
    --orthogroups "$orthogroups" \
    "${layout_meta_flag[@]}" \
    --model "$MODEL" \
    --window-bp $WINDOW_BP \
    --heap "$heap" \
    --timeout "$timeout" \
    --scratch "$cw/tmp" \
    --out "$CELLS_DIR" \
    --keep-outputs "$KEEP_OUTPUTS" > "$log/cell.log" 2>&1 || rc=$?
  printf '%s\n' "$rc" > "$log/exit-code.txt"
  date --iso-8601=seconds > "$log/finished-at.txt"
  [ "$rc" -ne 0 ] || rm -rf "$cw"
  return 0
}

POLL=${POLL:-15}
declare -A running=(); used=0
reap() {
  local p
  for p in "${!running[@]}"; do
    if ! kill -0 "$p" 2>/dev/null; then
      wait "$p" || true
      used=$((used - running[$p]))
      unset "running[$p]"
    fi
  done
}

for tok in "${cells[@]}"; do
  c=$(class_of_claims "${cell_claims[$tok]:-50000}")
  a=ALLOC_${c}_G
  need=${!a}
  while [ "${#running[@]}" -gt 0 ] && [ $((used + need)) -gt "$MEM_BUDGET_G" ]; do
    sleep "$POLL"
    reap
  done
  run_one_cell "$tok" &
  running[$!]=$need
  used=$((used + need))
  echo "$(date --iso-8601=seconds) started cell $tok (${need}G, in use ${used}G)"
done
while [ "${#running[@]}" -gt 0 ]; do
  sleep "$POLL"
  reap
done

fails=0
for tok in "${cells[@]}"; do
  lid=${tok%%:*}
  cs=${tok#*:}
  s=$(cat "$jobout/cell-${lid}-${cs}/exit-code.txt" 2>/dev/null || echo missing)
  echo "cell $tok: exit $s"
  [ "$s" = 0 ] || fails=$((fails + 1))
done
date --iso-8601=seconds > "$jobout/finished-at.txt"
printf '%s\n' "$fails" > "$jobout/failures.txt"
[ "$fails" -eq 0 ] && rm -rf "$work"
echo "done, failures=$fails"
[ "$fails" -eq 0 ]
