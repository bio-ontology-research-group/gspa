#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=gspa-row11-ec-go
#SBATCH --nodes=1
#SBATCH --time=06:00:00
#SBATCH --no-requeue
#SBATCH --exclude=nid00022,nid00497,nid01033
#SBATCH --output=/scratch/project/k10018/gspa-v3/out/row11-ec-go-agent-v1/slurm/%x-%j.log
#
# EC-aware row-11 comparison on the 22 Price genomes (plan tasks/v3-row11-ec-go-agent-2026-11-02.md,
# steps S3 and S4). REPO_ROOT and EXPECT_JAR_SHA256 must be passed explicitly. Nothing runs on login nodes.
#
# MODE=merge      merge 6 predictors' claims (raw_score >= MIN_RAW) for SET=ec6|ec6c across ASSEMBLIES (or all 22),
#                 up to PAR=6 assemblies in parallel inside the container.
# MODE=operons    derive distance-rule operons for all 22 layouts into $OUT_ROOT/operons.
# MODE=arms       run stage1, learned_logistic, learned_gbdt via row11_ec_go/run_arms.py for SET=ec6|ec6c.
# MODE=fullcheck  stage1 on full vs. filtered clean claims for one assembly (heap 168g), then compare_full_vs_filtered.py.
# MODE=py         run $PY_SCRIPT (relative to REPO_ROOT) with $PY_ARGS inside the container.
set -euo pipefail
if [ -n "${SLURM_JOB_ID:-}" ] && [ -z "${REPO_ROOT:-}" ]; then
  echo "REPO_ROOT must be passed explicitly inside Slurm" >&2; exit 2
fi
REPO_ROOT=${REPO_ROOT:?}
[ -f "$REPO_ROOT/benchmark/v3/row11_ec_go/run_arms.py" ] || { echo "REPO_ROOT has no code slice: $REPO_ROOT" >&2; exit 2; }

if [ -z "${EXPECT_JAR_SHA256:-}" ]; then
  echo "EXPECT_JAR_SHA256 must be set explicitly (no default hash hard-coded)" >&2; exit 2
fi

OUT_ROOT=${OUT_ROOT:-/scratch/project/k10018/gspa-v3/out/row11-ec-go-agent-v1}
CONTRACT_ROOT=${CONTRACT_ROOT:-/scratch/project/k10018/gspa-v3/out/price-contract-learned-v1}
CTX_ROOT=${CTX_ROOT:-/scratch/project/k10018/gspa-v3/out/price-context-v1}
RC1_PLAN=${RC1_PLAN:-$CONTRACT_ROOT/plan/plan.json}
PARTS_ROOT=${PARTS_ROOT:-/data/hohndor/gspa-v3/price-full-matrix-v2/parts}
PREDICTORS=${PREDICTORS:-clean dpz enzbert graphec mapred topec}
MIN_RAW=${MIN_RAW:-0.05}
EXTRA_ROOT=${EXTRA_ROOT:-$CONTRACT_ROOT/parts/swissprot_diamond_baseline}
LARGE_FILES_TSV=${LARGE_FILES_TSV:-$REPO_ROOT/benchmark/v3/results/price-contract-learned-v1/LARGE_FILES.tsv}
LAYOUT=${LAYOUT:-/data/hohndor/gspa-v3/price-full-matrix-v2/layout}
OPERONS=${OPERONS:-$OUT_ROOT/operons}
REFERENCE=${REFERENCE:-/data/hohndor/gspa-v3/price-production-v1/reference}
JAR=${JAR:-$OUT_ROOT/engine/gspa-3.0.0.jar}
MODELS=${MODELS:-$CONTRACT_ROOT/models}
MODEL=${MODEL:-$MODELS/model_logistic.json}
MODEL_GBDT=${MODEL_GBDT:-$MODELS/model_gbdt.json}
EXPECT_MODEL_SHA256=${EXPECT_MODEL_SHA256:-c9df93d60e7b290b7562fb8f0b25f69f72ecea5edcb9b17492251dfd1cd5512c}
COHORT_DIR=${COHORT_DIR:-$CTX_ROOT/price-cohort}
COHORT_METADATA=${COHORT_METADATA:-$COHORT_DIR/cohort_metadata.json}
PYTHON=${PYTHON:-/project/k10018/sw/envs/gspa-tools/bin/python3.12}
PYLIB=${PYLIB:-/data/hohndor/gspa-v3/price-production-v1/venv/lib/python3.12/site-packages}
FSROOT=${FSROOT:-/scratch/project/k10018/gspa-v3/shared/arm-fs-d04c9e11ee64db47}
IMAGE=${IMAGE:-/project/k10018/sw/images/debian-trixie}
WORK_ROOT=${WORK_ROOT:-/scratch/project/k10018/gspa-v3/work}
MODE=${MODE:?MODE required (merge|operons|arms|fullcheck|py)}
JOB=${SLURM_JOB_ID:-local$$}
logical_data=/data/hohndor/gspa-v3
jobout=$OUT_ROOT/jobs/$JOB
work=$WORK_ROOT/row11-ec-go-$JOB
unset JAVA_TOOL_OPTIONS JDK_JAVA_OPTIONS _JAVA_OPTIONS PYTHONPATH PYTHONHOME
module load singularity
launcher=$(readlink -f "$(command -v singularity)")
[ -f "$FSROOT/.complete" ] || { echo "shared package tree missing: $FSROOT" >&2; exit 2; }
mkdir -p "$jobout" "$work/tmp" "$OUT_ROOT/slurm"
date --iso-8601=seconds > "$jobout/started-at.txt"; hostname > "$jobout/hostname.txt"
cp "$0" "$jobout/executed.sbatch"
echo "$EXPECT_JAR_SHA256  $JAR" | sha256sum -c - > "$jobout/jar-check.log"
echo "$EXPECT_MODEL_SHA256  $MODEL" | sha256sum -c - >> "$jobout/jar-check.log"
sha256sum "$MODEL_GBDT" >> "$jobout/jar-check.log"

# Options (binds, env) before the image; every path the engine reads must be bound.
prefix=("$launcher" exec --no-home --no-mount tmp
  --bind /project/k10018/sw/envs/gspa-tools:/project/k10018/sw/envs/gspa-tools:ro
  --bind "/project/k10018/gspa-v3:$logical_data:ro"
  --bind "$OUT_ROOT:$OUT_ROOT" --bind "$work:$work"
  --bind "$(dirname "$JAR"):$(dirname "$JAR"):ro"
  --bind "$EXTRA_ROOT:$EXTRA_ROOT:ro"
  --bind "$MODELS:$MODELS:ro"
  --bind "$COHORT_DIR:$COHORT_DIR:ro"
  --bind "$REPO_ROOT:$REPO_ROOT:ro"
  --bind "$FSROOT/usr/lib/jvm/java-21-openjdk-amd64:/usr/lib/jvm/java-21-openjdk-amd64:ro"
  --bind "$FSROOT/home/leechuck:/home/leechuck:ro"
  --bind /usr/bin/time:/usr/bin/time:ro
  --env PATH=/usr/lib/jvm/java-21-openjdk-amd64/bin:/usr/bin:/bin
  --env LC_ALL=C.UTF-8 --env PYTHONDONTWRITEBYTECODE=1
  --env OPENBLAS_NUM_THREADS=1 --env OMP_NUM_THREADS=1 --env MKL_NUM_THREADS=1
  --env "PYTHONPATH=$REPO_ROOT/benchmark/v3:$PYLIB")
if [ -f "$RC1_PLAN" ]; then
  prefix+=(--bind "$(dirname "$RC1_PLAN"):$(dirname "$RC1_PLAN"):ro")
fi
printf '%s\n' "${prefix[@]}" > "$jobout/container-prefix.txt"

if [ -z "$LAYOUT" ] && [ -f "$RC1_PLAN" ]; then
  LAYOUT=$("${prefix[@]}" "$IMAGE" "$PYTHON" -c 'import json,sys; print(json.load(open(sys.argv[1]))["layout"])' "$RC1_PLAN")
fi
echo "mode=$MODE layout=$LAYOUT" >> "$jobout/started-at.txt"

if [ "$MODE" = "py" ]; then
  : "${PY_SCRIPT:?}"
  # shellcheck disable=SC2086
  "${prefix[@]}" --env "TMPDIR=$work/tmp" "$IMAGE" "$PYTHON" "$REPO_ROOT/$PY_SCRIPT" $PY_ARGS > "$jobout/py.log" 2>&1 || { tail -n 40 "$jobout/py.log"; exit 1; }
  tail -n 30 "$jobout/py.log"
  rm -rf "$work"; date --iso-8601=seconds > "$jobout/finished-at.txt"; exit 0
fi

if [ "$MODE" = "operons" ]; then
  mkdir -p "$OPERONS"
  "${prefix[@]}" --env "TMPDIR=$work/tmp" "$IMAGE" "$PYTHON" \
    "$REPO_ROOT/benchmark/v3/row11_ec_go/derive_operons.py" \
    --layouts "$LAYOUT" --out "$OPERONS" > "$jobout/operons.log" 2>&1 || { tail -n 40 "$jobout/operons.log"; exit 1; }
  tail -n 30 "$jobout/operons.log"
  rm -rf "$work"; date --iso-8601=seconds > "$jobout/finished-at.txt"; exit 0
fi

if [ "$MODE" = "merge" ]; then
  SET=${SET:?SET required (ec6|ec6c)}
  case "$SET" in
    ec6|ec6c) ;;
    *) echo "SET must be ec6 or ec6c: $SET" >&2; exit 2 ;;
  esac
  PAR=${PAR:-6}
  POLL=${POLL:-5}
  if [ -n "${ASSEMBLIES:-}" ]; then
    read -r -a assemblies <<< "${ASSEMBLIES//[,+]/ }"
  else
    mapfile -t assemblies < <("${prefix[@]}" "$IMAGE" "$PYTHON" -c 'import json,sys; print("\n".join(sorted(json.load(open(sys.argv[1]))["assemblies"])))' "$LAYOUT/report.json")
  fi
  [ "${#assemblies[@]}" -gt 0 ] || { echo "no assemblies resolved for merge" >&2; exit 2; }
  echo "job $JOB on $(hostname): mode=merge set=$SET par=$PAR assemblies=[${assemblies[*]}]"
  mkdir -p "$OUT_ROOT/parts/$SET"

  run_merge_one() {
    local asm=$1 mw log rc=0 extra_args=() extra_sha
    mw=$work/merge-$SET-$asm; log=$jobout/merge-$SET-$asm; mkdir -p "$mw/tmp" "$log" "$OUT_ROOT/parts/$SET/$asm"
    date --iso-8601=seconds > "$log/started-at.txt"
    if [ "$SET" = "ec6c" ]; then
      extra_sha=$(awk -F'\t' -v p="claims/${asm}.claims.jsonl" '$1 == p {print $2}' "$LARGE_FILES_TSV")
      if [ -z "$extra_sha" ]; then
        echo "missing sha256 for claims/${asm}.claims.jsonl in $LARGE_FILES_TSV" > "$log/merge.log"
        printf '2\n' > "$log/exit-code.txt"; date --iso-8601=seconds > "$log/finished-at.txt"
        return 0
      fi
      extra_args=(--extra-claims "$EXTRA_ROOT/$asm/claims.jsonl" --extra-sha256 "$extra_sha" --extra-name swissprot_diamond_baseline)
    fi
    # shellcheck disable=SC2086
    "${prefix[@]}" --env "TMPDIR=$mw/tmp" "$IMAGE" "$PYTHON" \
      "$REPO_ROOT/benchmark/v3/row11_ec_go/merge_claims.py" \
      --parts-root "$PARTS_ROOT" --predictors $PREDICTORS --assembly "$asm" \
      --min-raw "$MIN_RAW" --out "$OUT_ROOT/parts/$SET/$asm" "${extra_args[@]}" > "$log/merge.log" 2>&1 || rc=$?
    printf '%s\n' "$rc" > "$log/exit-code.txt"; date --iso-8601=seconds > "$log/finished-at.txt"
    [ "$rc" -ne 0 ] || rm -rf "$mw"
    return 0
  }

  declare -A running=()
  reap() {
    local p
    for p in "${!running[@]}"; do
      if ! kill -0 "$p" 2>/dev/null; then
        wait "$p" || true
        unset "running[$p]"
      fi
    done
  }
  for asm in "${assemblies[@]}"; do
    while [ "${#running[@]}" -ge "$PAR" ]; do sleep "$POLL"; reap; done
    run_merge_one "$asm" &
    running[$!]=$asm
    echo "$(date --iso-8601=seconds) started merge $SET/$asm (active ${#running[@]}/$PAR)"
  done
  while [ "${#running[@]}" -gt 0 ]; do sleep "$POLL"; reap; done

  fails=0
  for asm in "${assemblies[@]}"; do
    s=$(cat "$jobout/merge-$SET-$asm/exit-code.txt" 2>/dev/null || echo missing)
    echo "merge $SET/$asm: exit $s"
    [ "$s" = 0 ] || fails=$((fails + 1))
  done
  date --iso-8601=seconds > "$jobout/finished-at.txt"; printf '%s\n' "$fails" > "$jobout/failures.txt"
  [ "$fails" -eq 0 ] && rm -rf "$work"
  echo "merge done, failures=$fails"; [ "$fails" -eq 0 ]
  exit 0
fi

if [ "$MODE" = "arms" ]; then
  SET=${SET:?SET required (ec6|ec6c)}
  case "$SET" in
    ec6|ec6c) ;;
    *) echo "SET must be ec6 or ec6c: $SET" >&2; exit 2 ;;
  esac
  ARMS_HEAP=${ARMS_HEAP:-${ANOMALY_HEAP:-48g}}
  BATCH_THREADS=${BATCH_THREADS:-8}
  ANOMALY_FLOOR=${ANOMALY_FLOOR:-0.20}
  ANOMALY_MAX_CANDIDATES=${ANOMALY_MAX_CANDIDATES:-10}
  ANOMALY_MAX_LOCUS=${ANOMALY_MAX_LOCUS:-150}
  ARMS_OUT=${ARMS_OUT:-$OUT_ROOT/arms/$SET}
  mkdir -p "$ARMS_OUT"
  "${prefix[@]}" --env "TMPDIR=$work/tmp" --env "JAVA_TOOL_OPTIONS=-Djava.io.tmpdir=$work/tmp" "$IMAGE" "$PYTHON" \
    "$REPO_ROOT/benchmark/v3/row11_ec_go/run_arms.py" \
    --set "$SET" --parts "$OUT_ROOT/parts" --out "$ARMS_OUT" \
    --layouts "$LAYOUT" --operons "$OPERONS" \
    --cohort-metadata "$COHORT_METADATA" --jar "$JAR" --reference "$REFERENCE" \
    --models "$MODELS" --batch-threads "$BATCH_THREADS" --heap "$ARMS_HEAP" \
    --anomaly-floor "$ANOMALY_FLOOR" --anomaly-max-candidates "$ANOMALY_MAX_CANDIDATES" \
    --anomaly-max-locus "$ANOMALY_MAX_LOCUS" \
    > "$jobout/arms-$SET.log" 2>&1 || { tail -n 40 "$jobout/arms-$SET.log"; exit 1; }
  tail -n 20 "$jobout/arms-$SET.log"
  rm -rf "$work"; date --iso-8601=seconds > "$jobout/finished-at.txt"; exit 0
fi

if [ "$MODE" = "fullcheck" ]; then
  FULLCHECK_ASSEMBLY=${FULLCHECK_ASSEMBLY:-GCF_000015245.1}
  FULLCHECK_PREDICTOR=${FULLCHECK_PREDICTOR:-clean}
  FULLCHECK_HEAP=${FULLCHECK_HEAP:-168g}
  BATCH_THREADS=${BATCH_THREADS:-8}
  ANOMALY_FLOOR=${ANOMALY_FLOOR:-0.20}
  ANOMALY_MAX_CANDIDATES=${ANOMALY_MAX_CANDIDATES:-10}
  ANOMALY_MAX_LOCUS=${ANOMALY_MAX_LOCUS:-150}
  FULLCHECK_OUT=${FULLCHECK_OUT:-$OUT_ROOT/fullcheck/$FULLCHECK_PREDICTOR/$FULLCHECK_ASSEMBLY}
  filt_part=$FULLCHECK_OUT/part_filtered/$FULLCHECK_PREDICTOR/$FULLCHECK_ASSEMBLY
  full_part=$PARTS_ROOT/$FULLCHECK_PREDICTOR/$FULLCHECK_ASSEMBLY
  mkdir -p "$FULLCHECK_OUT" "$filt_part" "$OPERONS"

  if [ ! -f "$OPERONS/$FULLCHECK_ASSEMBLY.operons.tsv" ]; then
    "${prefix[@]}" --env "TMPDIR=$work/tmp" "$IMAGE" "$PYTHON" \
      "$REPO_ROOT/benchmark/v3/row11_ec_go/derive_operons.py" \
      --layouts "$LAYOUT" --out "$OPERONS" > "$jobout/fullcheck-operons.log" 2>&1 || { tail -n 40 "$jobout/fullcheck-operons.log"; exit 1; }
  fi

  "${prefix[@]}" --env "TMPDIR=$work/tmp" "$IMAGE" "$PYTHON" \
    "$REPO_ROOT/benchmark/v3/row11_ec_go/merge_claims.py" \
    --parts-root "$PARTS_ROOT" --predictors "$FULLCHECK_PREDICTOR" --assembly "$FULLCHECK_ASSEMBLY" \
    --min-raw "$MIN_RAW" --out "$filt_part" > "$jobout/fullcheck-merge.log" 2>&1 || { tail -n 40 "$jobout/fullcheck-merge.log"; exit 1; }

  "${prefix[@]}" --env "TMPDIR=$work/tmp" --env "JAVA_TOOL_OPTIONS=-Djava.io.tmpdir=$work/tmp" "$IMAGE" "$PYTHON" \
    "$REPO_ROOT/benchmark/v3/row11_ec_go/run_arms.py" \
    --fullcheck-part "$full_part" --out "$FULLCHECK_OUT/full" \
    --layouts "$LAYOUT" --operons "$OPERONS" \
    --cohort-metadata "$COHORT_METADATA" --jar "$JAR" --reference "$REFERENCE" \
    --models "$MODELS" --batch-threads "$BATCH_THREADS" --heap "$FULLCHECK_HEAP" \
    --anomaly-floor "$ANOMALY_FLOOR" --anomaly-max-candidates "$ANOMALY_MAX_CANDIDATES" \
    --anomaly-max-locus "$ANOMALY_MAX_LOCUS" \
    > "$jobout/fullcheck-full.log" 2>&1 || { tail -n 40 "$jobout/fullcheck-full.log"; exit 1; }

  "${prefix[@]}" --env "TMPDIR=$work/tmp" --env "JAVA_TOOL_OPTIONS=-Djava.io.tmpdir=$work/tmp" "$IMAGE" "$PYTHON" \
    "$REPO_ROOT/benchmark/v3/row11_ec_go/run_arms.py" \
    --fullcheck-part "$filt_part" --out "$FULLCHECK_OUT/filtered" \
    --layouts "$LAYOUT" --operons "$OPERONS" \
    --cohort-metadata "$COHORT_METADATA" --jar "$JAR" --reference "$REFERENCE" \
    --models "$MODELS" --batch-threads "$BATCH_THREADS" --heap "$FULLCHECK_HEAP" \
    --anomaly-floor "$ANOMALY_FLOOR" --anomaly-max-candidates "$ANOMALY_MAX_CANDIDATES" \
    --anomaly-max-locus "$ANOMALY_MAX_LOCUS" \
    > "$jobout/fullcheck-filtered.log" 2>&1 || { tail -n 40 "$jobout/fullcheck-filtered.log"; exit 1; }

  "${prefix[@]}" --env "TMPDIR=$work/tmp" "$IMAGE" "$PYTHON" \
    "$REPO_ROOT/benchmark/v3/row11_ec_go/compare_full_vs_filtered.py" \
    --full "$FULLCHECK_OUT/full/stage1/$FULLCHECK_ASSEMBLY.integrated.tsv" \
    --filtered "$FULLCHECK_OUT/filtered/stage1/$FULLCHECK_ASSEMBLY.integrated.tsv" \
    --min-raw "$MIN_RAW" \
    --claims-full "$full_part/claims.jsonl" \
    --out "$FULLCHECK_OUT/report.json" > "$jobout/fullcheck-compare.log" 2>&1 || { tail -n 40 "$jobout/fullcheck-compare.log"; exit 1; }

  tail -n 20 "$jobout/fullcheck-compare.log"
  rm -rf "$work"; date --iso-8601=seconds > "$jobout/finished-at.txt"; exit 0
fi

echo "Unknown MODE: $MODE" >&2
exit 2
