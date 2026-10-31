#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=gspa-price-rc-evaluate
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --exclude=nid00022,nid00497,nid01033
#SBATCH --time=12:00:00
#SBATCH --no-requeue
#SBATCH --output=/scratch/project/k10018/gspa-v3/out/price-reference-rc1/slurm/evaluate-%j.log
#
# Shaheen III evaluation of the release-candidate Price collection on 1 exclusive workq node.
# Equivalent of unimatrix price-reference-full-launch-v1/run.sbatch (job 10432) on the pinned
# Shaheen release-candidate collection ($OUT_ROOT/collection).
#
# Evaluation workflow:
#   1. Check archive package.sha256 and source/source.sha256 hashes byte-for-byte.
#   2. Check collection prerequisites (prerequisites_shaheen.py) verifying registration,
#      plan, and completed/verified collection artifacts.
#   3. Run metric and evaluator unit tests (test_price_reference_metrics.py).
#   4. Parallel per-predictor evaluation: runs the 6 predictors (topec, mapred, enzbert,
#      dpz, graphec, clean) concurrently across node cores into separate out dirs
#      ($jobout/evaluation-<predictor>) with per-predictor /usr/bin/time -v and independent
#      arithmetic checks ($jobout/independent-arithmetic-<predictor>.json).
#   5. Joint panel evaluation: runs evaluate_price_reference.py once across all six predictors
#      into $jobout/evaluation, creating the unified 132-cell report.json, followed by the
#      canonical independent arithmetic check ($jobout/independent-arithmetic.json).
#   6. Full receipts: executed.sbatch, started-at.txt, finished-at.txt, exit-code.txt, hostname.txt.
#
# Parameters (environment):
#   OUT_ROOT       - Root of collection output (default: /scratch/project/k10018/gspa-v3/out/price-production-rc1)
#   EVAL_ROOT      - Root of evaluation staging (default: /scratch/project/k10018/gspa-v3/out/price-reference-rc1)
#   COLLECTION     - Path to collection (default: $OUT_ROOT/collection)
#   FINALIZER_DIR  - Optional path to finalize job directory (default: auto-detected under $OUT_ROOT/jobs)
#   RUN_PARALLEL   - 1 = run 6 individual predictor evaluations in parallel (default: 1)
#   RUN_JOINT      - 1 = run full joint 6-predictor evaluation (default: 1)
#   DRY_RUN        - 1 = print container commands without execution (default: 0)
#   PYTHON         - Python interpreter (default: /project/k10018/sw/envs/gspa-tools/bin/python3.12)
#   PYLIB          - Python site-packages (default: /data/hohndor/gspa-v3/price-production-v1/venv/lib/python3.12/site-packages)
#   FSROOT         - Shared package tree (default: /scratch/project/k10018/gspa-v3/shared/arm-fs-d04c9e11ee64db47)
#   IMAGE          - Singularity image (default: /project/k10018/sw/images/debian-trixie)
#   WORK_ROOT      - Scratch work root (default: /scratch/project/k10018/gspa-v3/work)
set -euo pipefail
unset JAVA_TOOL_OPTIONS JDK_JAVA_OPTIONS _JAVA_OPTIONS PYTHONPATH PYTHONHOME

OUT_ROOT=${OUT_ROOT:-/scratch/project/k10018/gspa-v3/out/price-production-rc1}
EVAL_ROOT=${EVAL_ROOT:-/scratch/project/k10018/gspa-v3/out/price-reference-rc1}
COLLECTION=${COLLECTION:-$OUT_ROOT/collection}
PYTHON=${PYTHON:-/project/k10018/sw/envs/gspa-tools/bin/python3.12}
PYLIB=${PYLIB:-/data/hohndor/gspa-v3/price-production-v1/venv/lib/python3.12/site-packages}
FSROOT=${FSROOT:-/scratch/project/k10018/gspa-v3/shared/arm-fs-d04c9e11ee64db47}
IMAGE=${IMAGE:-/project/k10018/sw/images/debian-trixie}
WORK_ROOT=${WORK_ROOT:-/scratch/project/k10018/gspa-v3/work}
FINALIZER_DIR=${FINALIZER_DIR:-}
RUN_PARALLEL=${RUN_PARALLEL:-1}
RUN_JOINT=${RUN_JOINT:-1}
DRY_RUN=${DRY_RUN:-0}
KEEP_WORK=${KEEP_WORK:-0}
PREDICTORS=(topec mapred enzbert dpz graphec clean)
logical_data=/data/hohndor/gspa-v3

JOB=${SLURM_JOB_ID:-dry$$}
TMPDIR=/dev/shm/$JOB
export TMPDIR

if [ "$DRY_RUN" = 1 ]; then
  launcher=/usr/bin/singularity
  jobout=${DRY_OUT:-${TMPDIR:?}/price-rc-dry/$JOB}
  work=$jobout/work
else
  module load singularity
  launcher=$(readlink -f "$(command -v singularity)")
  [ -d "$EVAL_ROOT" ] || { echo "EVAL_ROOT missing: $EVAL_ROOT" >&2; exit 2; }
  [ -d "$COLLECTION" ] || { echo "Collection missing: $COLLECTION" >&2; exit 2; }
  [ -f "$FSROOT/.complete" ] || { echo "shared package tree missing: $FSROOT" >&2; exit 2; }
  jobout=$EVAL_ROOT/jobs/$JOB
  work=$WORK_ROOT/$JOB
fi
mkdir -p "$jobout" "$TMPDIR" "$work/tmp"

finish() {
  code=$?
  trap - EXIT
  printf '%s\n' "$code" > "$jobout/exit-code.txt"
  date --iso-8601=seconds > "$jobout/finished-at.txt"
  [ -d "$TMPDIR" ] && rm -rf "$TMPDIR"
  [ "$KEEP_WORK" = 1 ] || [ "$DRY_RUN" = 1 ] || rm -rf "$work"
  exit "$code"
}
trap finish EXIT

date --iso-8601=seconds > "$jobout/started-at.txt"
hostname > "$jobout/hostname.txt"
cp "$0" "$jobout/executed.sbatch"
if [ "$DRY_RUN" != 1 ]; then
  singularity --version > "$jobout/singularity-version.txt"
fi

# Container launcher matching price-rc-register and price-rc-cells
container=("$launcher" exec --no-home --no-mount tmp
  --bind /project/k10018/sw/envs/gspa-tools:/project/k10018/sw/envs/gspa-tools:ro
  --bind "/project/k10018/gspa-v3:$logical_data:ro"
  --bind "$OUT_ROOT:$OUT_ROOT"
  --bind "$EVAL_ROOT:$EVAL_ROOT"
  --bind "$work:$work"
  --bind "$TMPDIR:$TMPDIR"
  --bind "$FSROOT/usr/lib/jvm/java-21-openjdk-amd64:/usr/lib/jvm/java-21-openjdk-amd64:ro"
  --bind "$FSROOT/home/leechuck:/home/leechuck:ro"
  --bind /usr/bin/time:/usr/bin/time:ro
  --env PATH=/usr/lib/jvm/java-21-openjdk-amd64/bin:/usr/bin:/bin
  --env LC_ALL=C.UTF-8
  --env PYTHONDONTWRITEBYTECODE=1
  --env "PYTHONPATH=$PYLIB"
  --env "TMPDIR=$TMPDIR")
# Options (e.g. --pwd, --env) must come before the image: use "${cbase[@]}" <options> "$IMAGE" <command>.
cbase=("${container[@]}")
container=("${cbase[@]}" "$IMAGE")
printf '%s\n' "${container[@]}" > "$jobout/container-command.txt"

run() {
  if [ "$DRY_RUN" = 1 ]; then
    printf '%q ' "$@" >> "$jobout/dry-run.txt"
    echo >> "$jobout/dry-run.txt"
  else
    "$@"
  fi
}

# Step 1: Package and source checksum verification
run "${cbase[@]}" --pwd "$EVAL_ROOT" "$IMAGE" sha256sum -c package.sha256 > "$jobout/package-check.log"
run "${cbase[@]}" --pwd "$EVAL_ROOT/source" "$IMAGE" sha256sum -c "$EVAL_ROOT/source.sha256" > "$jobout/source-check.log"

# Step 2: Prerequisites check
# We use prerequisites_shaheen.py which accepts parameterised collection and finalizer paths.
# PREREQUISITES_SCRIPT environment variable allows overriding if needed.
if [ -n "${PREREQUISITES_SCRIPT:-}" ]; then
  run "${container[@]}" "$PYTHON" "$PREREQUISITES_SCRIPT" \
    --collection "$COLLECTION" \
    ${FINALIZER_DIR:+--finalizer "$FINALIZER_DIR"} \
    --out "$jobout/prerequisites.json" > "$jobout/prerequisites.log" 2>&1
elif [ -f "$EVAL_ROOT/prerequisites.py" ] && [ ! -f "$EVAL_ROOT/prerequisites_shaheen.py" ] && [ "$DRY_RUN" != 1 ]; then
  run "${container[@]}" "$PYTHON" "$EVAL_ROOT/prerequisites.py" \
    "$jobout/prerequisites.json" > "$jobout/prerequisites.log" 2>&1
else
  run "${container[@]}" "$PYTHON" "$EVAL_ROOT/prerequisites_shaheen.py" \
    --collection "$COLLECTION" \
    ${FINALIZER_DIR:+--finalizer "$FINALIZER_DIR"} \
    --out "$jobout/prerequisites.json" > "$jobout/prerequisites.log" 2>&1
fi

# Step 3: Unittest metrics suite
run "${cbase[@]}" --pwd "$EVAL_ROOT/source" \
  --env PYTHONPATH="$EVAL_ROOT/source/benchmark/v3:$PYLIB" "$IMAGE" \
  "$PYTHON" -m unittest discover -s "$EVAL_ROOT/source/benchmark/v3/tests" \
  -p test_price_reference_metrics.py -v > "$jobout/tests.log" 2>&1

# Step 4: Parallel per-predictor evaluation and arithmetic check
# Running the 6 predictors concurrently across node cores:
# - Independent per-predictor evaluation and arithmetic receipts
# - Faster wall-clock execution by reading disjoint cell subsets across Lustre
# - Clear failure isolation if a specific predictor or cell has issues
if [ "$RUN_PARALLEL" = 1 ]; then
  for p in "${PREDICTORS[@]}"; do
    eval_p_out="$jobout/evaluation-$p"
    if [ "$DRY_RUN" = 1 ]; then
      run /usr/bin/time -v -o "$jobout/resources-$p.txt" "${container[@]}" "$PYTHON" \
        "$EVAL_ROOT/source/benchmark/v3/evaluate_price_reference.py" \
        --collection "$COLLECTION" \
        --manifest "$EVAL_ROOT/source/benchmark/manifests/suite6_price149_manifest.tsv" \
        --ec-audit "$EVAL_ROOT/ec-audit.json" \
        --predictors "$p" \
        --out "$eval_p_out"
      run "${container[@]}" "$PYTHON" "$EVAL_ROOT/check_retained.py" \
        --root "$eval_p_out" \
        --audit "$EVAL_ROOT/ec-audit.json" \
        --manifest "$EVAL_ROOT/source/benchmark/manifests/suite6_price149_manifest.tsv"
    else
      (
        /usr/bin/time -v -o "$jobout/resources-$p.txt" "${container[@]}" "$PYTHON" \
          "$EVAL_ROOT/source/benchmark/v3/evaluate_price_reference.py" \
          --collection "$COLLECTION" \
          --manifest "$EVAL_ROOT/source/benchmark/manifests/suite6_price149_manifest.tsv" \
          --ec-audit "$EVAL_ROOT/ec-audit.json" \
          --predictors "$p" \
          --out "$eval_p_out" > "$jobout/evaluation-$p.log" 2>&1
        "${container[@]}" "$PYTHON" "$EVAL_ROOT/check_retained.py" \
          --root "$eval_p_out" \
          --audit "$EVAL_ROOT/ec-audit.json" \
          --manifest "$EVAL_ROOT/source/benchmark/manifests/suite6_price149_manifest.tsv" \
          > "$jobout/independent-arithmetic-$p.json" 2> "$jobout/check-retained-$p.log"
      ) &
    fi
  done
  if [ "$DRY_RUN" != 1 ]; then
    for pid in $(jobs -p); do
      wait "$pid"
    done
  fi
fi

# Step 5: Full joint panel evaluation and arithmetic check
# Evaluates all six predictors together, mirroring the unimatrix benchmark job:
# produces unified $jobout/evaluation/report.json and $jobout/independent-arithmetic.json
# covering all 132 cells under a single registration.
if [ "$RUN_JOINT" = 1 ]; then
  eval_out="$jobout/evaluation"
  run /usr/bin/time -v -o "$jobout/resources.txt" "${container[@]}" "$PYTHON" \
    "$EVAL_ROOT/source/benchmark/v3/evaluate_price_reference.py" \
    --collection "$COLLECTION" \
    --manifest "$EVAL_ROOT/source/benchmark/manifests/suite6_price149_manifest.tsv" \
    --ec-audit "$EVAL_ROOT/ec-audit.json" \
    --predictors "${PREDICTORS[@]}" \
    --out "$eval_out" > "$jobout/evaluation.log" 2>&1
  run "${container[@]}" "$PYTHON" "$EVAL_ROOT/check_retained.py" \
    --root "$eval_out" \
    --audit "$EVAL_ROOT/ec-audit.json" \
    --manifest "$EVAL_ROOT/source/benchmark/manifests/suite6_price149_manifest.tsv" \
    > "$jobout/independent-arithmetic.json" 2> "$jobout/check-retained.log"
fi

if [ "$DRY_RUN" = 1 ]; then
  echo "=== DRY RUN: container commands logged to $jobout/dry-run.txt ==="
  cat "$jobout/dry-run.txt"
  exit 0
fi
