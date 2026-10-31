#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=gspa-price-stage2b-cells
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --no-requeue
#SBATCH --exclude=nid00022,nid00497,nid01033
#SBATCH --output=/scratch/project/k10018/gspa-v3/out/price-stage2b-rc1/slurm/cells-%j.log
#
# Shaheen III Stage 2B (flux-certified refiner) Price 132-cell collection script.
# Arms:
#   - stage0:  unchanged stage0 command (--enable-priors none)
#   - stage2b: stage1 command (default priors) + --enable-fba-refiner --fba-data-dir <dir>
#              --fba-gram neg --fba-medium minimal --fba-python-bin /project/k10018/sw/envs/gspa-v3/bin/python
#              --fba-sidecar $OUT_ROOT/code/scripts/stoichiometric_locus_refiner.py
#
# Behavior when invoked outside Slurm (SLURM_JOB_ID unset, RUN_PAYLOAD!=1, DRY_RUN!=1):
#   - Prints the Slurm submission schedule (MODE=register, 51 packed cell jobs, MODE=finalize) by default.
#   - Submits only when `--submit` is explicitly passed.
#
# Modes inside Slurm (or with RUN_PAYLOAD=1 / DRY_RUN=1):
#   - MODE=register: builds $OUT_ROOT/fba_data from the pinned gapsmith reference archive, snapshots code/JAR,
#                    runs unit tests, and registers $OUT_ROOT/collection via run_price_stage2b_collection.py prepare.
#   - MODE=cells:    executes registered cells in CELLS using the 3-step Singularity external launch
#                    (prepare-external -> host exec of recorded locator argv -> publish-external) with
#                    memory-class admission matching price-rc-cells.sbatch.
#   - MODE=finalize: verifies preflight hashes and runs run_price_stage2b_collection.py finalize.

set -euo pipefail

# Under Slurm, $0 is the spooled copy (/var/spool/...), so the checkout is passed explicitly (REPO_ROOT, exported by the
# launcher) and only derived from the script location when run directly.
SCRIPT_PATH=$(readlink -f "$0")
if [ -z "${REPO_ROOT:-}" ]; then
  HERE=$(cd "$(dirname "$SCRIPT_PATH")" && pwd)
  REPO_ROOT=$(cd "$HERE/../../.." && pwd)
fi
HERE=$REPO_ROOT/benchmark/v3/shaheen
SCRIPT_PATH=$HERE/price-stage2b-cells.sbatch
[ -f "$REPO_ROOT/benchmark/v3/run_price_local.py" ] || { echo "REPO_ROOT has no checkout: $REPO_ROOT" >&2; exit 2; }

OUT_ROOT=${OUT_ROOT:-/scratch/project/k10018/gspa-v3/out/price-stage2b-rc1}
JAR=${JAR:-/project/k10018/gspa-v3/release-candidate-1/gspa-3.0.0.jar}
FROZEN=${FROZEN:-/project/k10018/gspa-v3/price-production-v3}
PYTHON=${PYTHON:-/project/k10018/sw/envs/gspa-tools/bin/python3.12}
PYLIB=${PYLIB:-/data/hohndor/gspa-v3/price-production-v1/venv/lib/python3.12/site-packages}
FBA_PYTHON_BIN=${FBA_PYTHON_BIN:-/project/k10018/sw/envs/gspa-v3/bin/python}
GAPSMITH_SOURCE=${GAPSMITH_SOURCE:-/project/k10018/sw/gapsmith-eq-ibex/data_merged.tar}
FBA_DATA_DIR=${FBA_DATA_DIR:-$OUT_ROOT/fba_data}
FBA_GRAM=${FBA_GRAM:-neg}
FBA_MEDIUM=${FBA_MEDIUM:-minimal}
FSROOT=${FSROOT:-/scratch/project/k10018/gspa-v3/shared/arm-fs-d04c9e11ee64db47}
IMAGE=${IMAGE:-/project/k10018/sw/images/debian-trixie}
WORK_ROOT=${WORK_ROOT:-/scratch/project/k10018/gspa-v3/work}
MEM_BUDGET_G=${MEM_BUDGET_G:-340}
HEAP_SMALL=${HEAP_SMALL:-24g} HEAP_MEDIUM=${HEAP_MEDIUM:-96g} HEAP_LARGE=${HEAP_LARGE:-168g}
ALLOC_SMALL_G=${ALLOC_SMALL_G:-36} ALLOC_MEDIUM_G=${ALLOC_MEDIUM_G:-112} ALLOC_LARGE_G=${ALLOC_LARGE_G:-192}
TIMEOUT_SMALL=${TIMEOUT_SMALL:-7200} TIMEOUT_MEDIUM=${TIMEOUT_MEDIUM:-14400} TIMEOUT_LARGE=${TIMEOUT_LARGE:-21600}
PIN=${PIN:-1} DRY_RUN=${DRY_RUN:-0} KEEP_WORK=${KEEP_WORK:-0}

if [ -z "${SLURM_JOB_ID:-}" ] && [ "${RUN_PAYLOAD:-0}" != "1" ] && [ "$DRY_RUN" != "1" ]; then
  submit=0
  case "${1:-}" in
    --submit) submit=1 ;;
    "") ;;
    *) echo "Usage: $0 [--submit]" >&2; exit 2 ;;
  esac

  LARGE_PER_NODE=${LARGE_PER_NODE:-1}
  MEDIUM_PER_NODE=${MEDIUM_PER_NODE:-3}
  SMALL_PER_NODE=${SMALL_PER_NODE:-9}
  CELL_TIME=${CELL_TIME:-12:00:00}
  SMALL_TIME=${SMALL_TIME:-06:00:00}
  FINAL_TIME=${FINAL_TIME:-24:00:00}
  common=(--account=k10018 --partition=workq --nodes=1 --exclude=nid00022,nid00497,nid01033 --no-requeue)
  pass=(REPO_ROOT OUT_ROOT JAR FROZEN PYTHON PYLIB FBA_PYTHON_BIN GAPSMITH_SOURCE FBA_DATA_DIR FBA_GRAM FBA_MEDIUM
        FSROOT IMAGE WORK_ROOT MEM_BUDGET_G HEAP_SMALL HEAP_MEDIUM HEAP_LARGE ALLOC_SMALL_G ALLOC_MEDIUM_G
        ALLOC_LARGE_G TIMEOUT_SMALL TIMEOUT_MEDIUM TIMEOUT_LARGE PIN KEEP_WORK)
  envs=()
  for v in "${pass[@]}"; do [ -z "${!v:-}" ] || envs+=("$v=${!v}"); done

  lines=()
  emit() { lines+=("$1"); }
  q() { [ "$#" -eq 0 ] || printf '%q ' "$@"; }
  groups() {
    local per=$1 n=0 g=(); shift
    for i in "$@"; do g+=("$i"); n=$((n+1)); if [ "$n" -eq "$per" ]; then echo "${g[*]}"; g=(); n=0; fi; done
    [ "${#g[@]}" -eq 0 ] || echo "${g[*]}"
  }

  emit "set -euo pipefail"
  emit "mkdir -p $(q "$OUT_ROOT/slurm")"
  emit "test ! -e $(q "$OUT_ROOT/collection")"
  emit "REGISTER=\$(env $(q "${envs[@]}")MODE=register sbatch --parsable $(q "${common[@]}") --time=02:00:00 --job-name=gspa-price-stage2b-register --output=$(q "$OUT_ROOT/slurm/register-%j.log") $(q "$SCRIPT_PATH"))"
  emit "echo register=\$REGISTER"
  emit "CELL_JOBS=''"
  k=0
  while read -r cls per t; do
    case $cls in
      large) idx=$(seq 131 -1 101) ;;
      medium) idx=$(seq 100 -1 66) ;;
      small) idx=$(seq 65 -1 0) ;;
    esac
    # shellcheck disable=SC2086
    while read -r group; do
      k=$((k+1))
      emit "J=\$(env $(q "${envs[@]}")MODE=cells CELLS=$(q "$group") sbatch --parsable $(q "${common[@]}") --time=$t --dependency=afterok:\$REGISTER --job-name=gspa-price-stage2b-$cls --output=$(q "$OUT_ROOT/slurm/cells-%j.log") $(q "$SCRIPT_PATH")); CELL_JOBS=\$CELL_JOBS:\$J; echo \"$cls cells [$group] job=\$J\""
    done < <(groups "$per" $idx)
  done <<EOF
large $LARGE_PER_NODE $CELL_TIME
medium $MEDIUM_PER_NODE $CELL_TIME
small $SMALL_PER_NODE $SMALL_TIME
EOF
  emit "FINAL=\$(env $(q "${envs[@]}")MODE=finalize sbatch --parsable $(q "${common[@]}") --time=$FINAL_TIME --dependency=afterok\$CELL_JOBS --job-name=gspa-price-stage2b-finalize --output=$(q "$OUT_ROOT/slurm/finalize-%j.log") $(q "$SCRIPT_PATH"))"
  emit "echo finalize=\$FINAL"
  emit "printf '%s\\n' \"register \$REGISTER\" \"cells \${CELL_JOBS#:}\" \"finalize \$FINAL\" > \"$OUT_ROOT/slurm/schedule-\$REGISTER.txt\""

  if [ "$submit" = "1" ]; then
    command -v sbatch >/dev/null || { echo "sbatch not found: run on a Shaheen login node" >&2; exit 2; }
    printf '%s\n' "${lines[@]}" | bash
  else
    echo "# Dry-run mode (default): $k cell jobs + register + finalize for Price Stage 2B; pass --submit to submit"
    printf '%s\n' "${lines[@]}"
  fi
  exit 0
fi

MODE=${MODE:-cells}
JOB=${SLURM_JOB_ID:-dry$$}
logical_data=/data/hohndor/gspa-v3
collection=$OUT_ROOT/collection
jobout=$OUT_ROOT/jobs/$JOB
work=$WORK_ROOT/$JOB
unset JAVA_TOOL_OPTIONS JDK_JAVA_OPTIONS _JAVA_OPTIONS PYTHONPATH PYTHONHOME

FROZEN_CODE_HASHES_SHA256=cfb8ac69520d2f1df07849cd0825ffe0224e97b9f29a275cde1779490968a7c4
ORIGINAL_CELLS=4bc45a18801dc7bc2d1ef4df46aad44f82fd267be12d2312d6ebf5a41aa486bb
ORIGINAL_DATA_INPUTS=bcadd3959fd0ba06896df6eef4c88ebeb8f446ae6d51b295e5a88fb1aa989937
ENGINE=(run_price_matrix_capacity.py run_price_production_collection.py stage_ledger.py v3_metrics.py v3_runner.py
        verified_batch.py verify_candidate_universe.py requirements-streaming.txt)

if [ "$DRY_RUN" = "1" ]; then
  launcher=/usr/bin/singularity
  jobout=${DRY_OUT:-${TMPDIR:?}/price-stage2b-dry/$JOB}; work=$jobout/work
else
  module load singularity
  launcher=$(readlink -f "$(command -v singularity)")
  [ -f "$FSROOT/.complete" ] || { echo "shared package tree missing: $FSROOT" >&2; exit 2; }
fi
mkdir -p "$jobout" "$work/tmp"
date --iso-8601=seconds > "$jobout/started-at.txt"; hostname > "$jobout/hostname.txt"
cp "$SCRIPT_PATH" "$jobout/executed.sbatch"

# Host launcher argv up to (excluding) the image. Binds /project/k10018/sw/envs/gspa-v3 for --fba-python-bin
# and enforces single-threaded BLAS inside the container.
prefix=("$launcher" exec --no-home --no-mount tmp
  --bind /project/k10018/sw/envs/gspa-tools:/project/k10018/sw/envs/gspa-tools:ro
  --bind /project/k10018/sw/envs/gspa-v3:/project/k10018/sw/envs/gspa-v3:ro
  --bind "/project/k10018/gspa-v3:$logical_data:ro" --bind "$OUT_ROOT:$OUT_ROOT" --bind "$work:$work"
  --bind "$FSROOT/usr/lib/jvm/java-21-openjdk-amd64:/usr/lib/jvm/java-21-openjdk-amd64:ro"
  --bind "$FSROOT/home/leechuck:/home/leechuck:ro"
  --bind /usr/bin/time:/usr/bin/time:ro
  --env PATH=/usr/lib/jvm/java-21-openjdk-amd64/bin:/usr/bin:/bin
  --env LC_ALL=C.UTF-8 --env PYTHONDONTWRITEBYTECODE=1
  --env OPENBLAS_NUM_THREADS=1 --env OMP_NUM_THREADS=1 --env MKL_NUM_THREADS=1
  --env "PYTHONPATH=$PYLIB")
wrap=("${prefix[@]}" --bind "$launcher:$launcher:ro" --bind "$IMAGE:$IMAGE:ro")
/usr/bin/python3 -c 'import json,sys; print(json.dumps(sys.argv[1:]))' "${prefix[@]}" > "$work/prefix.json"
cp "$work/prefix.json" "$jobout/prefix.json"
run() { if [ "$DRY_RUN" = "1" ]; then printf '%q ' "$@" >> "$jobout/dry-run.txt"; echo >> "$jobout/dry-run.txt"; else "$@"; fi; }

if [ "$MODE" = "register" ]; then
  [ ! -e "$OUT_ROOT/collection" ] || { echo "refusing: $OUT_ROOT/collection already exists" >&2; exit 2; }
  [ ! -e "$OUT_ROOT/code" ] || { echo "refusing: $OUT_ROOT/code already exists" >&2; exit 2; }
  mkdir -p "$OUT_ROOT/code/stage2b" "$OUT_ROOT/code/scripts" "$OUT_ROOT/code/tests"

  if [ "$DRY_RUN" != "1" ]; then
    echo "$FROZEN_CODE_HASHES_SHA256  $FROZEN/code-hashes.json" | sha256sum -c - > "$jobout/frozen-code-check.log"
    for f in "${ENGINE[@]}"; do cp "$FROZEN/code/$f" "$OUT_ROOT/code/$f"; done
    cp "$FROZEN/code-hashes.json" "$OUT_ROOT/code/frozen-code-hashes.json"
    cp "$REPO_ROOT/benchmark/v3/run_price_local.py" "$OUT_ROOT/code/run_price_local.py"
    cp "$REPO_ROOT/benchmark/v3/price_reference_metrics.py" "$OUT_ROOT/code/price_reference_metrics.py"
    cp "$REPO_ROOT/benchmark/v3/evaluate_price_reference.py" "$OUT_ROOT/code/evaluate_price_reference.py"
    cp "$REPO_ROOT/benchmark/v3/ec_identity.py" "$OUT_ROOT/code/ec_identity.py"
    cp "$REPO_ROOT/benchmark/v3/stage2b/__init__.py" "$OUT_ROOT/code/stage2b/__init__.py"
    cp "$REPO_ROOT/benchmark/v3/stage2b/build_fba_data.py" "$OUT_ROOT/code/stage2b/build_fba_data.py"
    cp "$REPO_ROOT/benchmark/v3/stage2b/run_price_stage2b_collection.py" "$OUT_ROOT/code/stage2b/run_price_stage2b_collection.py"
    cp "$REPO_ROOT/benchmark/v3/stage2b/evaluate_price_stage2b.py" "$OUT_ROOT/code/stage2b/evaluate_price_stage2b.py"
    cp "$REPO_ROOT/scripts/stoichiometric_locus_refiner.py" "$OUT_ROOT/code/scripts/stoichiometric_locus_refiner.py"
    cp "$REPO_ROOT/benchmark/v3/tests/test_stage2b_build_fba_data.py" "$OUT_ROOT/code/tests/test_stage2b_build_fba_data.py"
    cp "$REPO_ROOT/benchmark/v3/tests/test_stage2b_collection_and_eval.py" "$OUT_ROOT/code/tests/test_stage2b_collection_and_eval.py"
    cp "$JAR" "$OUT_ROOT/gspa-3.0.0.jar"
    chmod -R a-w "$OUT_ROOT/code" "$OUT_ROOT/gspa-3.0.0.jar"
  fi

  if [ ! -f "$FBA_DATA_DIR/universal_metadata.json" ]; then
    run /usr/bin/time -v -o "$jobout/build-fba-data.resources.txt" \
      env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 TMPDIR="$work/tmp" \
      "$FBA_PYTHON_BIN" "$OUT_ROOT/code/stage2b/build_fba_data.py" \
      --gapsmith-source "$GAPSMITH_SOURCE" \
      --output-dir "$FBA_DATA_DIR" \
      --verify-pinned-sha256 > "$jobout/build-fba-data.log" 2>&1
  fi

  testtmp=/dev/shm/price-stage2b-tests-$JOB
  run mkdir -p "$testtmp"
  # The tests expect the repository layout (manifests, scripts/); run them from the checkout the snapshot was copied from.
  run env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 TMPDIR="$testtmp" PYTHONDONTWRITEBYTECODE=1 \
    PYTHONPATH="$REPO_ROOT/benchmark/v3" \
    "$FBA_PYTHON_BIN" -m unittest discover -s "$REPO_ROOT/benchmark/v3/tests" -p "test_stage2b_*.py" -v \
    > "$jobout/tests.log" 2>&1

  run /usr/bin/time -v -o "$jobout/prepare.resources.txt" \
    "${prefix[@]}" --pwd "$OUT_ROOT" --env "TMPDIR=$work/tmp" "$IMAGE" "$PYTHON" \
    "$OUT_ROOT/code/stage2b/run_price_stage2b_collection.py" prepare \
    --out "$OUT_ROOT/collection" \
    --export "$logical_data/price-full-matrix-v2/plan.json" \
    --audit "$logical_data/price-full-collection-audit-v1/verification.json" \
    --jar "$OUT_ROOT/gspa-3.0.0.jar" \
    --reference "$logical_data/price-production-v1/reference" \
    --fba-data-dir "$FBA_DATA_DIR" \
    --fba-gram "$FBA_GRAM" \
    --fba-medium "$FBA_MEDIUM" \
    --fba-python-bin "$FBA_PYTHON_BIN" \
    --fba-sidecar "$OUT_ROOT/code/scripts/stoichiometric_locus_refiner.py" \
    > "$jobout/prepare.log" 2>&1

  export OUT_ROOT ORIGINAL_CELLS ORIGINAL_DATA_INPUTS JOB
  run "${prefix[@]}" --pwd "$OUT_ROOT" --env "TMPDIR=$work/tmp" "$IMAGE" "$PYTHON" - > "$jobout/compare-original.log" 2>&1 <<'PY'
import hashlib, json, os
from pathlib import Path
out = Path(os.environ['OUT_ROOT']); data = '/data/hohndor/gspa-v3/'
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
digest = lambda x: hashlib.sha256(json.dumps(x, sort_keys=True).encode()).hexdigest()
plan = json.loads((out/'collection/plan.json').read_text())
registration = json.loads((out/'collection/registration.json').read_text())
inputs = registration['inputs']
assert len(plan['cells']) == 132 and plan['expected_claims'] == 1526919081
assert plan['arms'] == ['stage0', 'stage2b']
assert plan['taxonomy_verification']['all_gram_negative'] is True
assert plan['taxonomy_verification']['all_pseudomonadota'] is True
assert digest(plan['cells']) == os.environ['ORIGINAL_CELLS'], 'plan cells differ from original collection'
prefixes = tuple(data+p for p in ('price-full-matrix-v2/', 'price-full-collection-audit-v1/', 'price-production-v1/reference/'))
assert digest({k: v for k, v in inputs.items() if k.startswith(prefixes)}) == os.environ['ORIGINAL_DATA_INPUTS'], 'data inputs differ'
jar = str(out/'gspa-3.0.0.jar')
assert plan['jar'] == jar and inputs.get(jar) == sha(jar)
receipt = dict(
    schema='gspa.price-stage2b-registration/1',
    job=os.environ['JOB'],
    plan_sha256=registration['plan_sha256'],
    registration_sha256=sha(out/'collection/registration.json'),
    jar=plan['jar'],
    jar_sha256=inputs.get(jar),
    fba_data_dir=plan['fba_data_dir'],
    fba_gram=plan['fba_gram'],
    fba_medium=plan['fba_medium'],
    fba_python_bin=plan['fba_python_bin'],
    fba_sidecar=plan['fba_sidecar'],
    original_plan_cells_identical=True,
    original_data_inputs_identical=True,
)
with (out/'register-receipt.json').open('x') as stream:
    json.dump(receipt, stream, indent=2, sort_keys=True); stream.write('\n')
print(json.dumps(receipt, indent=2, sort_keys=True))
PY
  [ "$DRY_RUN" = "1" ] && exit 0
  cp "$OUT_ROOT/register-receipt.json" "$jobout/"
  rm -rf "$work" "$testtmp"
  exit 0
fi

if [ "$DRY_RUN" != "1" ]; then
  [ -f "$collection/registration.json" ] && [ -f "$OUT_ROOT/register-receipt.json" ] || { echo "collection not registered: $collection" >&2; exit 2; }
fi

if [ "$MODE" = "finalize" ]; then
  export OUT_ROOT
  run "${prefix[@]}" --env "TMPDIR=$work/tmp" "$IMAGE" "$PYTHON" - > "$jobout/preflight.log" 2>&1 <<'PY'
import hashlib, json, os
from pathlib import Path
out = Path(os.environ['OUT_ROOT']); sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
receipt = json.loads((out/'register-receipt.json').read_text())
registration = json.loads((out/'collection/registration.json').read_text())
assert sha(out/'collection/registration.json') == receipt['registration_sha256']
assert sha(out/'collection/plan.json') == registration['plan_sha256'] == receipt['plan_sha256']
assert not (out/'collection/collection-verification.json').exists()
print('Stage 2B registration and plan verified; jar', receipt['jar_sha256'])
PY
  code=0
  run /usr/bin/time -v -o "$jobout/resources.txt" "${prefix[@]}" --env "TMPDIR=$work/tmp" "$IMAGE" "$PYTHON" \
    "$OUT_ROOT/code/stage2b/run_price_stage2b_collection.py" finalize --out "$collection" > "$jobout/finalize.log" 2>&1 || code=$?
  printf '%s\n' "$code" > "$jobout/exit-code.txt"; date --iso-8601=seconds > "$jobout/finished-at.txt"
  [ "$code" -eq 0 ] || exit "$code"
  [ "$DRY_RUN" = "1" ] && exit 0
  cp "$collection/collection-verification.json" "$jobout/"
  sha256sum "$jobout/collection-verification.json" > "$jobout/receipt.sha256"
  rm -rf "$work"; exit 0
fi

: "${CELLS:?cell indices required}"
read -r -a cells <<< "${CELLS//[,+]/ }"
class_of() { if [ "$1" -le 65 ]; then echo SMALL; elif [ "$1" -le 100 ]; then echo MEDIUM; else echo LARGE; fi; }
total=0
for i in "${cells[@]}"; do
  [[ "$i" =~ ^[0-9]+$ ]] && [ "$i" -le 131 ] || { echo "invalid cell index: $i" >&2; exit 2; }
  c=$(class_of "$i"); a=ALLOC_${c}_G; total=$((total + ${!a}))
done
declare -A cpus=()
if [ "$PIN" = "1" ] && [ "$total" -le "$MEM_BUDGET_G" ] && [ "${#cells[@]}" -gt 1 ]; then
  mapfile -t cores < <(lscpu -p=CORE | grep -v '^#' | sort -n -u)
  per=$(( ${#cores[@]} / ${#cells[@]} ))
  for k in "${!cells[@]}"; do
    lo=${cores[$((k*per))]}; hi=${cores[$(((k+1)*per-1))]}
    cpus[${cells[$k]}]=$(lscpu -p=CPU,CORE | grep -v '^#' | awk -F, -v lo="$lo" -v hi="$hi" '$2>=lo && $2<=hi {print $1}' | paste -sd,)
  done
fi
echo "job $JOB on $(hostname): cells=[${cells[*]}] alloc=${total}G budget=${MEM_BUDGET_G}G pin=${#cpus[@]}"

run_cell() {
  local i=$1 c heap timeout attempt cw log locator rc=0 prc=0 affinity=()
  c=$(class_of "$i"); heap=HEAP_$c; timeout=TIMEOUT_$c; heap=${!heap}; timeout=${!timeout}
  attempt=job-$JOB; cw=$work/cell-$i; log=$jobout/cell-$i
  mkdir -p "$cw/tmp" "$log"
  date --iso-8601=seconds > "$log/started-at.txt"
  echo "cell $i class=$c heap=$heap timeout=$timeout cpus=${cpus[$i]:-all}" > "$log/cell.txt"
  [ -z "${cpus[$i]:-}" ] || affinity=(taskset -c "${cpus[$i]}")
  local args=(--out "$collection" --index "$i" --attempt "$attempt" --heap "$heap" --timeout "$timeout"
              --scratch "$cw" --container "$IMAGE")
  run "${wrap[@]}" --env "TMPDIR=$cw/tmp" "$IMAGE" "$PYTHON" \
    "$OUT_ROOT/code/stage2b/run_price_stage2b_collection.py" prepare-external "${args[@]}" \
    --launch-prefix "$work/prefix.json" > "$log/prepare.log" 2>&1 || prc=$?
  if [ "$prc" -ne 0 ]; then
    printf 'prepare %s\n' "$prc" > "$log/exit-code.txt"; date --iso-8601=seconds > "$log/finished-at.txt"; return 0
  fi
  locator=$OUT_ROOT/local-executions/$i-$attempt.json
  [ "$DRY_RUN" = "1" ] || [ "$(tail -n 1 "$log/prepare.log")" = "$locator" ] || { echo 'prepare 99' > "$log/exit-code.txt"; return 0; }
  run "${affinity[@]}" /usr/bin/time -v -o "$log/engine.resources.txt" /usr/bin/python3 -c \
    'import json,os,sys; c=json.load(open(sys.argv[1]))["command"]; os.execv(c[0],c)' "$locator" > "$log/engine.log" 2>&1 || rc=$?
  if [ "$DRY_RUN" = "1" ]; then sleep "0.$((RANDOM % 5))"; fi
  run "${wrap[@]}" --env "TMPDIR=$cw/tmp" "$IMAGE" "$PYTHON" \
    "$OUT_ROOT/code/stage2b/run_price_stage2b_collection.py" publish-external "${args[@]}" \
    --returncode "$rc" > "$log/publish.log" 2>&1 || prc=$?
  printf 'engine %s publish %s\n' "$rc" "$prc" > "$log/exit-code.txt"; date --iso-8601=seconds > "$log/finished-at.txt"
  if [ "$rc" -eq 0 ] && [ "$prc" -eq 0 ] && [ "$KEEP_WORK" != "1" ]; then rm -rf "$cw"; fi
  return 0
}

POLL=${POLL:-10}
declare -A running=()
used=0
reap() { local p; for p in "${!running[@]}"; do if ! kill -0 "$p" 2>/dev/null; then wait "$p" || true; used=$((used - running[$p])); unset "running[$p]"; fi; done; }
for i in "${cells[@]}"; do
  c=$(class_of "$i"); a=ALLOC_${c}_G; need=${!a}
  while [ "${#running[@]}" -gt 0 ] && [ $((used + need)) -gt "$MEM_BUDGET_G" ]; do sleep "$POLL"; reap; done
  run_cell "$i" & running[$!]=$need; used=$((used + need))
  echo "$(date --iso-8601=seconds) started cell $i (${need}G, in use ${used}G)"
done
while [ "${#running[@]}" -gt 0 ]; do sleep "$POLL"; reap; done

fails=0
for i in "${cells[@]}"; do
  status=$(cat "$jobout/cell-$i/exit-code.txt" 2>/dev/null || echo missing)
  echo "cell $i: $status"
  [ "$status" = 'engine 0 publish 0' ] || fails=$((fails + 1))
done
date --iso-8601=seconds > "$jobout/finished-at.txt"; printf '%s\n' "$fails" > "$jobout/failures.txt"
if [ "$fails" -eq 0 ] && [ "$KEEP_WORK" != "1" ]; then rm -rf "$work"; fi
echo "done, failures=$fails"; [ "$fails" -eq 0 ]
