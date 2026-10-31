#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=gspa-price-rc-cells
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --no-requeue
#SBATCH --exclude=nid00022
#SBATCH --output=/scratch/project/k10018/gspa-v3/out/price-production-rc1/slurm/cells-%j.log
# Shaheen III execution of registered Price cells (price-rc-register.sbatch) on one exclusive node, or the final
# collection audit. Shaheen equivalent of unimatrix price-production-v3/{small,medium,large}.sbatch and
# price-finalize-v2/finalize.sbatch, with run_price_local.py semantics unchanged: private output on job scratch bound
# at the logical attempt path, mount guard, both arms, original receipts, byte-checked atomic publication.
# The wrapper runs inside the container, where Singularity cannot be nested, so it uses the opt-in external launch:
#   1. container: run_price_local.py --phase prepare  checks and reserves the cell, records the engine command;
#   2. host:      executes exactly that recorded argv (os.execv of the JSON list), keeps its return code;
#   3. container: run_price_local.py --phase publish --returncode N  checks runtime identity and publishes.
# MODE=cells (default): CELLS = cell indices (space, comma or + separated). Cells start in the given order whenever
#   their memory class fits MEM_BUDGET_G; when all fit at once, each gets a disjoint block of physical cores.
#   Classes are the original resource waves: small 0-65 (heap 24g, 36 GiB, 7200 s/arm), medium 66-100 (96g,
#   112 GiB, 14400 s), large 101-131 (168g, 192 GiB, 21600 s). Overrides: HEAP_*, ALLOC_*_G, TIMEOUT_*.
# MODE=finalize: run_price_production_collection.py finalize (unchanged original finalizer) after a preflight.
# DRY_RUN=1 prints every command and simulates the cells without containers (for review; never on login nodes).
set -euo pipefail
MODE=${MODE:-cells}
OUT_ROOT=${OUT_ROOT:-/scratch/project/k10018/gspa-v3/out/price-production-rc1}
PYTHON=${PYTHON:-/project/k10018/sw/envs/gspa-tools/bin/python3.12}
PYLIB=${PYLIB:-/data/hohndor/gspa-v3/price-production-v1/venv/lib/python3.12/site-packages}
FSROOT=${FSROOT:-/scratch/project/k10018/gspa-v3/shared/arm-fs-d04c9e11ee64db47}
IMAGE=${IMAGE:-/project/k10018/sw/images/debian-trixie}
WORK_ROOT=${WORK_ROOT:-/scratch/project/k10018/gspa-v3/work}
MEM_BUDGET_G=${MEM_BUDGET_G:-340}     # of 370,000 MB (361 GiB) per workq node
HEAP_SMALL=${HEAP_SMALL:-24g} HEAP_MEDIUM=${HEAP_MEDIUM:-96g} HEAP_LARGE=${HEAP_LARGE:-168g}
ALLOC_SMALL_G=${ALLOC_SMALL_G:-36} ALLOC_MEDIUM_G=${ALLOC_MEDIUM_G:-112} ALLOC_LARGE_G=${ALLOC_LARGE_G:-192}
TIMEOUT_SMALL=${TIMEOUT_SMALL:-7200} TIMEOUT_MEDIUM=${TIMEOUT_MEDIUM:-14400} TIMEOUT_LARGE=${TIMEOUT_LARGE:-21600}
PIN=${PIN:-1} DRY_RUN=${DRY_RUN:-0} KEEP_WORK=${KEEP_WORK:-0}
JOB=${SLURM_JOB_ID:-dry$$}
logical_data=/data/hohndor/gspa-v3
collection=$OUT_ROOT/collection
jobout=$OUT_ROOT/jobs/$JOB
work=$WORK_ROOT/$JOB
unset JAVA_TOOL_OPTIONS JDK_JAVA_OPTIONS _JAVA_OPTIONS PYTHONPATH PYTHONHOME

if [ "$DRY_RUN" = 1 ]; then
  launcher=/usr/bin/singularity
  jobout=${DRY_OUT:-${TMPDIR:?}/price-rc-dry/$JOB}; work=$jobout/work
else
  module load singularity
  launcher=$(readlink -f "$(command -v singularity)")
  [ -f "$collection/registration.json" ] && [ -f "$OUT_ROOT/register-receipt.json" ] || { echo "collection not registered: $collection" >&2; exit 2; }
  [ -f "$FSROOT/.complete" ] || { echo "shared package tree missing: $FSROOT" >&2; exit 2; }
fi
mkdir -p "$jobout" "$work/tmp"
date --iso-8601=seconds > "$jobout/started-at.txt"; hostname > "$jobout/hostname.txt"
cp "$0" "$jobout/executed.sbatch"

# Host launcher argv up to (excluding) the image. It is recorded verbatim in every execution locator.
prefix=("$launcher" exec --no-home --no-mount tmp
  --bind /project/k10018/sw/envs/gspa-tools:/project/k10018/sw/envs/gspa-tools:ro
  --bind "/project/k10018/gspa-v3:$logical_data:ro" --bind "$OUT_ROOT:$OUT_ROOT" --bind "$work:$work"
  --bind "$FSROOT/usr/lib/jvm/java-21-openjdk-amd64:/usr/lib/jvm/java-21-openjdk-amd64:ro" --bind "$FSROOT/home/leechuck:/home/leechuck:ro"
  --bind /usr/bin/time:/usr/bin/time:ro
  --env PATH=/usr/lib/jvm/java-21-openjdk-amd64/bin:/usr/bin:/bin --env LC_ALL=C.UTF-8 --env PYTHONDONTWRITEBYTECODE=1
  --env "PYTHONPATH=$PYLIB")
# The wrapper additionally sees the launcher and the image it hashes (runtime_files, container_files).
wrap=("${prefix[@]}" --bind "$launcher:$launcher:ro" --bind "$IMAGE:$IMAGE:ro")
/usr/bin/python3 -c 'import json,sys; print(json.dumps(sys.argv[1:]))' "${prefix[@]}" > "$work/prefix.json"
cp "$work/prefix.json" "$jobout/prefix.json"
run() { if [ "$DRY_RUN" = 1 ]; then printf '%q ' "$@" >> "$jobout/dry-run.txt"; echo >> "$jobout/dry-run.txt"; else "$@"; fi; }

if [ "$MODE" = finalize ]; then
  export OUT_ROOT
  run "${prefix[@]}" --env "TMPDIR=$work/tmp" "$IMAGE" "$PYTHON" - > "$jobout/preflight.log" 2>&1 <<'PY'
import hashlib, json, os
from pathlib import Path
out = Path(os.environ['OUT_ROOT']); sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
receipt = json.loads((out/'register-receipt.json').read_text())
registration = json.loads((out/'collection/registration.json').read_text())
assert sha(out/'collection/registration.json') == receipt['registration_sha256']
assert sha(out/'collection/plan.json') == registration['plan_sha256'] == receipt['plan_sha256']
files = set((out/'code').glob('*.py'))
expected = {Path(k) for k in registration['inputs'] if Path(k).parent == out/'code' and Path(k).suffix == '.py'}
assert expected and files == expected
for p in files: assert sha(p) == registration['inputs'][str(p)], p
assert not (out/'collection/collection-verification.json').exists()
print('registration, plan and all', len(files), 'registered Python sources verified; jar', receipt['jar_sha256'])
PY
  code=0
  run /usr/bin/time -v -o "$jobout/resources.txt" "${prefix[@]}" --env "TMPDIR=$work/tmp" "$IMAGE" "$PYTHON" \
    "$OUT_ROOT/code/run_price_production_collection.py" finalize --out "$collection" > "$jobout/finalize.log" 2>&1 || code=$?
  printf '%s\n' "$code" > "$jobout/exit-code.txt"; date --iso-8601=seconds > "$jobout/finished-at.txt"
  [ "$code" -eq 0 ] || exit "$code"
  [ "$DRY_RUN" = 1 ] && exit 0
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
# Disjoint physical-core blocks (with their hyperthread siblings) when every cell of this node runs at once.
declare -A cpus=()
if [ "$PIN" = 1 ] && [ "$total" -le "$MEM_BUDGET_G" ] && [ "${#cells[@]}" -gt 1 ]; then
  mapfile -t cores < <(lscpu -p=CORE | grep -v '^#' | sort -n -u)
  per=$(( ${#cores[@]} / ${#cells[@]} ))
  for k in "${!cells[@]}"; do
    lo=${cores[$((k*per))]}; hi=${cores[$(((k+1)*per-1))]}
    cpus[${cells[$k]}]=$(lscpu -p=CPU,CORE | grep -v '^#' | awk -F, -v lo="$lo" -v hi="$hi" '$2>=lo && $2<=hi {print $1}' | paste -sd,)
  done
fi
echo "job $JOB on $(hostname): cells=[${cells[*]}] alloc=${total}G budget=${MEM_BUDGET_G}G pin=${#cpus[@]}"

run_cell() {   # $1 = cell index; never returns nonzero (the summary below decides)
  local i=$1 c heap timeout attempt cw log locator rc=0 prc=0 affinity=()
  c=$(class_of "$i"); heap=HEAP_$c; timeout=TIMEOUT_$c; heap=${!heap}; timeout=${!timeout}
  attempt=job-$JOB; cw=$work/cell-$i; log=$jobout/cell-$i
  mkdir -p "$cw/tmp" "$log"
  date --iso-8601=seconds > "$log/started-at.txt"
  echo "cell $i class=$c heap=$heap timeout=$timeout cpus=${cpus[$i]:-all}" > "$log/cell.txt"
  [ -z "${cpus[$i]:-}" ] || affinity=(taskset -c "${cpus[$i]}")
  local args=(--out "$collection" --index "$i" --attempt "$attempt" --heap "$heap" --timeout "$timeout"
              --scratch "$cw" --container "$IMAGE")
  run "${wrap[@]}" --env "TMPDIR=$cw/tmp" "$IMAGE" "$PYTHON" "$OUT_ROOT/code/run_price_local.py" "${args[@]}" \
    --phase prepare --launch-prefix "$work/prefix.json" > "$log/prepare.log" 2>&1 || prc=$?
  if [ "$prc" -ne 0 ]; then
    printf 'prepare %s\n' "$prc" > "$log/exit-code.txt"; date --iso-8601=seconds > "$log/finished-at.txt"; return 0
  fi
  locator=$OUT_ROOT/local-executions/$i-$attempt.json
  [ "$DRY_RUN" = 1 ] || [ "$(tail -n 1 "$log/prepare.log")" = "$locator" ] || { echo 'prepare 99' > "$log/exit-code.txt"; return 0; }
  # Host step: the recorded argv, verbatim (it contains a multi-line -c program, so no shell word splitting).
  run "${affinity[@]}" /usr/bin/time -v -o "$log/engine.resources.txt" /usr/bin/python3 -c \
    'import json,os,sys; c=json.load(open(sys.argv[1]))["command"]; os.execv(c[0],c)' "$locator" > "$log/engine.log" 2>&1 || rc=$?
  if [ "$DRY_RUN" = 1 ]; then sleep "0.$((RANDOM % 5))"; fi
  run "${wrap[@]}" --env "TMPDIR=$cw/tmp" "$IMAGE" "$PYTHON" "$OUT_ROOT/code/run_price_local.py" "${args[@]}" \
    --phase publish --returncode "$rc" > "$log/publish.log" 2>&1 || prc=$?
  printf 'engine %s publish %s\n' "$rc" "$prc" > "$log/exit-code.txt"; date --iso-8601=seconds > "$log/finished-at.txt"
  if [ "$rc" -eq 0 ] && [ "$prc" -eq 0 ] && [ "$KEEP_WORK" != 1 ]; then rm -rf "$cw"; fi
  return 0
}

# Memory-class admission: start the next cell when it fits beside the running ones (a lone cell always starts).
# Polling, not wait -n: bash reaps finished children itself, so kill -0 fails once a cell has exited.
POLL=${POLL:-10}
declare -A running=()   # pid -> GiB
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
if [ "$fails" -eq 0 ] && [ "$KEEP_WORK" != 1 ]; then rm -rf "$work"; fi
echo "done, failures=$fails"; [ "$fails" -eq 0 ]
