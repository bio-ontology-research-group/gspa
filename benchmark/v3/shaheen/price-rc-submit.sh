#!/bin/bash
# Print (default) or submit (--submit) the Shaheen III schedule of the release-candidate Price collection:
#   1. price-rc-register.sbatch                         one node, fresh registration under the candidate JAR;
#   2. price-rc-cells.sbatch, afterok:register           all 132 cells packed onto nodes by memory class;
#   3. price-rc-cells.sbatch MODE=finalize, afterok:all  the unchanged collection finalizer.
# Run on a Shaheen login node from a checkout that holds this directory. The printed text is itself a shell
# script (job ids are captured into variables), so it can be reviewed and then run unchanged with --submit.
# Packing (original resource classes, all cells of a node start together, largest first):
#   large  101-131  LARGE_PER_NODE  (default 1: 168g heap, 192 GiB)  -> 31 nodes
#   medium  66-100  MEDIUM_PER_NODE (default 3: 96g heap, 112 GiB)   -> 12 nodes
#   small    0-65   SMALL_PER_NODE  (default 9: 24g heap, 36 GiB)    ->  8 nodes
# Environment overrides: OUT_ROOT JAR EXPECT_JAR_SHA256 EXPECT_WRAPPER_SHA256 FROZEN WRAPPER WRAPPER_TEST PYTHON PYLIB FSROOT IMAGE,
# *_PER_NODE, CELL_TIME / SMALL_TIME / FINAL_TIME, and anything price-rc-cells.sbatch reads (HEAP_*, ...).
set -euo pipefail
submit=0
case "${1:-}" in --submit) submit=1;; '') ;; *) echo "usage: $0 [--submit]" >&2; exit 2;; esac
here=$(cd "$(dirname "$0")" && pwd)
OUT_ROOT=${OUT_ROOT:-/scratch/project/k10018/gspa-v3/out/price-production-rc1}
JAR=${JAR:-/project/k10018/gspa-v3/release-candidate-1/gspa-3.0.0.jar}
EXPECT_JAR_SHA256=${EXPECT_JAR_SHA256:-}
LARGE_PER_NODE=${LARGE_PER_NODE:-1} MEDIUM_PER_NODE=${MEDIUM_PER_NODE:-3} SMALL_PER_NODE=${SMALL_PER_NODE:-9}
CELL_TIME=${CELL_TIME:-12:00:00} SMALL_TIME=${SMALL_TIME:-04:00:00} FINAL_TIME=${FINAL_TIME:-24:00:00}
common=(--account=k10018 --partition=workq --nodes=1 --exclude=nid00022 --no-requeue)
pass=(OUT_ROOT JAR EXPECT_JAR_SHA256 EXPECT_WRAPPER_SHA256 FROZEN WRAPPER WRAPPER_TEST PYTHON PYLIB FSROOT IMAGE WORK_ROOT MEM_BUDGET_G
      HEAP_SMALL HEAP_MEDIUM HEAP_LARGE ALLOC_SMALL_G ALLOC_MEDIUM_G ALLOC_LARGE_G TIMEOUT_SMALL TIMEOUT_MEDIUM
      TIMEOUT_LARGE PIN KEEP_WORK)
envs=()
for v in "${pass[@]}"; do [ -z "${!v:-}" ] || envs+=("$v=${!v}"); done

lines=()
emit() { lines+=("$1"); }
q() { [ "$#" -eq 0 ] || printf '%q ' "$@"; }
groups() {  # $1 = per node, then indices (largest first) -> one space-joined group per line
  local per=$1 n=0 g=(); shift
  for i in "$@"; do g+=("$i"); n=$((n+1)); if [ "$n" -eq "$per" ]; then echo "${g[*]}"; g=(); n=0; fi; done
  [ "${#g[@]}" -eq 0 ] || echo "${g[*]}"
}
emit "set -euo pipefail"
emit "mkdir -p $(q "$OUT_ROOT/slurm")"
emit "test ! -e $(q "$OUT_ROOT/collection")"
emit "REGISTER=\$(env $(q "${envs[@]}")sbatch --parsable $(q "${common[@]}") --time=02:00:00 --job-name=gspa-price-rc-register --output=$(q "$OUT_ROOT/slurm/register-%j.log") $(q "$here/price-rc-register.sbatch"))"
emit "echo register=\$REGISTER"
emit "CELL_JOBS=''"
k=0
while read -r cls per t; do
  case $cls in large) idx=$(seq 131 -1 101);; medium) idx=$(seq 100 -1 66);; small) idx=$(seq 65 -1 0);; esac
  # shellcheck disable=SC2086
  while read -r group; do
    k=$((k+1))
    emit "J=\$(env $(q "${envs[@]}")CELLS=$(q "$group") sbatch --parsable $(q "${common[@]}") --time=$t --dependency=afterok:\$REGISTER --job-name=gspa-price-rc-$cls --output=$(q "$OUT_ROOT/slurm/cells-%j.log") $(q "$here/price-rc-cells.sbatch")); CELL_JOBS=\$CELL_JOBS:\$J; echo \"$cls cells [$group] job=\$J\""
  done < <(groups "$per" $idx)
done <<EOF
large $LARGE_PER_NODE $CELL_TIME
medium $MEDIUM_PER_NODE $CELL_TIME
small $SMALL_PER_NODE $SMALL_TIME
EOF
emit "FINAL=\$(env $(q "${envs[@]}")MODE=finalize sbatch --parsable $(q "${common[@]}") --time=$FINAL_TIME --dependency=afterok\$CELL_JOBS --job-name=gspa-price-rc-finalize --output=$(q "$OUT_ROOT/slurm/finalize-%j.log") $(q "$here/price-rc-cells.sbatch"))"
emit "echo finalize=\$FINAL"
emit "printf '%s\\n' \"register \$REGISTER\" \"cells \${CELL_JOBS#:}\" \"finalize \$FINAL\" > \"$OUT_ROOT/slurm/schedule-\$REGISTER.txt\""

if [ "$submit" = 1 ]; then
  command -v sbatch > /dev/null || { echo 'sbatch not found: run on a Shaheen login node' >&2; exit 2; }
  printf '%s\n' "${lines[@]}" | bash
else
  echo "# $k cell jobs + register + finalize; review, then rerun with --submit"
  printf '%s\n' "${lines[@]}"
fi
