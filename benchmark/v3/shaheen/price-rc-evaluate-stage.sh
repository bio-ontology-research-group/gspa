#!/bin/bash
# Stage evaluation inputs for the release-candidate Price collection on Shaheen III.
# Copies archive inputs and source repository files into EVAL_ROOT, then verifies
# them using sha256sum -c against package.sha256 and source.sha256.
#
# Defaults:
#   EVAL_ROOT = /scratch/project/k10018/gspa-v3/out/price-reference-rc1
#
# Modes:
#   ./price-rc-evaluate-stage.sh        # dry-run: prints commands only
#   ./price-rc-evaluate-stage.sh --run  # executes staging and verification
set -euo pipefail

run_mode=0
case "${1:-}" in
  --run) run_mode=1;;
  '') run_mode=0;;
  *) echo "usage: $0 [--run]" >&2; exit 2;;
esac

here=$(cd "$(dirname "$0")" && pwd)
REPO_ROOT=${REPO_ROOT:-$(cd "$here/../../.." && pwd)}
ARCHIVE_DIR=${ARCHIVE_DIR:-$REPO_ROOT/benchmark/v3/results/price-reference-full-launch-v1}
EVAL_ROOT=${EVAL_ROOT:-/scratch/project/k10018/gspa-v3/out/price-reference-rc1}

[ -d "$REPO_ROOT" ] || { echo "REPO_ROOT not found: $REPO_ROOT" >&2; exit 2; }
[ -d "$ARCHIVE_DIR" ] || { echo "ARCHIVE_DIR not found: $ARCHIVE_DIR" >&2; exit 2; }
[ -f "$ARCHIVE_DIR/source.sha256" ] || { echo "source.sha256 not found: $ARCHIVE_DIR/source.sha256" >&2; exit 2; }
[ -f "$ARCHIVE_DIR/package.sha256" ] || { echo "package.sha256 not found: $ARCHIVE_DIR/package.sha256" >&2; exit 2; }

cmd() {
  if [ "$run_mode" = 1 ]; then
    "$@"
  else
    printf '%q ' "$@"
    echo
  fi
}

# 1. Ensure target directory structure exists
cmd mkdir -p "$EVAL_ROOT" "$EVAL_ROOT/slurm" "$EVAL_ROOT/jobs" "$EVAL_ROOT/source"

# 2. Copy evaluation archive files
ARCHIVE_FILES=(
  package.sha256
  source.sha256
  archive.sha256
  check_retained.py
  ec-audit.json
  metadata_preflight.py
  prerequisites.py
  run.sbatch
  README.md
  metadata-preflight.json
  review-prompt.txt
  review.err
  review.json
  scheduler.txt
  tests.log
)

for f in "${ARCHIVE_FILES[@]}"; do
  if [ -f "$ARCHIVE_DIR/$f" ]; then
    cmd cp "$ARCHIVE_DIR/$f" "$EVAL_ROOT/$f"
  fi
done

# 3. Copy Shaheen evaluate companion files
if [ -f "$here/prerequisites_shaheen.py" ]; then
  cmd cp "$here/prerequisites_shaheen.py" "$EVAL_ROOT/prerequisites_shaheen.py"
fi
if [ -f "$here/price-rc-evaluate.sbatch" ]; then
  cmd cp "$here/price-rc-evaluate.sbatch" "$EVAL_ROOT/price-rc-evaluate.sbatch"
fi

# 4. Copy repo source files listed in source.sha256
while read -r digest rel; do
  [ -n "$rel" ] || continue
  src="$REPO_ROOT/$rel"
  dst="$EVAL_ROOT/source/$rel"
  cmd mkdir -p "$(dirname "$dst")"
  cmd cp "$src" "$dst"
done < "$ARCHIVE_DIR/source.sha256"

# 5. Verification
if [ "$run_mode" = 1 ]; then
  echo "Verifying staged archive package files in $EVAL_ROOT..."
  ( cd "$EVAL_ROOT" && sha256sum -c package.sha256 )
  echo "Verifying staged source files in $EVAL_ROOT/source..."
  ( cd "$EVAL_ROOT/source" && sha256sum -c "$EVAL_ROOT/source.sha256" )
  echo "Staging complete and verified: $EVAL_ROOT"
else
  printf 'cd %q && sha256sum -c package.sha256\n' "$EVAL_ROOT"
  printf 'cd %q && sha256sum -c %q\n' "$EVAL_ROOT/source" "$EVAL_ROOT/source.sha256"
fi
