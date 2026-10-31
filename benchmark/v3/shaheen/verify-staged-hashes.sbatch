#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --exclude=nid00022
#SBATCH --time=02:00:00
#SBATCH --no-requeue
#SBATCH --job-name=gspa-verify-stage
#SBATCH --output=/scratch/project/k10018/gspa-v3/out/stage-verify/%j.log
# Check a staged tree on /project against a "sha256  relative/path" list made at the source, in parallel.
# Environment: ROOT (directory the paths are relative to), LIST (hash list), OUT (result file).
set -uo pipefail
: "${ROOT:?}" "${LIST:?}" "${OUT:?}"
cd "$ROOT"
split -n l/64 -d "$LIST" "$TMPDIR/part." 2>/dev/null || { mkdir -p /scratch/project/k10018/gspa-v3/work/$SLURM_JOB_ID; export TMPDIR=/scratch/project/k10018/gspa-v3/work/$SLURM_JOB_ID; split -n l/64 -d "$LIST" "$TMPDIR/part."; }
ls "$TMPDIR"/part.* | xargs -P 64 -I{} sh -c 'sha256sum -c {} 2>&1' > "$OUT.raw"
ok=$(grep -c ': OK$' "$OUT.raw"); bad=$(grep -v ': OK$' "$OUT.raw" | grep -c .)
{ echo "list=$LIST root=$ROOT lines=$(wc -l < "$LIST") ok=$ok not_ok=$bad"; grep -v ': OK$' "$OUT.raw" | head -50; } > "$OUT"
cat "$OUT"; [ "$bad" -eq 0 ]
