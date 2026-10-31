#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=gspa-wp4-fenced-ref
#SBATCH --nodes=1
#SBATCH --time=06:00:00
#SBATCH --no-requeue
#SBATCH --output=/scratch/project/k10018/gspa-v3/out/wp4-fenced-reference-v1/slurm/%j.log
# WP4 design §6: build the label-blind fenced gapsmith references (bacteria and archaea) and verify them
# (benchmark/v3/build_fenced_gapsmith_reference.py, Task P). Runs natively (conda python, no container):
# the builder is stdlib Python plus the DIAMOND identity screen (2.1.11, the claims' version).
set -euo pipefail
source /project/k10018/sw/eq-env.sh
PY=/project/k10018/sw/envs/gspa-v3/bin/python
DIAMOND=/project/k10018/sw/envs/diamond-2.1.11/bin/diamond
SRC=/project/k10018/code/gspa-v3            # worktree slice (read-only on compute nodes)
IN=/project/k10018/gspa-v3/wp4-fenced-inputs-v1
OUT=/scratch/project/k10018/gspa-v3/out/wp4-fenced-reference-v1
mkdir -p "$OUT/slurm" "$OUT/jobs/$SLURM_JOB_ID"; cd "$SRC"
export TMPDIR=/scratch/project/k10018/gspa-v3/work/$SLURM_JOB_ID/tmp; mkdir -p "$TMPDIR"
export PYTHONPATH=$SRC/benchmark/v3 OMP_NUM_THREADS=1
( cd "$IN" && sha256sum -c SHA256SUMS ) > "$OUT/jobs/$SLURM_JOB_ID/inputs-check.log"
"$DIAMOND" version > "$OUT/jobs/$SLURM_JOB_ID/diamond-version.txt"
run_kind() {  # $1 kind, $2 reference dir
  local kind=$1 ref=$2
  local out=$OUT/$kind log=$OUT/jobs/$SLURM_JOB_ID/$kind
  [ -e "$out" ] && { echo "$out exists; refusing to overwrite" >&2; return 2; }
  date --iso-8601=seconds > "$log.build.started"
  /usr/bin/time -v -o "$log.build.resources.txt" "$PY" benchmark/v3/build_fenced_gapsmith_reference.py build \
    --reference "$ref" --claim-fence "$IN/reference_report.json" --fence-provenance "$IN/provenance.json" \
    --targets "$IN/all_proteins.faa" --out "$out" --tool diamond --executable "$DIAMOND" --threads 90 > "$log.build.log" 2>&1
  date --iso-8601=seconds > "$log.build.finished"
  /usr/bin/time -v -o "$log.verify.resources.txt" "$PY" benchmark/v3/build_fenced_gapsmith_reference.py verify \
    --out "$out" --tool diamond --executable "$DIAMOND" --threads 90 --targets "$IN/all_proteins.faa" --check-inputs > "$log.verify.log" 2>&1
  date --iso-8601=seconds > "$log.verify.finished"
  echo "$kind done"
}
run_kind bacteria /project/k10018/gspa-v3/gapsmith-reference-bacteria-v1/reference &
run_kind archaea /project/k10018/gspa-v3/archaea-reference-v1/reference &
wait
rm -rf "/scratch/project/k10018/gspa-v3/work/$SLURM_JOB_ID"
echo done
