#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=gspa-wp4-batches
#SBATCH --nodes=1
#SBATCH --time=08:00:00
#SBATCH --no-requeue
#SBATCH --output=/scratch/project/k10018/gspa-v3/out/wp4-full-v1-out/slurm/%j.log
# WP4 design §10 steps 3–5 on Shaheen III: gap records (gaps), fenced context inputs (context), the arm plans
# (plan) and the leakage audit (audit), run inside the unimatrix container recipe so that every path in the
# receipts and arm commands is a /data/hohndor/gspa-v3 path that the arm jobs (arm-pack.sbatch) resolve.
#
# Binds (nested bind targets must exist in the read-only mirror):
#   /project/k10018/gspa-v3                     -> /data/hohndor/gspa-v3 (ro)      the mirror
#   /scratch/.../out/wp4-search-v1              -> /data/hohndor/gspa-v3/wp4-search-v1 (ro)   search roots with observations
#   /scratch/.../out/wp4-full-v1-out            -> /data/hohndor/gspa-v3/wp4-full-v1-out      writable: batches/, jobs/, slurm/
#   package tree (relocated test env)           -> /home/leechuck                 Python 3.10 + harness deps
#   /project/k10018/code/gspa-v3                -> /home/leechuck/Public/software/gspa (ro)   the code slice at the pilot's repo path
# Arm commands written by `plan` map wp4-full-v1-out -> wp4-full-v1 (--path-map): after this job, batches/ is copied into
# the mirror bundle /project/k10018/gspa-v3/wp4-full-v1/batches and the observations into .../wp4-search-v1/<set>/<sk>/,
# so that arm-pack.sbatch (BUNDLE=wp4-full-v1) finds everything under the logical bundle paths.
#
# Environment: STAGES (default "tests gaps context plan audit"), GENOMES (targets; default: perturbations/removal/*),
# JAR_MAP=OLD_DIR=NEW_DIR adds a plan --path-map that rewrites the JAR in every arm command (2026-10-30: cohort
# permutation-seed fix, 10402 -> ec-prefix-identity-v4). KEEP_WORK=1 keeps the extracted package. Each stage refuses to overwrite its output (the preparers create their
# directories exclusively), so a rerun needs a fresh batches/ or explicit STAGES.
set -euo pipefail
module load singularity
image=/project/k10018/sw/images/debian-trixie
data_root=/project/k10018/gspa-v3; L=/data/hohndor/gspa-v3
search_out=/scratch/project/k10018/gspa-v3/out/wp4-search-v1
out_root=/scratch/project/k10018/gspa-v3/out/wp4-full-v1-out; root=$L/wp4-full-v1-out
code=/project/k10018/code/gspa-v3; repo=/home/leechuck/Public/software/gspa
work=/scratch/project/k10018/gspa-v3/work/$SLURM_JOB_ID
jobout=$out_root/jobs/$SLURM_JOB_ID
stages=${STAGES:-tests gaps context plan audit}
mkdir -p "$work/fs" "$jobout" "$out_root/slurm" "$out_root/batches"
export TMPDIR=$work/tmp; mkdir -p "$TMPDIR"
finish() {
  code=$?; trap - EXIT; set +e
  printf '%s\n' "$code" > "$jobout/exit-code.txt"; date --iso-8601=seconds > "$jobout/finished-at.txt"
  if [ "$code" -eq 0 ] && [ "${KEEP_WORK:-0}" != 1 ]; then rm -rf "$work"; fi
  exit "$code"
}
trap 'exit 143' TERM INT   # scancel: let the EXIT trap record 143, not the last command's status
trap finish EXIT
date --iso-8601=seconds > "$jobout/started-at.txt"; hostname > "$jobout/hostname.txt"; cp "$0" "$jobout/executed.sbatch"
zstd -dc "$data_root/observed-gap-context-modes-v1-relocation/package.tar.zst" | tar -xf - -C "$work/fs"
mkdir -p "$work/fs/home/leechuck/tmp" "$work/fs$repo"
container=(singularity exec --no-home --no-mount tmp --pwd "$repo"
 --bind "$data_root:$L:ro" --bind "$search_out:$L/wp4-search-v1:ro" --bind "$out_root:$root"
 --bind "$work/fs/home/leechuck:/home/leechuck" --bind "$code:$repo:ro"
 --env TMPDIR=/home/leechuck/tmp --env LC_ALL=C.UTF-8 --env PYTHONDONTWRITEBYTECODE=1
 --env XDG_CONFIG_HOME=/home/leechuck/.config --env XDG_CACHE_HOME=/home/leechuck/.cache --env XDG_DATA_HOME=/home/leechuck/.local/share
 --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env PYTHONPATH="$repo/benchmark/v3"
 "$image")
printf '%s\n' "${container[@]}" > "$jobout/container-command.txt"
python=/home/leechuck/tmp/gspa-v3-testenv/bin/python
src=$repo/benchmark/v3
B=$root/batches
S=$L/wp4-search-v1; PERT=$L/wp4-full-v1/perturbations; PEERS=$L/wp4-full-v1/peers; SEQ=$L/wp4-full-v1/sequence
EC2GO=$repo/benchmark/whole_genome_cafa/data/reference/ec2go.txt
genomes=${GENOMES:-$(ls "$data_root/wp4-full-v1/perturbations/removal")}
"${container[@]}" "$python" -c 'import sys; print(sys.version)' > "$jobout/python-version.txt" 2>&1
superkingdom() {  # superkingdom <set> <genome>: the (set, superkingdom) root whose manifest names the genome's variants
  local set=$1 g=$2 sk
  for sk in archaea bacteria; do
    grep -q "\"$g-" "$search_out/$set/$sk/manifest.json" && { echo "$sk"; return 0; }
  done
  echo "no $set search root names variants of $g" >&2; return 1
}
run() {  # run <stage> <args...>: one containerized stage with resources and a log
  local stage=$1; shift
  date --iso-8601=seconds > "$jobout/$stage.started"
  /usr/bin/time -v -o "$jobout/$stage.resources.txt" "${container[@]}" "$@" > "$jobout/$stage.log" 2>&1
  date --iso-8601=seconds > "$jobout/$stage.finished"
}
for stage in $stages; do
  case "$stage" in
  tests)
    run tests "$python" -m unittest discover -s "$src/tests" -p 'test_*wp4*.py' ;;
  gaps)
    pert=(); for g in $genomes; do pert+=(--perturbations "$PERT/removal/$g" --perturbations "$PERT/fragment/$g"); done
    run gaps "$python" "$src/prepare_wp4_batches.py" gaps --ec2go "$EC2GO" --out "$B/gaps" \
      --search-root "cohort-archaea=$S/cohort/archaea" --search-root "cohort-bacteria=$S/cohort/bacteria" \
      --search-root "peers-archaea=$S/peers/archaea" --search-root "peers-bacteria=$S/peers/bacteria" \
      --search-root "removal-archaea=$S/removal/archaea" --search-root "removal-bacteria=$S/removal/bacteria" \
      --search-root "fragment-archaea=$S/fragment/archaea" --search-root "fragment-bacteria=$S/fragment/bacteria" \
      --peers-dir "$PEERS" "${pert[@]}" \
      --cohort-root "$L/cohort-pathway-search-v1" --cohort-root "$L/cohort-pathway-search-archaea-v2" \
      --cohort-replay "$L/cohort-pathway-replay-v2" ;;
  context)
    # cohort + peers once (sequence metadata), then one restricted context per removal/fragment output (Task T).
    # Idempotent: a context whose receipt exists is kept (the preparer refuses to overwrite). A perturbation output
    # without copies (a target without peers has no removal copies) gets its context from the cohort + peer roots
    # alone: its metadata population is exactly the sequence population.
    if [ ! -f "$out_root/batches/context/wp4-context-receipt.json" ]; then
      run context "$python" "$src/prepare_wp4_batches.py" context --cohort-inputs "$L/cohort-context-modes-v2-inputs" \
        --cohort-root "$S/cohort/archaea" --cohort-root "$S/cohort/bacteria" \
        --search-root "peers-archaea=$S/peers/archaea" --search-root "peers-bacteria=$S/peers/bacteria" \
        --sequence "$SEQ" --out "$B/context"
    fi
    for g in $genomes; do
      for set in removal fragment; do
        [ -f "$out_root/batches/context-$set-$g/wp4-context-receipt.json" ] && continue
        copies=$(tail -n +2 "$data_root/wp4-full-v1/perturbations/$set/$g/gapsmith_inputs.tsv" 2>/dev/null | wc -l)
        extra=()
        if [ "$copies" -gt 0 ]; then
          sk=$(superkingdom "$set" "$g")
          extra=(--search-root "$set-$sk=$S/$set/$sk" --restrict-to-metadata)
        fi
        run "context-$set-$g" "$python" "$src/prepare_wp4_batches.py" context --cohort-inputs "$L/cohort-context-modes-v2-inputs" \
          --cohort-root "$S/cohort/archaea" --cohort-root "$S/cohort/bacteria" \
          --search-root "peers-archaea=$S/peers/archaea" --search-root "peers-bacteria=$S/peers/bacteria" \
          "${extra[@]}" --metadata "$PERT/$set/$g/cohort_metadata.json" \
          --sequence "$SEQ" --out "$B/context-$set-$g"
      done
    done ;;
  plan)
    pert=(--perturbations "$PERT/shuffle")
    for g in $genomes; do
      pert+=(--perturbations "$PERT/removal/$g=$B/context-removal-$g" --perturbations "$PERT/fragment/$g=$B/context-fragment-$g")
      [ -d "$data_root/wp4-full-v1/perturbations/duplicate/$g" ] && pert+=(--perturbations "$PERT/duplicate/$g")
    done
    run plan "$python" "$src/prepare_wp4_batches.py" plan --blocks all \
      --peers-tsv "$repo/benchmark/v3/results/peer-selection-v1/peers.tsv" \
      --plans-10406 "$L/role-genome-arm-v1/local/10406-plans" --peers-dir "$PEERS" --sequence "$SEQ" \
      --layouts "$L/training-cohort/layout" \
      --models-root "$L/source-context-full-v1/jobs/10406/checkpoint/training/diamond_homology" \
      --context "$B/context" --gaps "$B/gaps" "${pert[@]}" \
      --cohort-claims-dir "$L/source-context-full-v1/jobs/10406/checkpoint/training/claims/diamond_homology" \
      --path-map "$root=$L/wp4-full-v1" ${JAR_MAP:+--path-map "$JAR_MAP"} --out "$B/plans" ;;
  audit)
    # Leakage gate per superkingdom: every receipt written above, the fenced builder report + verify receipt, and all
    # proteins of the 18 targets (the 10406 cohort proteomes searched unchanged).
    # Globs expand on the host: list the receipt dirs under $out_root and name them by their container path.
    targets=(); for f in "$data_root"/cohort-pathway-search-v1/proteomes/*.faa; do   # all 18 targets (10406 cohort proteomes)
      targets+=("$L/${f#$data_root/}"); done
    receipts=("$B/gaps" "$B/plans"); for c in "$out_root"/batches/context*; do receipts+=("$B/$(basename "$c")"); done
    for sk in bacteria archaea; do
      run "audit-$sk" "$python" "$src/audit_wp4_leakage.py" --receipts "${receipts[@]}" \
        --fenced-report "$L/wp4-fenced-reference-v1/$sk/fenced_reference_report.json" \
        --fenced-verify "$L/wp4-fenced-reference-v1/$sk/verify/verify_report.json" \
        --targets "${targets[@]}" \
        --out "$B/leakage-$sk.json"
    done ;;
  seal|score|report)
    # Analysis (design §10 step 5) over the arm runs: out/wp4-full-v1 (written by arm-pack.sbatch, BUNDLE=wp4-full-v1)
    # is bound read-only over the mirror bundle so the runs resolve at the paths the arms tables name.
    acont=("${container[@]:0:${#container[@]}-1}" --bind /scratch/project/k10018/gspa-v3/out/wp4-full-v1:$L/wp4-full-v1:ro "$image")
    A=$B/analysis; mkdir -p "$out_root/batches/analysis"
    ain=$L/wp4-full-v1/analysis-inputs
    case "$stage" in
    seal)
      arms=(); for f in /scratch/project/k10018/gspa-v3/out/wp4-full-v1/batches/plans/*/arms.tsv; do
        arms+=("$L/wp4-full-v1/batches/plans/$(basename "$(dirname "$f")")/arms.tsv"); done
      date --iso-8601=seconds > "$jobout/seal.started"
      /usr/bin/time -v -o "$jobout/seal.resources.txt" "${acont[@]}" "$python" "$src/analyze_wp4_recurrence.py" seal \
        --arms "${arms[@]}" --runs "$L/wp4-full-v1/batches/plans" --out "$A/seal" --function-types EC,GO \
        --orthogroups "$SEQ/cohort_orthogroups.tsv" --peer-tiers "$ain/peer_tiers.tsv" \
        --claims-dir "$L/source-context-full-v1/jobs/10406/checkpoint/training/claims/diamond_homology" > "$jobout/seal.log" 2>&1
      date --iso-8601=seconds > "$jobout/seal.finished" ;;
    score)
      date --iso-8601=seconds > "$jobout/score.started"
      /usr/bin/time -v -o "$jobout/score.resources.txt" "${acont[@]}" "$python" "$src/analyze_wp4_recurrence.py" score \
        --labels "$L/training-cohort/labels/labels.tsv" --seal "$A/seal" --out "$A/scores" \
        --removals "$ain/removals.tsv" --go-obo "$L/training-cohort/go-basic.obo" > "$jobout/score.log" 2>&1
      date --iso-8601=seconds > "$jobout/score.finished" ;;
    report)
      date --iso-8601=seconds > "$jobout/report.started"
      "${acont[@]}" "$python" "$src/analyze_wp4_recurrence.py" report --scores "$A/scores" --out "$A/report" > "$jobout/report.log" 2>&1
      date --iso-8601=seconds > "$jobout/report.finished" ;;
    esac ;;
  *)
    echo "unknown stage: $stage" >&2; exit 2 ;;
  esac
done
echo done > "$jobout/done.txt"
