#!/bin/bash
# Off-growth enzyme discovery, second round (plan tasks/v3-row12-offgrowth-improvement-2026-11-02.md):
# one exclusive node per (genome, candidate-table kind). Runs the masked test superset for every case of the
# genome and the unmasked deployment run, packs <genome>.<kind>.outputs.tar on /scratch.
#
# Usage from a Shaheen login node (dry run prints the sbatch commands; --submit submits):
#   TIER=A KIND=six  GENOMES="GCA_000006965.1 GCF_001633105.1" bash benchmark/v3/shaheen/offgrowth2-genome.sbatch --submit
#   TIER=A KIND=clean ... ; TIER=B KIND=clean ... ; TIER=C KIND=clean ...
# Environment (defaults below): REPO_ROOT, OUT_ROOT, MODELS_TAR, CANDIDATES_DIR, REF_DIR, PANEL_TSV, KEGG_MODULES,
# LAYOUT_DIR, LABEL_CASES_TSV, EXTRA_CASES_TSV, PRICE_MANIFEST, PYTHON_BIN, WALLTIME, MAX_CASES, CASE_IDS, SKIP_UNMASKED.
set -euo pipefail

REPO_ROOT=${REPO_ROOT:-/lustre2/project/k10018/code/gspa-v3}
OUT_ROOT=${OUT_ROOT:-/scratch/project/k10018/gspa-v3/out/offgrowth2-v1}
INPUTS=${INPUTS:-/project/k10018/gspa-v3/offgrowth2-inputs-v1}
TIER=${TIER:-A}
KIND=${KIND:-six}
MODELS_TAR=${MODELS_TAR:-/scratch/project/k10018/gspa-v3/benchmark-models-v1/batches/batch_000/outputs.tar}
CANDIDATES_DIR=${CANDIDATES_DIR:-$INPUTS/candidates-$KIND}
REF_DIR=${REF_DIR:-/project/k10018/gspa-v3/gapsmith-reference-bacteria-v1/reference}
PANEL_TSV=${PANEL_TSV:-$REPO_ROOT/benchmark/v3/offgrowth/data/biolog_pm1_pm4_panel.tsv}
KEGG_MODULES=${KEGG_MODULES:-$INPUTS/kegg_modules.tsv}
LAYOUT_DIR=${LAYOUT_DIR:-$INPUTS/layout}
PYTHON_BIN=${PYTHON_BIN:-/project/k10018/sw/envs/gspa-v3/bin/python3}
WALLTIME=${WALLTIME:-08:00:00}
WORKERS=${WORKERS:-192}

# Case and label sources per tier (plan §3, §5). Tier A: Price-148 + Swiss-Prot (L0, re-run) + Fitness Browser (L1).
# Tier B: Swiss-Prot experimental ECs of the cohort bacteria (L2). Tier C: Fitness Browser reannotations (L3).
case "$TIER" in
    A) PRICE_MANIFEST=${PRICE_MANIFEST:-$REPO_ROOT/benchmark/manifests/suite6_price149_manifest.tsv}
       EXTRA_CASES_TSV=${EXTRA_CASES_TSV:-$REPO_ROOT/benchmark/v3/offgrowth/data/swissprot_experimental_ecs_price22.tsv}
       LABEL_CASES_TSV=${LABEL_CASES_TSV:-$REPO_ROOT/benchmark/v3/offgrowth2/data/fitness_browser_ec_reannotations_price22.tsv} ;;
    B) PRICE_MANIFEST=${PRICE_MANIFEST:-}
       EXTRA_CASES_TSV=${EXTRA_CASES_TSV:-$REPO_ROOT/benchmark/v3/offgrowth2/data/swissprot_experimental_ecs_cohort8.tsv}
       LABEL_CASES_TSV=${LABEL_CASES_TSV:-} ;;
    C) PRICE_MANIFEST=${PRICE_MANIFEST:-}
       EXTRA_CASES_TSV=${EXTRA_CASES_TSV:-}
       LABEL_CASES_TSV=${LABEL_CASES_TSV:-$REPO_ROOT/benchmark/v3/offgrowth2/data/fitness_browser_ec_reannotations_tierc.tsv}
       MODELS_TAR=${MODELS_TAR_C:-/scratch/project/k10018/gspa-v3/tierc-models-v1/batches/batch_000/outputs.tar} ;;
    *) echo "TIER must be A, B or C" >&2; exit 2 ;;
esac

script_path="$(readlink -f "${BASH_SOURCE[0]}")"
if [ -z "${SLURM_JOB_ID:-}" ]; then
    submit=0
    [ "${1:-}" = "--submit" ] && submit=1
    : "${GENOMES:?GENOMES must list the genome ids to run}"
    mkdir -p "$OUT_ROOT/slurm" "$OUT_ROOT/genomes" "$OUT_ROOT/failed"
    for gid in $GENOMES; do
        job_name="gspa-og2-${TIER}-${KIND}-${gid}"
        log_file="$OUT_ROOT/slurm/${job_name}-%j.log"
        export_vars="ALL,TIER=$TIER,KIND=$KIND,GENOME_ID=$gid,REPO_ROOT=$REPO_ROOT,OUT_ROOT=$OUT_ROOT,INPUTS=$INPUTS,MODELS_TAR=$MODELS_TAR,CANDIDATES_DIR=$CANDIDATES_DIR,REF_DIR=$REF_DIR,PANEL_TSV=$PANEL_TSV,KEGG_MODULES=$KEGG_MODULES,LAYOUT_DIR=$LAYOUT_DIR,PYTHON_BIN=$PYTHON_BIN,WORKERS=$WORKERS,PRICE_MANIFEST=$PRICE_MANIFEST,EXTRA_CASES_TSV=$EXTRA_CASES_TSV,LABEL_CASES_TSV=$LABEL_CASES_TSV,MAX_CASES=${MAX_CASES:-},CASE_IDS=${CASE_IDS:-},SKIP_UNMASKED=${SKIP_UNMASKED:-}"
        cmd="sbatch --parsable --account=k10018 --partition=workq --nodes=1 --ntasks=1 --cpus-per-task=192 --exclusive --no-requeue --exclude=nid00022,nid00497,nid01033 --time=$WALLTIME --job-name=$job_name --output=$log_file --export=$export_vars $script_path"
        echo "$cmd"
        if [ "$submit" -eq 1 ]; then
            jid=$(eval "$cmd")
            echo "submitted $job_name as $jid"
        fi
    done
    exit 0
fi

: "${GENOME_ID:?GENOME_ID must be set}"
WORK_DIR="/dev/shm/offgrowth2-${SLURM_JOB_ID}-${GENOME_ID}"
STAGE_OUT="$WORK_DIR/out/${GENOME_ID}"
mkdir -p "$WORK_DIR/models" "$WORK_DIR/candidates" "$WORK_DIR/ref" "$STAGE_OUT" "$OUT_ROOT/genomes" "$OUT_ROOT/slurm"

cleanup() {
    local code=$?
    trap - EXIT
    set +e
    if [ "$code" -ne 0 ] && [ -d "$STAGE_OUT" ]; then
        mkdir -p "$OUT_ROOT/failed/${SLURM_JOB_ID}-${GENOME_ID}-${KIND}"
        find "$STAGE_OUT" -maxdepth 1 -type f \( -name '*.txt' -o -name '*.log' -o -name '*.sha256' -o -name '*.json' -o -name '*.tsv' \) \
            -exec cp {} "$OUT_ROOT/failed/${SLURM_JOB_ID}-${GENOME_ID}-${KIND}/" \;
        echo "exit $code; logs kept in $OUT_ROOT/failed/${SLURM_JOB_ID}-${GENOME_ID}-${KIND}" >&2
    fi
    rm -rf "$WORK_DIR"
    exit "$code"
}
trap 'exit 143' TERM INT
trap cleanup EXIT

date --iso-8601=seconds > "$STAGE_OUT/started_at.txt"
hostname > "$STAGE_OUT/hostname.txt"
cp "$0" "$STAGE_OUT/executed.sbatch"

# 1. Reference tables, model members, candidate table and case files into node memory (one extraction per node)
cp -f "$PANEL_TSV" "$WORK_DIR/ref/biolog_pm1_pm4_panel.tsv"
cp -f "$REF_DIR/seed_metabolites_edited.tsv" "$REF_DIR/seed_reactions_corrected.tsv" \
      "$REF_DIR/seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv" "$WORK_DIR/ref/"
cp -f "$KEGG_MODULES" "$WORK_DIR/ref/kegg_modules.tsv"
tar -xf "$MODELS_TAR" -C "$WORK_DIR/models" "${GENOME_ID}/"
cp -f "$CANDIDATES_DIR/${GENOME_ID}.candidates.tsv" "$WORK_DIR/candidates/${GENOME_ID}.candidates.tsv"
CASE_ARGS=()
HASH_FILES=()
if [ -n "${PRICE_MANIFEST:-}" ]; then
    cp -f "$PRICE_MANIFEST" "$WORK_DIR/ref/price_manifest.tsv"; CASE_ARGS+=(--price-manifest "$WORK_DIR/ref/price_manifest.tsv"); HASH_FILES+=("$WORK_DIR/ref/price_manifest.tsv")
fi
if [ -n "${EXTRA_CASES_TSV:-}" ]; then
    cp -f "$EXTRA_CASES_TSV" "$WORK_DIR/ref/extra_cases.tsv"; CASE_ARGS+=(--extra-cases-tsv "$WORK_DIR/ref/extra_cases.tsv"); HASH_FILES+=("$WORK_DIR/ref/extra_cases.tsv")
fi
if [ -n "${LABEL_CASES_TSV:-}" ]; then
    cp -f "$LABEL_CASES_TSV" "$WORK_DIR/ref/label_cases.tsv"; CASE_ARGS+=(--label-cases-tsv "$WORK_DIR/ref/label_cases.tsv" --label-cases-prefix "FB:"); HASH_FILES+=("$WORK_DIR/ref/label_cases.tsv")
fi
LAYOUT_ARGS=()
if [ -f "$LAYOUT_DIR/${GENOME_ID}.layout.tsv" ]; then
    cp -f "$LAYOUT_DIR/${GENOME_ID}.layout.tsv" "$WORK_DIR/ref/layout.tsv"; LAYOUT_ARGS=(--layout-tsv "$WORK_DIR/ref/layout.tsv"); HASH_FILES+=("$WORK_DIR/ref/layout.tsv")
fi
EXTRA_ARGS=(--workers "$WORKERS" --candidate-table-kind "$KIND")
[ -n "${MAX_CASES:-}" ] && EXTRA_ARGS+=(--max-cases "$MAX_CASES")
[ -n "${CASE_IDS:-}" ] && EXTRA_ARGS+=(--case-ids "$CASE_IDS")
[ -n "${SKIP_UNMASKED:-}" ] && EXTRA_ARGS+=(--skip-unmasked)

sha256sum "$WORK_DIR"/ref/*.tsv "$WORK_DIR/models/${GENOME_ID}/${GENOME_ID}-draft.xml" \
    "$WORK_DIR/models/${GENOME_ID}/${GENOME_ID}-filled.xml" "$WORK_DIR/models/${GENOME_ID}/${GENOME_ID}-medium.csv" \
    "$WORK_DIR/models/${GENOME_ID}/${GENOME_ID}-Transporter.tbl" "$WORK_DIR/candidates/${GENOME_ID}.candidates.tsv" \
    > "$STAGE_OUT/inputs.sha256"
(cd "$REPO_ROOT/benchmark/v3" && find offgrowth2 offgrowth followups -name '*.py' -print0 | sort -z | xargs -0 sha256sum) > "$STAGE_OUT/code.sha256"

cd "$REPO_ROOT/benchmark/v3"
export PYTHONPATH="$REPO_ROOT/benchmark/v3"
"$PYTHON_BIN" -c 'import sys, cobra, optlang; from offgrowth.arms import detect_solvers; print("python:", sys.version); print("cobra:", cobra.__version__); print("optlang:", optlang.__version__); print("solvers:", detect_solvers())' > "$STAGE_OUT/versions.txt" 2>&1

# 2. Run the superset (masked cases + unmasked run) under /usr/bin/time -v
/usr/bin/time -v -o "$STAGE_OUT/time_v.txt" \
    "$PYTHON_BIN" -m offgrowth2.runner genome \
    --genome-id "$GENOME_ID" \
    --genome-model-dir "$WORK_DIR/models/${GENOME_ID}" \
    --candidates-tsv "$WORK_DIR/candidates/${GENOME_ID}.candidates.tsv" \
    --seed-metabolites "$WORK_DIR/ref/seed_metabolites_edited.tsv" \
    --seed-reactions "$WORK_DIR/ref/seed_reactions_corrected.tsv" \
    --ec-aliases "$WORK_DIR/ref/seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv" \
    --panel-tsv "$WORK_DIR/ref/biolog_pm1_pm4_panel.tsv" \
    --kegg-modules "$WORK_DIR/ref/kegg_modules.tsv" \
    "${CASE_ARGS[@]}" "${LAYOUT_ARGS[@]}" "${EXTRA_ARGS[@]}" \
    --out-dir "$WORK_DIR/out" \
    2>&1 | tee "$STAGE_OUT/runner.log" > "$OUT_ROOT/slurm/runner-${SLURM_JOB_ID}-${GENOME_ID}-${KIND}.log"

date --iso-8601=seconds > "$STAGE_OUT/finished_at.txt"
TAR_DEST="$OUT_ROOT/genomes/${GENOME_ID}.${KIND}.outputs.tar"
tar -cf "${TAR_DEST}.partial" -C "$WORK_DIR/out" "${GENOME_ID}"
mv "${TAR_DEST}.partial" "$TAR_DEST"
sha256sum "$TAR_DEST" > "${TAR_DEST}.sha256"
