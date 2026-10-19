#!/bin/bash
# P02 -- Stage 1/2A ablation grid + nine-genome cross-clade perturbation (both designs).
# Host: unimatrix01 (SLURM, partition debug). Never run compute on the login node.
#
# Usage (from the unimatrix01 login node):
#   bash run.sh submit_idfix                 # CURRENT: v1+v2 with current jar on ID-normalised claims, + ablation ID check
#   bash run.sh submit                       # historical: three runs on raw claims IDs (v1/current, v2/current, v1/old jar)
#   bash run.sh collect                      # build the results tables from all five run summaries -> $REPRO/results/
#   sbatch ... run.sh perturb <design> <jar> <label> <out_dir> [normalised|raw]   # one run
#   sbatch ... run.sh ablation               # ablation grid (01 -> 02 -> 03 -> 05); NOT re-run for the paper
#
# Copy this folder to $PIPE (below) before submitting.
#SBATCH -p debug
#SBATCH --exclude=unimatrix,node007   # node007 kills our batch jobs at launch (signal 53; 5 occurrences on 2026-10-19)
#SBATCH -c 16
#SBATCH --mem=110G
#SBATCH -t 04:00:00

set -euo pipefail
PY=/storage/miniforge3/envs/metagenomics/bin/python3
REPRO=/data/hohndor/gspa/proteomes/bench9_unified/repro_2026-10-19
PIPE=$REPRO/P02_ablation_perturbation
JAR_CURRENT=/data/hohndor/gspa/bin/gspa-unified.jar                       # v2.1.0, Oct 17 12:48, constants 0.68/0.68
JAR_OLD=/data/hohndor/gspa/proteomes/bench9_unified/gspa-unified-bench.jar # Oct 14, constants 3.0/3.5 (no --batch-manifest)

mode=${1:-help}
case "$mode" in
  submit)
    mkdir -p "$REPRO/slurm"
    submit_one() {  # design jar label dir  (historical runs, raw claims IDs)
      sbatch -J "p02_$4" -o "$REPRO/slurm/p02_$4.out" -e "$REPRO/slurm/p02_$4.err" \
        "$PIPE/run.sh" perturb "$1" "$2" "$3" "$REPRO/$4" raw
    }
    [ "${2:-all}" = all ] && submit_one v1 "$JAR_CURRENT" gspa-unified.jar      v1_current_jar
    submit_one v2 "$JAR_CURRENT" gspa-unified.jar      v2_current_jar
    submit_one v1 "$JAR_OLD"     gspa-unified-bench.jar v1_old_jar
    ;;
  submit_idfix)
    # Current runs: claims ID-normalised + de-duplicated (run_perturbation.py default since 2026-10-19).
    # The three runs above predate the fix; `submit` passes claims_ids=raw to reproduce them.
    mkdir -p "$REPRO/slurm"
    for d in v1 v2; do
      sbatch --mem=60G -J "p02_${d}_idfix" -o "$REPRO/slurm/p02_${d}_current_jar_idfix.out" -e "$REPRO/slurm/p02_${d}_current_jar_idfix.err" \
        "$PIPE/run.sh" perturb "$d" "$JAR_CURRENT" gspa-unified.jar "$REPRO/${d}_current_jar_idfix" normalised
    done
    sbatch --mem=60G -J p02_abl_idfix -o "$REPRO/slurm/p02_ablation_idfix_check.out" -e "$REPRO/slurm/p02_ablation_idfix_check.err" \
      "$PIPE/run.sh" ablation_idfix_check
    ;;
  ablation_idfix_check)
    $PY "$PIPE/ablation/ablation_idfix_check.py" --out-dir "$REPRO/ablation_idfix_check" --jar "$JAR_CURRENT" --jobs 16
    ;;
  perturb)
    design=$2; jar=$3; label=$4; out=$5; claims_ids=${6:-normalised}
    echo "host=$(hostname) job=${SLURM_JOB_ID:-none} start=$(date -Is)"
    # per-job runner (old jar) uses 8g heaps -> 12 concurrent JVMs fit in 110G
    if unzip -p "$jar" gspa/cli/IntegrateCommand.class | grep -aq -- --batch-manifest; then J=16; else J=12; fi
    $PY "$PIPE/run_perturbation.py" --design "$design" --jar "$jar" --jar-label "$label" \
        --out-dir "$out" --jobs "$J" --claims-ids "$claims_ids"
    echo "end=$(date -Is)"
    ;;
  collect)  # writes pillar1_2_perturbation_{designs,summary}.tsv into $REPRO/results
    mkdir -p "$REPRO/results"
    $PY "$PIPE/collect_results.py" --results-dir "$REPRO/results" --runs \
      "$REPRO/v1_old_jar" "$REPRO/v1_current_jar" "$REPRO/v2_current_jar" \
      "$REPRO/v1_current_jar_idfix" "$REPRO/v2_current_jar_idfix"
    ;;
  ablation)
    # Ablation grid as originally run (job script bench9_unified/run_pillar1_2_slurm.sh, Oct 17).
    # The scripts hard-code WORK_DIR=/data/hohndor/gspa/proteomes/bench9_unified and JAR
    # /data/hohndor/gspa/bin/gspa-unified.jar; they overwrite integrated/, cafa_eval/, gaef_eval/.
    cd /data/hohndor/gspa/proteomes/bench9_unified
    $PY "$PIPE/ablation/01_run_ablation_grid.py" --jobs 16 --force
    $PY "$PIPE/ablation/02_run_cafa_eval.py" --jobs 16 --force
    $PY "$PIPE/ablation/03_run_gaef_eval.py" --jobs 16 --force
    $PY "$PIPE/ablation/05_compile_results.py" --out-dir "$REPRO/ablation_results"
    ;;
  *)
    sed -n '2,12p' "$0"; exit 1;;
esac
