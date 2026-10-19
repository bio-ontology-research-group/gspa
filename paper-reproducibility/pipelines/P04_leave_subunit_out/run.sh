#!/bin/bash
# P04 -- Leave-subunit-out (LSO) heteromeric complex benchmark, two decoy designs.
# Host: unimatrix01 (SLURM, partition debug). Never run compute on the login node.
#
# Usage
#   (local checkout)      bash run.sh deploy     # rsync entry points to $PIPE on unimatrix01
#   (unimatrix01 login)   bash run.sh submit     # sbatch the full pipeline (prepare -> integrate -> agent -> evaluate)
#   (local checkout)      bash run.sh fetch      # copy result TSVs back into results/
#   sbatch ... run.sh job [step]                 # what `submit` runs on the compute node (step default: all)
#
#SBATCH -p debug
#SBATCH --exclude=unimatrix,node007   # node007 fails jobs at launch (signal 53) on 2026-10-19
#SBATCH -c 16
#SBATCH --mem=110G
#SBATCH -t 08:00:00

set -euo pipefail
PY=/storage/miniforge3/envs/metagenomics/bin/python3
REPRO=/data/hohndor/gspa/proteomes/pillar3_unified/repro_2026-10-19
PIPE=$REPRO/P04_leave_subunit_out
JAR=/data/hohndor/gspa/proteomes/pillar3_unified/gspa-pillar3.jar   # v2.1.0 build used for the Oct-15 LSO run
GO_OBO=/data/hohndor/gspa/reference/go.obo
API_BASE=http://unimatrix01.kaust.edu.sa:8000/v1                     # vLLM, served model qwen3.8-27b (Qwen3.8-27B-FP8)
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

mode=${1:-help}
case "$mode" in
  deploy)
    REPO=$(cd "$HERE/../../.." && pwd)
    ssh unimatrix01 "mkdir -p $PIPE $REPRO/slurm"
    rsync -a "$REPO/benchmark/unified/pillar3/run_pillar3_complex_lso.py" \
             "$HERE/p04_agent_wrapper.py" "$HERE/run.sh" "unimatrix01:$PIPE/"
    ;;
  submit)
    mkdir -p "$REPRO/slurm"
    sbatch -J p04_lso -o "$REPRO/slurm/p04_lso.out" -e "$REPRO/slurm/p04_lso.err" "$PIPE/run.sh" job "${2:-all}"
    ;;
  job)
    step=${2:-all}
    echo "host=$(hostname) job=${SLURM_JOB_ID:-none} step=$step start=$(date -Is)"
    # API key: the sidecar reads $BORG_LLM_API_KEY / $OPENAI_API_KEY or ~/.config/borg-llm/api-key (never printed)
    $PY "$PIPE/run_pillar3_complex_lso.py" --out-dir "$REPRO" --step "$step" --jar "$JAR" \
        --workers 8 --agent-wrapper "$PIPE/p04_agent_wrapper.py" --go-obo "$GO_OBO" --api-base "$API_BASE"
    echo "end=$(date -Is)"
    ;;
  fetch)
    REPO=$(cd "$HERE/../../.." && pwd)
    rsync -a "unimatrix01:$REPRO/pillar3_complex_lso.tsv" "unimatrix01:$REPRO/pillar3_complex_lso_mcnemar.tsv" "$REPO/results/"
    ;;
  *)
    sed -n '2,12p' "$0"
    ;;
esac
