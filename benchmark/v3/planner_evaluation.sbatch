#!/bin/bash
# Equal-budget condition-planner benchmark: CPU-only, no GPU. Prepared for the root; NOT submitted by the author.
#
# What needs a cluster: nothing computational (the exhaustive sweep and all deterministic arms take minutes on 2 CPUs).
# What is slow: actual Qwen calls to the unimatrix endpoint (observed ~110 s per request under load, ~6.6k prompt tokens,
# 8 requests per gene at budget 8). Four shard workers keep the endpoint load modest; re-run the same command to resume
# (finished (gene, arm, budget, seed) rows are skipped).
#
# Secrets: the endpoint key is read by metabolic_experiments.load_api_key() from BORG_LLM_API_KEY / OPENAI_API_KEY or
# ~/.config/borg-llm/api-key. Do not put it in this file.
#
#   sbatch benchmark/v3/planner_evaluation.sbatch
#
#SBATCH --job-name=gspa-planner-eval
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=planner-eval-%j.log

set -euo pipefail
REPO=${REPO:-$PWD}
DATA=${DATA:-/tmp/gspa-v3-metabolic-data}                       # prepare_metabolic_benchmarks.py output (see docs)
OUT=${OUT:-$REPO/benchmark/v3/results/planner-evaluation-cluster}
PY=${PY:-python3}                                               # needs numpy, scipy, python-libsbml, cobra not required
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 GSPA_LLM_URL=${GSPA_LLM_URL:-http://unimatrix01.kaust.edu.sa:8000/v1}

CASES=$DATA/derived/cases_manuscript_regulatory.jsonl
MODEL=$DATA/raw/iML1515.xml
RUN="$PY -I $REPO/benchmark/v3/run_planner_benchmark.py --cases $CASES --model $MODEL --out $OUT"

mkdir -p "$OUT"
$RUN --stages sweep,deterministic
$PY -I $REPO/benchmark/v3/run_metabolic_benchmark.py --cases "$CASES" --model "$MODEL" --out "$OUT/phenotype-reference"

pids=()
for i in 0 1 2 3; do
  $RUN --stages qwen --shard $i/4 > "$OUT/qwen.shard$i.log" 2>&1 &
  pids+=($!)
done
status=0
for p in "${pids[@]}"; do wait "$p" || status=1; done

$RUN --stages report --phenotype-reference "$OUT/phenotype-reference/report.json"
exit $status
