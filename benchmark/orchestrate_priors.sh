#!/bin/bash
# Orchestrator: wait for gapseq_master.log to report "ALL GAPSEQ DONE",
# then parse gaps, run full-priors integration, compute F-max, compile.
set -u
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics

ROOT=/data/hohndor/gspa/proteomes
BENCH=${ROOT}/bench9
REF=/data/hohndor/gspa/reference
BP=/data/hohndor/gspa/benchmark-py

echo "[orchestrator] waiting for gapseq to finish..."
while true; do
  if grep -q 'ALL GAPSEQ DONE' ${ROOT}/gapseq_master.log 2>/dev/null; then
    break
  fi
  sleep 60
done
echo "[orchestrator] gapseq done. Parsing gaps..."

mkdir -p ${BENCH}/gapseq
for tag in ecoli hpylori mgenitalium mjannaschii ecolo157 bsubtilis mtb synechocystis paeruginosa; do
  pt=${ROOT}/gapseq/${tag}/${tag}-all-Pathways.tbl
  rt=${ROOT}/gapseq/${tag}/${tag}-all-Reactions.tbl
  if [[ -s ${pt} && -s ${rt} ]]; then
    python3 ${BP}/parse_gapseq_gaps.py \
      --pathways-tbl ${pt} --reactions-tbl ${rt} \
      --ec2go ${REF}/ec2go.txt \
      --tag ${tag} \
      --out ${BENCH}/gapseq/${tag}_gaps.jsonl
  else
    echo "  ${tag}: no gapseq output, skip"
  fi
done

echo "[orchestrator] running full-priors integration..."
bash ${ROOT}/run_integrate_full_priors.sh > ${BENCH}/full_priors_integrate.log 2>&1

echo "[orchestrator] computing F-max (combined vs full-priors vs PGAP)..."
bash ${ROOT}/run_fmax_full_priors.sh > ${BENCH}/full_priors_fmax.log 2>&1

echo "[orchestrator] compiling comparison..."
python3 ${BP}/compile_ablation.py \
  --dir ${BENCH}/full_priors \
  --genomes ecoli hpylori mgenitalium mjannaschii ecolo157 bsubtilis mtb synechocystis paeruginosa \
  --configs diamond pfam combined priors \
  > ${BENCH}/FULL_PRIORS_REPORT.txt 2>&1 || true

# Pretty-print the F-max JSONs from the full_priors dir directly
python3 ${BP}/print_fmax.py ${BENCH}/full_priors/*_fmax.json > ${BENCH}/FULL_PRIORS_FMAX_SUMMARY.txt 2>&1 || true

echo "[orchestrator] DONE"
