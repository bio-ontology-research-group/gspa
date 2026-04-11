#!/bin/bash
# F-max comparison: full-priors GSPA vs combined (no priors) vs PGAP.
set -u
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics
ROOT=/data/hohndor/gspa/proteomes
BENCH=${ROOT}/bench9
FULL=${BENCH}/full_priors
BP=/data/hohndor/gspa/benchmark-py

GENOMES=(ecoli hpylori mgenitalium mjannaschii ecolo157 bsubtilis mtb synechocystis paeruginosa)

for tag in "${GENOMES[@]}"; do
  echo "=== ${tag} ==="
  truth_exp=${ROOT}/truth_dual/${tag}_truth_exp.tsv
  truth_all=${ROOT}/truth_dual/${tag}_truth_all.tsv
  combined=${BENCH}/${tag}_integrated.tsv
  full=${FULL}/${tag}_integrated.tsv
  [[ -s ${combined} && -s ${full} ]] || { echo "  skip: missing files"; continue; }

  pgap_arg=""
  if [[ -s ${ROOT}/${tag}_pgap.tsv && $(wc -l <${ROOT}/${tag}_pgap.tsv) -gt 1 ]]; then
    pgap_arg="--pgap ${ROOT}/${tag}_pgap.tsv"
  fi
  map_arg=""
  [[ -s ${ROOT}/${tag}.refseq_to_uniprot.tsv ]] && map_arg="--gspa-key-map ${ROOT}/${tag}.refseq_to_uniprot.tsv"

  python3 ${BP}/benchmark_pgap_v2.py \
    --truth exp:${truth_exp} --truth all:${truth_all} \
    --gspa ${combined} --gspa-priors ${full} ${pgap_arg} ${map_arg} \
    --tag ${tag} --n-bootstrap 200 \
    > ${FULL}/${tag}_fmax.json 2> ${FULL}/${tag}_fmax.err
  echo "  rc=$?"
done
echo DONE
