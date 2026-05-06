#!/bin/bash
# Phase D scoring: F-max micro + CAFA for each panel genome's mdF
# predictions, against the matching truth_*.tsv. Emits SUMMARY.tsv.
set -euo pipefail

BP=/data/hohndor/gspa/benchmark-py
RUN=/data/hohndor/mdf-runs

# Truth files: bench10 has them under .../bench10/truth/ for the 10 PGAP
# genomes; bench9 (hpylori, mgenitalium, mjannaschii) has them under a
# parallel layout at /data/hohndor/gspa/proteomes/bench9/truth_dual or
# /data/hohndor/gspa/proteomes/bench9/<tag>_truth_*.tsv. Try both.
truth_for() {
  local tag=$1
  for cand in \
      /data/hohndor/gspa/proteomes/bench10/truth/${tag}_truth_all.tsv \
      /data/hohndor/gspa/proteomes/truth_dual/${tag}_truth_all.tsv \
      /data/hohndor/gspa/proteomes/truth/${tag}_truth.tsv \
      /data/hohndor/gspa/proteomes/bench_gtdb30/truth/${tag}_truth_all.tsv; do
    if [ -s "$cand" ]; then echo "$cand"; return 0; fi
  done
  return 0
}

mkdir -p ${RUN}/scoring
SUMMARY=${RUN}/scoring/SUMMARY.tsv
echo -e "tag\tn_truth\tn_pred\tfmax_micro\tfmax_cafa\ttruth_path" > ${SUMMARY}

for tag in hpylori mgenitalium mjannaschii \
           vcholerae saureus spneumoniae ccrescentus rprowazekii \
           tpallidum tthermophilus dradiodurans scoelicolor pfuriosus; do
  pred=${RUN}/${tag}/${tag}_mdf_gspa.tsv
  truth=$(truth_for ${tag})
  [[ -s ${pred}  ]] || { echo -e "${tag}\t-\t-\t-\t-\t<no pred>" >> ${SUMMARY}; continue; }
  [[ -s ${truth} ]] || { echo -e "${tag}\t-\t-\t-\t-\t<no truth>" >> ${SUMMARY}; continue; }

  json=${RUN}/scoring/${tag}_fmax.json
  python3 ${BP}/benchmark_pgap_v2.py \
    --gspa ${pred} --truth all:${truth} --tag ${tag} --n-bootstrap 200 \
    > ${json} 2> ${RUN}/scoring/${tag}_fmax.err || true

  n_pred=$(tail -n +2 ${pred} | wc -l)
  n_truth=$(tail -n +2 ${truth} | awk '{print $1}' | sort -u | wc -l)
  fmax_micro=$(python3 -c "import json; d=json.load(open('${json}')); print(d['by_truth']['all']['results'][0]['fmax_overall'])" 2>/dev/null || echo '')
  fmax_cafa=$(python3 -c "import json; d=json.load(open('${json}')); print(d['by_truth']['all']['results'][0]['fmax_cafa_overall'])" 2>/dev/null || echo '')
  echo -e "${tag}\t${n_truth}\t${n_pred}\t${fmax_micro}\t${fmax_cafa}\t${truth}" >> ${SUMMARY}
done

echo "==SUMMARY=="
column -t -s $'\t' ${SUMMARY}
