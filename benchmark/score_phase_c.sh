#!/bin/bash
# Phase C scoring: F-max micro + CAFA per (config, tag) for the
# 10 PGAP genomes, against bench10/truth/*_truth_all.tsv.
set -euo pipefail

BP=/data/hohndor/gspa/benchmark-py
ROOT=/data/hohndor/gspa/proteomes/bench10
cd ${ROOT}

CONFIGS=(C1_baseline C2_q050 C2_q070 C2_q075)
TAGS=(vcholerae saureus spneumoniae ccrescentus rprowazekii tpallidum tthermophilus dradiodurans scoelicolor pfuriosus)

SUMMARY=phase10_retune/SUMMARY.tsv
echo -e "config\ttag\tfmax_micro\tfmax_cafa\tn_truth\touter_iters\tpromoted_per_iter\tgap_source" > ${SUMMARY}

mkdir -p phase10_retune/results
for cfg in "${CONFIGS[@]}"; do
  out_dir=phase10_retune/${cfg}
  res_dir=phase10_retune/results/${cfg}
  mkdir -p ${res_dir}
  for tag in "${TAGS[@]}"; do
    integrated=${out_dir}/${tag}_integrated.tsv
    truth_all=truth/${tag}_truth_all.tsv
    [[ -s ${integrated} ]] || continue
    [[ -s ${truth_all}  ]] || continue

    json=${res_dir}/${tag}_fmax.json
    python3 ${BP}/benchmark_pgap_v2.py \
      --gspa ${integrated} --truth all:${truth_all} --tag ${tag} --n-bootstrap 200 \
      > ${json} 2> ${res_dir}/${tag}_fmax.err || true

    log=${out_dir}/${tag}.log
    outer_iters=$(grep -oE 'iter=[0-9]+' ${log} 2>/dev/null | head -1 | sed 's/iter=//' || true)
    promoted=$(grep -oE 'promoted_per_iter=\[[^]]*\]' ${log} 2>/dev/null | head -1 || true)
    gap_src='none'
    [[ -s gaps/${tag}_real.jsonl ]] && gap_src='real' || { [[ -s gaps/${tag}_gaps.jsonl ]] && gap_src='synthetic'; }

    fmax_micro=$(python3 -c "import json; d=json.load(open('${json}')); print(d['by_truth']['all']['results'][0]['fmax_overall'])" 2>/dev/null || echo '')
    fmax_cafa=$(python3 -c "import json; d=json.load(open('${json}')); print(d['by_truth']['all']['results'][0]['fmax_cafa_overall'])" 2>/dev/null || echo '')
    n_truth=$(python3 -c "import json; d=json.load(open('${json}')); print(d['by_truth']['all']['truth_proteins'])" 2>/dev/null || echo '')

    echo -e "${cfg}\t${tag}\t${fmax_micro}\t${fmax_cafa}\t${n_truth}\t${outer_iters:-0}\t${promoted}\t${gap_src}" >> ${SUMMARY}
  done
done

echo "===SUMMARY==="
column -t -s $'\t' ${SUMMARY}
