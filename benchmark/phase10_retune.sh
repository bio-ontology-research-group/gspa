#!/bin/bash
# phase10_retune.sh — v1.5.0 Phase C retune of the iterative dark-matter loop.
#
# Variables: qBase ∈ {0.50 (legacy default), 0.70, 0.75}.
# Held constant: claims (bench10 dpi_merged), priors, operons, gaps, jar.
#
# Per-genome gap source:
#   - real gapseq output if gaps/${tag}_real.jsonl exists (ccrescentus,
#     rprowazekii, dradiodurans, tthermophilus)
#   - fall back to synthetic 400-gap output gaps/${tag}_gaps.jsonl
#     for genomes that hit the documented zero-byte Reactions.tbl bug
#     (vcholerae, saureus, spneumoniae, tpallidum, scoelicolor, pfuriosus)
#
# Plus the C1 baseline (no --iterate-gapseq) for direct comparison
# against the existing benchmark/RESULTS.md numbers.
#
# Outputs:
#   phase10_retune/<config>/<tag>_integrated.tsv
#   phase10_retune/<config>/<tag>_suggestions.tsv
#   phase10_retune/<config>/<tag>.log
#   phase10_retune/results/<config>/<tag>_fmax.json
#   phase10_retune/SUMMARY.tsv
set -euo pipefail

JAVA=/data/hohndor/envs/java21/bin/java
JAR=/data/hohndor/gspa-v1.5.0/gspa-cli/build/libs/gspa-1.5.0-SNAPSHOT.jar
REF=/data/hohndor/gspa/reference
BP=/data/hohndor/gspa/benchmark-py
ROOT=/data/hohndor/gspa/proteomes/bench10
cd ${ROOT}

TAGS=(vcholerae saureus spneumoniae ccrescentus rprowazekii tpallidum tthermophilus dradiodurans scoelicolor pfuriosus)
declare -A KIN=(
  [vcholerae]=bacteria [saureus]=bacteria [spneumoniae]=bacteria [ccrescentus]=bacteria
  [rprowazekii]=bacteria [tpallidum]=bacteria [tthermophilus]=bacteria [dradiodurans]=bacteria
  [scoelicolor]=bacteria [pfuriosus]=archaea
)

declare -A FLAGS
FLAGS[C1_baseline]=''
FLAGS[C2_q050]='--dark-matter --iterate-gapseq --gapseq-q-base 0.50'
FLAGS[C2_q070]='--dark-matter --iterate-gapseq --gapseq-q-base 0.70'
FLAGS[C2_q075]='--dark-matter --iterate-gapseq --gapseq-q-base 0.75'

CONFIGS=(C1_baseline C2_q050 C2_q070 C2_q075)

mkdir -p phase10_retune/results

# STEP 1 — integrate per (config, tag)
for cfg in "${CONFIGS[@]}"; do
  flags="${FLAGS[$cfg]}"
  out_dir=phase10_retune/${cfg}
  mkdir -p ${out_dir}
  echo "============================="
  echo "CONFIG ${cfg}: ${flags:-(no Phase 10 flags)}"
  echo "============================="

  for tag in "${TAGS[@]}"; do
    k=${KIN[$tag]}
    claims=claims/${tag}_dpi_merged.jsonl
    out_tsv=${out_dir}/${tag}_integrated.tsv
    sug_tsv=${out_dir}/${tag}_suggestions.tsv
    log=${out_dir}/${tag}.log

    [[ -s ${claims} ]] || { echo "  ${tag}: no claims; skip"; continue; }
    [[ -s ${out_tsv} ]] && { echo "  ${tag}: cached"; continue; }

    op_arg=""
    [[ -s ops/${tag}_operons_up.tsv ]] && op_arg="--operons ops/${tag}_operons_up.tsv"

    # Real gapseq if available, else synthetic 400-gap fallback.
    gp_arg=""
    if   [[ -s gaps/${tag}_real.jsonl ]]; then gp_arg="--gaps gaps/${tag}_real.jsonl"
    elif [[ -s gaps/${tag}_gaps.jsonl ]]; then gp_arg="--gaps gaps/${tag}_gaps.jsonl"
    fi

    echo "  ${tag} (${k}) op=${op_arg:-none} gaps=${gp_arg:-none}"

    ${JAVA} -jar ${JAR} integrate \
      --claims ${claims} --out ${out_tsv} \
      --suggestions-out ${sug_tsv} \
      --go-owl ${REF}/go.owl --lite \
      --essential-profile ${k} \
      --pathways ${REF}/kegg_pathways.tsv --ec2go ${REF}/ec2go.txt \
      --enable-priors essentiality,coherence,gap_filling,genomic_context \
      ${op_arg} ${gp_arg} \
      ${flags} > ${log} 2>&1 \
      || echo "  ${tag}: integrate failed (see ${log})"
  done
done

# STEP 2 — F-max (micro + CAFA, 200-bootstrap) per (config, tag)
echo "============================="
echo "F-max scoring"
echo "============================="

SUMMARY=phase10_retune/SUMMARY.tsv
echo -e "config\ttag\tfmax_micro\tfmax_cafa\tn_truth\touter_iters\tpromoted_per_iter\tgap_source" > ${SUMMARY}

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
    outer_iters=$(grep -oE 'iter=[0-9]+' ${log} | head -1 | sed 's/iter=//')
    promoted=$(grep -oE 'promoted_per_iter=\[[^]]*\]' ${log} | head -1)
    gap_src='none'
    [[ -s gaps/${tag}_real.jsonl ]] && gap_src='real' || { [[ -s gaps/${tag}_gaps.jsonl ]] && gap_src='synthetic'; }

    fmax_micro=$(python3 -c "import json; d=json.load(open('${json}')); print(d['by_truth']['all']['results'][0]['fmax_overall'])" 2>/dev/null || echo '')
    fmax_cafa=$(python3 -c "import json; d=json.load(open('${json}')); print(d['by_truth']['all']['results'][0]['fmax_cafa_overall'])" 2>/dev/null || echo '')
    n_truth=$(python3 -c "import json; d=json.load(open('${json}')); print(d['by_truth']['all']['truth_proteins'])" 2>/dev/null || echo '')

    echo -e "${cfg}\t${tag}\t${fmax_micro}\t${fmax_cafa}\t${n_truth}\t${outer_iters:-0}\t${promoted}\t${gap_src}" >> ${SUMMARY}
  done
done

echo "============================="
echo "SUMMARY → ${SUMMARY}"
echo "============================="
column -t -s $'\t' ${SUMMARY} | head -45
