#!/bin/bash
# bench10_phase10.sh — Phase 10 Part 1 benchmark on the existing 10-genome PGAP set.
#
# Runs 5 configurations on the same cached predictor outputs:
#   C1 baseline                       — Phase 7 + one-shot Phase 8, gapseq on genome
#   C2 iterate                        — outer fixed-point, no clustering
#   C3 iterate + cluster 0.9          — with intra-genome MMseqs2 clustering
#   C4 iterate + cluster + blastp     — gapseq blastp against full proteome
#   C5 iterate + cluster + gapseq reps — gapseq blastp against cluster reps (Phase 11 preview)
#
# Plus one sensitivity row: C2 with --gapseq-pin-promotions false.
#
# Assumes predictor outputs + claims + truth files already exist under
# /data/hohndor/gspa/proteomes/bench10 (from bench10_pipeline.sh).
#
# Per-config subdirectories:
#   integrated_phase10/<config>/<tag>_integrated.tsv
#   integrated_phase10/<config>/<tag>_suggestions.tsv
#   integrated_phase10/<config>/<tag>.log
#
# Outputs:
#   results_phase10/<config>/<tag>_fmax.json
#   results_phase10/SUMMARY.tsv            — one row per (config, tag)
set -euo pipefail

source /storage/miniforge3/etc/profile.d/conda.sh
conda activate bench10

JAVA=/data/hohndor/envs/jdk-21.0.10+7/bin/java
JAR=/data/hohndor/gspa/bin/gspa.jar
REF=/data/hohndor/gspa/reference
BP=/data/hohndor/gspa/benchmark-py
ROOT=/data/hohndor/gspa/proteomes/bench10
cd ${ROOT}

TAGS=(vcholerae saureus spneumoniae ccrescentus rprowazekii tpallidum tthermophilus dradiodurans scoelicolor pfuriosus)
declare -A GENOME_KINGDOM=(
  [vcholerae]=bacteria [saureus]=bacteria [spneumoniae]=bacteria [ccrescentus]=bacteria
  [rprowazekii]=bacteria [tpallidum]=bacteria [tthermophilus]=bacteria [dradiodurans]=bacteria
  [scoelicolor]=bacteria [pfuriosus]=archaea
)

# ------------------------------------------------------------------
# Per-config flag sets (baseline is v1.0.0 behavior: no outer loop,
# no clustering, gapseq on genome — same as bench10_pipeline.sh STEP 9).
# ------------------------------------------------------------------
declare -A CONFIG_FLAGS
CONFIG_FLAGS[C1_baseline]=''
CONFIG_FLAGS[C2_iterate]='--dark-matter --iterate-gapseq'
CONFIG_FLAGS[C2_iterate_nopin]='--dark-matter --iterate-gapseq --gapseq-pin-promotions false'
CONFIG_FLAGS[C3_iter_cluster]='--dark-matter --iterate-gapseq --intragenome-cluster 0.9'
CONFIG_FLAGS[C4_iter_cluster_blastp]='--dark-matter --iterate-gapseq --intragenome-cluster 0.9 --gapseq-target proteome'
CONFIG_FLAGS[C5_iter_cluster_reps]='--dark-matter --iterate-gapseq --intragenome-cluster 0.9 --gapseq-target reps'

CONFIGS=(C1_baseline C2_iterate C3_iter_cluster C4_iter_cluster_blastp C5_iter_cluster_reps C2_iterate_nopin)

# ------------------------------------------------------------------
# STEP 1: Run integrate for each (config, genome)
# ------------------------------------------------------------------
mkdir -p integrated_phase10 results_phase10

for cfg in "${CONFIGS[@]}"; do
  flags="${CONFIG_FLAGS[$cfg]}"
  out_dir=integrated_phase10/${cfg}
  mkdir -p ${out_dir}
  echo "========================================"
  echo "CONFIG ${cfg}: ${flags:-(no Phase 10 flags)}"
  echo "========================================"

  for tag in "${TAGS[@]}"; do
    echo "--- ${cfg} / ${tag} ---"
    k=${GENOME_KINGDOM[$tag]}
    # Pick whichever merged claims file actually exists.
    claims=""
    for cand in claims/${tag}_dpi_merged.jsonl claims/${tag}_merged.jsonl claims/${tag}_dp_claims.jsonl; do
      [[ -s ${cand} ]] && { claims=${cand}; break; }
    done
    out_tsv=${out_dir}/${tag}_integrated.tsv
    sug_tsv=${out_dir}/${tag}_suggestions.tsv
    log=${out_dir}/${tag}.log

    [[ -n ${claims} ]] || { echo "  no claims file found; skip"; continue; }
    [[ -s ${out_tsv} ]] && { echo "  already done (${out_tsv})"; continue; }

    # Operon / gap files are best-effort: run without them if missing.
    op_arg=""
    [[ -s benchmark/operons/${tag}_operons.tsv ]] && op_arg="--operons benchmark/operons/${tag}_operons.tsv"
    gp_arg=""
    [[ -s benchmark/gaps/${tag}_gaps.jsonl ]] && gp_arg="--gaps benchmark/gaps/${tag}_gaps.jsonl"

    ${JAVA} -jar ${JAR} integrate \
      --claims ${claims} --out ${out_tsv} \
      --suggestions-out ${sug_tsv} \
      --go-owl ${REF}/go.owl --lite \
      --essential-profile ${k} \
      --pathways ${REF}/kegg_pathways.tsv --ec2go ${REF}/ec2go.txt \
      --enable-priors essentiality,coherence,gap_filling,genomic_context \
      ${op_arg} ${gp_arg} \
      ${flags} \
      > ${log} 2>&1 || echo "  integrate failed (see ${log})"

    if [[ -s ${out_tsv} ]]; then
      echo "  wrote $(wc -l < ${out_tsv}) rows (suggestions: $( [[ -s ${sug_tsv} ]] && wc -l < ${sug_tsv} || echo 0 ))"
    else
      echo "  [NO OUTPUT]"
    fi
  done
done

# ------------------------------------------------------------------
# STEP 2: F-max per (config, tag) against GOA truth
# ------------------------------------------------------------------
echo "========================================"
echo "F-max evaluation"
echo "========================================"

SUMMARY=results_phase10/SUMMARY.tsv
echo -e "config\ttag\ttruth\tfmax_micro\tfmax_cafa\tcoverage\tn_proteins\touter_iters\tpromoted_total" > ${SUMMARY}

for cfg in "${CONFIGS[@]}"; do
  out_dir=integrated_phase10/${cfg}
  res_dir=results_phase10/${cfg}
  mkdir -p ${res_dir}

  for tag in "${TAGS[@]}"; do
    integrated=${out_dir}/${tag}_integrated.tsv
    truth_all=truth/${tag}_truth_all.tsv
    [[ -s ${integrated} ]] || continue
    [[ -s ${truth_all}  ]] || continue

    json=${res_dir}/${tag}_fmax.json
    python3 ${BP}/benchmark_pgap_v2.py \
      --gspa ${integrated} \
      --truth ${truth_all} \
      --genome ${tag} \
      --json-out ${json} \
      > ${res_dir}/${tag}_fmax.log 2>&1 || true

    # Also scrape outer-loop trace from the integrate log (best effort).
    log=${out_dir}/${tag}.log
    outer_iters=$(grep -oE 'iter=[0-9]+' ${log} | head -1 | sed 's/iter=//')
    promoted_total=$(grep -oE 'promoted_per_iter=\[[^]]*\]' ${log} | head -1)

    fmax_micro=$(python3 -c "import json,sys; d=json.load(open('${json}')); print(d.get('fmax_overall',{}).get('fmax',''))" 2>/dev/null || echo '')
    fmax_cafa=$(python3 -c "import json,sys; d=json.load(open('${json}')); print(d.get('fmax_cafa_overall',{}).get('fmax',''))" 2>/dev/null || echo '')
    cov=$(python3 -c "import json,sys; d=json.load(open('${json}')); print(d.get('coverage',''))" 2>/dev/null || echo '')
    npro=$(python3 -c "import json,sys; d=json.load(open('${json}')); print(d.get('n_proteins',''))" 2>/dev/null || echo '')

    echo -e "${cfg}\t${tag}\tall\t${fmax_micro}\t${fmax_cafa}\t${cov}\t${npro}\t${outer_iters:-0}\t${promoted_total:-}" >> ${SUMMARY}
  done
done

echo "========================================"
echo "Summary written to ${SUMMARY}"
echo "========================================"
column -t -s $'\t' ${SUMMARY} | head -20
