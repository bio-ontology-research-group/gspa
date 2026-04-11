#!/bin/bash
# Run integration with the full prior stack using:
#   - go-plus.owl for taxon constraints
#   - KEGG pathway DB
#   - per-genome operons
#   - per-genome gapseq gaps (if available)
#   - essentiality, coherence, consistency, gap_filling, genomic_context
set -u
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics
JAVA=/data/hohndor/envs/java21/bin/java
JAR=/data/hohndor/gspa/bin/gspa.jar
ROOT=/data/hohndor/gspa/proteomes
REF=/data/hohndor/gspa/reference
BENCH=${ROOT}/bench9
OUT=${BENCH}/full_priors
mkdir -p ${OUT}

GENOMES=(ecoli hpylori mgenitalium mjannaschii ecolo157 bsubtilis mtb synechocystis paeruginosa)
KINGDOMS=(bacteria bacteria bacteria archaea bacteria bacteria bacteria bacteria bacteria)

for i in "${!GENOMES[@]}"; do
  tag=${GENOMES[$i]}
  k=${KINGDOMS[$i]}
  claims=${BENCH}/${tag}_claims.jsonl
  [[ -s ${claims} ]] || continue
  echo "=== ${tag} (${k}) ==="

  gaps_arg=""
  if [[ -s ${BENCH}/gapseq/${tag}_gaps.jsonl ]]; then
    gaps_arg="--gaps ${BENCH}/gapseq/${tag}_gaps.jsonl"
  fi
  operons_arg=""
  if [[ -s ${ROOT}/operons/${tag}_operons.tsv ]]; then
    operons_arg="--operons ${ROOT}/operons/${tag}_operons.tsv"
  fi

  ${JAVA} -jar ${JAR} integrate \
    --claims ${claims} \
    --out ${OUT}/${tag}_integrated.tsv \
    --go-owl ${REF}/go-plus.owl --lite \
    --essential-profile ${k} \
    --pathways ${REF}/kegg_pathways.tsv \
    --ec2go ${REF}/ec2go.txt \
    --reasoner-cache ${REF}/reasoner-cache \
    ${operons_arg} ${gaps_arg} \
    --enable-priors essentiality,coherence,consistency,gap_filling,genomic_context \
    > ${OUT}/${tag}_integrate.log 2>&1
  rc=$?
  echo "  rc=${rc}  $(wc -l <${OUT}/${tag}_integrated.tsv 2>/dev/null || echo 0) lines"
  grep -E 'Pathways|Operons|gaps|taxon|priors_fired|SAT|Converged' ${OUT}/${tag}_integrate.log | tail -5
done
echo DONE
