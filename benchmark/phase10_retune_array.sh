#!/bin/bash
# Per-array-task driver: ONE (config, tag) pair.
# Reads task spec from a manifest file; uses SLURM_ARRAY_TASK_ID as line number.
set -euo pipefail

JAVA=/data/hohndor/envs/java21/bin/java
JAR=/data/hohndor/gspa-v1.5.0/gspa-cli/build/libs/gspa-1.5.0-SNAPSHOT.jar
REF=/data/hohndor/gspa/reference
ROOT=/data/hohndor/gspa/proteomes/bench10
MANIFEST=${ROOT}/phase10_retune/manifest.tsv

cd ${ROOT}

idx=${SLURM_ARRAY_TASK_ID:-1}
line=$(awk -v n=${idx} 'NR==n' ${MANIFEST})
[[ -n ${line} ]] || { echo "no manifest entry for idx=${idx}"; exit 1; }

cfg=$(echo "${line}" | cut -f1)
tag=$(echo "${line}" | cut -f2)
k=$(echo "${line}" | cut -f3)
flags=$(echo "${line}" | cut -f4)

claims=claims/${tag}_dpi_merged.jsonl
out_dir=phase10_retune/${cfg}
mkdir -p ${out_dir}
out_tsv=${out_dir}/${tag}_integrated.tsv
sug_tsv=${out_dir}/${tag}_suggestions.tsv
log=${out_dir}/${tag}.log

if [[ -s ${out_tsv} ]]; then
  echo "[${idx}] ${cfg}/${tag}: already cached"
  exit 0
fi

op_arg=""
[[ -s ops/${tag}_operons_up.tsv ]] && op_arg="--operons ops/${tag}_operons_up.tsv"
gp_arg=""
if   [[ -s gaps/${tag}_real.jsonl ]]; then gp_arg="--gaps gaps/${tag}_real.jsonl"
elif [[ -s gaps/${tag}_gaps.jsonl ]]; then gp_arg="--gaps gaps/${tag}_gaps.jsonl"
fi

echo "[${idx}] ${cfg}/${tag} (${k}) flags=${flags}"
${JAVA} -jar ${JAR} integrate \
  --claims ${claims} --out ${out_tsv} \
  --suggestions-out ${sug_tsv} \
  --go-owl ${REF}/go.owl --lite \
  --essential-profile ${k} \
  --pathways ${REF}/kegg_pathways.tsv --ec2go ${REF}/ec2go.txt \
  --enable-priors essentiality,coherence,gap_filling,genomic_context \
  ${op_arg} ${gp_arg} \
  ${flags} > ${log} 2>&1

echo "[${idx}] done; rows=$(wc -l < ${out_tsv})"
