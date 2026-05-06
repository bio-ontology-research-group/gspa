#!/bin/bash
# Run integration with the full prior stack:
#   - go-plus.owl for taxon constraints
#   - KEGG pathway DB
#   - per-genome operons (heuristic intergenic-distance caller)
#   - per-genome gapseq gaps (if available)
#   - essentiality, coherence, consistency, gap_filling, genomic_context
#
# Usage:
#   run_integrate_full_priors.sh [--dry-run]
#
# Env-var overrides (useful for testing or ad-hoc layouts):
#   JAVA, JAR, ROOT, REF, BENCH, OUT, GENOMES, KINGDOMS
set -u

# ---- arg parsing --------------------------------------------------------
DRY_RUN=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run)       DRY_RUN=1; shift;;
    -h|--help)
      echo "Usage: $0 [--dry-run]"
      echo "Env: JAVA JAR ROOT REF BENCH OUT GENOMES KINGDOMS"
      exit 0;;
    *) echo "Unknown arg: $1" >&2; exit 2;;
  esac
done

# ---- environment (overrideable) -----------------------------------------
# Conda activation is optional so the script works on dev hosts where the
# cluster's conda installation is absent (e.g. when running --dry-run).
if [[ -f /storage/miniforge3/etc/profile.d/conda.sh ]]; then
  source /storage/miniforge3/etc/profile.d/conda.sh
  conda activate metagenomics
fi

JAVA=${JAVA:-/data/hohndor/envs/java21/bin/java}
JAR=${JAR:-/data/hohndor/gspa/bin/gspa.jar}
ROOT=${ROOT:-/data/hohndor/gspa/proteomes}
REF=${REF:-/data/hohndor/gspa/reference}
BENCH=${BENCH:-${ROOT}/bench9}
OUT=${OUT:-${BENCH}/full_priors}

# Default 9-genome benchmark set; override via env for ad-hoc runs.
GENOMES_DEFAULT='ecoli hpylori mgenitalium mjannaschii ecolo157 bsubtilis mtb synechocystis paeruginosa'
KINGDOMS_DEFAULT='bacteria bacteria bacteria archaea bacteria bacteria bacteria bacteria bacteria'
read -r -a GENOMES <<<"${GENOMES:-$GENOMES_DEFAULT}"
read -r -a KINGDOMS <<<"${KINGDOMS:-$KINGDOMS_DEFAULT}"

mkdir -p "${OUT}"

echo "out-dir: ${OUT}"
[[ ${DRY_RUN} -eq 1 ]] && echo "(dry-run: command will be printed, not executed)"

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

  op_path=${ROOT}/operons/${tag}_operons.tsv
  operons_arg=""
  if [[ -s ${op_path} ]]; then
    operons_arg="--operons ${op_path}"
  fi
  echo "  operons=${op_path}"

  cmd=( "${JAVA}" -jar "${JAR}" integrate
        --claims "${claims}"
        --out "${OUT}/${tag}_integrated.tsv"
        --go-owl "${REF}/go-plus.owl" --lite
        --essential-profile "${k}"
        --pathways "${REF}/kegg_pathways.tsv"
        --ec2go "${REF}/ec2go.txt"
        --reasoner-cache "${REF}/reasoner-cache"
        ${operons_arg} ${gaps_arg}
        --enable-priors essentiality,coherence,consistency,gap_filling,genomic_context )

  if [[ ${DRY_RUN} -eq 1 ]]; then
    echo "  cmd: ${cmd[*]}"
    continue
  fi

  "${cmd[@]}" > "${OUT}/${tag}_integrate.log" 2>&1
  rc=$?
  echo "  rc=${rc}  $(wc -l <"${OUT}/${tag}_integrated.tsv" 2>/dev/null || echo 0) lines"
  grep -E 'Pathways|Operons|gaps|taxon|priors_fired|SAT|Converged' "${OUT}/${tag}_integrate.log" | tail -5
done
echo DONE
