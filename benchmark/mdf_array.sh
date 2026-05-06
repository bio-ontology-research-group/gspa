#!/bin/bash
# Per-array-task driver: run mDeepFRI on ONE panel genome, sequence-only mode.
#
# 13-genome PGAP-comparison panel:
#   1=hpylori 2=mgenitalium 3=mjannaschii
#   4=vcholerae 5=saureus 6=spneumoniae 7=ccrescentus 8=rprowazekii
#   9=tpallidum 10=tthermophilus 11=dradiodurans 12=scoelicolor 13=pfuriosus
set -euo pipefail

PANEL=(hpylori mgenitalium mjannaschii vcholerae saureus spneumoniae ccrescentus rprowazekii tpallidum tthermophilus dradiodurans scoelicolor pfuriosus)

idx=${SLURM_ARRAY_TASK_ID:-1}
tag="${PANEL[$((idx-1))]}"
[[ -n ${tag} ]] || { echo "no panel entry for idx=${idx}"; exit 1; }

# Source FAA: bench9 layout for bench9 genomes, bench10 layout otherwise.
faa=""
for cand in /data/hohndor/gspa/proteomes/${tag}.faa \
            /data/hohndor/gspa/proteomes/bench10/${tag}.faa; do
  [[ -s ${cand} ]] && { faa=${cand}; break; }
done
[[ -n ${faa} ]] || { echo "no fasta for ${tag}"; exit 2; }

out_dir=/data/hohndor/mdf-runs/${tag}
mkdir -p ${out_dir}

if [[ -s ${out_dir}/results.tsv ]]; then
  echo "[${idx}] ${tag}: already cached"
  exit 0
fi

echo "[${idx}] ${tag}: input=${faa} ($(grep -c '^>' ${faa}) proteins)"
/data/hohndor/envs/mdf-venv/bin/mDeepFRI predict-function \
  --skip-pdb \
  -i ${faa} \
  -o ${out_dir} \
  -w /data/hohndor/mdf-models-v1 \
  -t ${SLURM_CPUS_PER_TASK:-4} \
  > ${out_dir}/run.log 2>&1

n=$(wc -l < ${out_dir}/results.tsv 2>/dev/null || echo 0)
echo "[${idx}] ${tag}: done; ${n} result rows"

# Adapt to GSPA shape for benchmark_pgap_v2.py.
python3 /data/hohndor/gspa-v1.5.0/benchmark/parse_mdf_predictions.py \
  --mdf-csv ${out_dir}/results.tsv \
  --out ${out_dir}/${tag}_mdf_gspa.tsv
