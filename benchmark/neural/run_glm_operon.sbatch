#!/usr/bin/env bash
#SBATCH --job-name=gspa-glm-operon
#SBATCH --partition=pi-hohndor
#SBATCH --qos=pi-hohndor
#SBATCH --gres=gpu:1
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH -t 06:00:00
#SBATCH -o /mnt/data/u/hohndor/gspa-glm/logs/glm-operon-%A_%a.out
#SBATCH -e /mnt/data/u/hohndor/gspa-glm/logs/glm-operon-%A_%a.err

# Run gLM operon caller on a single genome, selected by SLURM_ARRAY_TASK_ID
# (1-based row index into the manifest). For the freecycle partition,
# replace the partition/qos lines with:
#     #SBATCH --partition=freecycle
#     #SBATCH --qos=freecycle
#     #SBATCH --gres=gpu:h200:1
#     #SBATCH --requeue
# (per orix-workbench/README.md).

set -euo pipefail

# ---- env ----------------------------------------------------------------
# Pixi is the cluster-preferred env manager (see orix-workbench/CLAUDE.md).
# Until a pixi recipe is committed, fall back to conda + a project env that
# carries gLM, esm-2, torch, numpy.
if command -v pixi >/dev/null && [[ -n "${PIXI_PROJECT_ROOT:-}" ]]; then
    eval "$(pixi shell-hook)"
elif [[ -f /home/hohndor/miniforge3/etc/profile.d/conda.sh ]]; then
    source /home/hohndor/miniforge3/etc/profile.d/conda.sh
    conda activate gspa-glm
fi
export PYTHONUNBUFFERED=1

# ---- paths --------------------------------------------------------------
ROOT=/mnt/data/u/hohndor/gspa-glm
NEURAL=$(realpath "$(dirname "$0")")        # benchmark/neural/ on whatever node
WEIGHTS=/mnt/data/u/hohndor/gLM/weights     # download once; pin a commit SHA
TORCH_HOME=/mnt/data/u/hohndor/esm2         # ESM2-650M cache
MANIFEST=${MANIFEST:-$ROOT/manifest.tsv}    # tag<TAB>fasta_path<TAB>gff_path<TAB>output_dir
MODE=${MODE:-real}                          # real | mock | self-test

mkdir -p "$ROOT/logs"
export TORCH_HOME

# ---- pick a row ---------------------------------------------------------
ROW=$(( ${SLURM_ARRAY_TASK_ID:-1} + 1 ))    # +1 to skip header
read -r TAG FASTA GFF OUT < <(awk -F'\t' -v r=$ROW 'NR==r{print $1, $2, $3, $4}' "$MANIFEST")

if [[ -z "${TAG:-}" || ! -s "$FASTA" || ! -s "$GFF" ]]; then
    echo "manifest row $ROW invalid: tag=$TAG fasta=$FASTA gff=$GFF" >&2
    exit 2
fi
mkdir -p "$OUT"

N_PROT=$(grep -c '^>' "$FASTA")
echo "=== $TAG ($N_PROT proteins; mode=$MODE) ==="
date

python "$NEURAL/run_glm_operon.py" \
    --mode "$MODE" \
    --fasta "$FASTA" --gff "$GFF" \
    --weights "$WEIGHTS" \
    --operons-out "$OUT/operons.tsv" \
    --confidence-out "$OUT/operons_confidence.tsv" \
    --centroids-out "$OUT/operons_centroids.npz" \
    --protein-embeddings-out "$OUT/protein_embeddings.npz" \
    --device cuda

echo "=== $TAG output ==="
wc -l "$OUT/operons.tsv" "$OUT/operons_confidence.tsv"
ls -la "$OUT"
echo DONE "$TAG"; date
