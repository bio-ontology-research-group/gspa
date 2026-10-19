#!/bin/bash
#SBATCH --job-name=p07_extract
#SBATCH --partition=batch
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --array=0-107%40
#SBATCH --output=/ibex/scratch/projects/c2014/hohndor/meteor108_repro_2026-10-19/logs/extract_%A_%a.out
# P07: per-genome draft / METEOR / Grow-Support extraction (read-only on the METEOR tree).
set -o pipefail
D=/ibex/scratch/projects/c2014/hohndor/meteor108_repro_2026-10-19
R=/ibex/scratch/projects/c2014/kexin/funcarve/psb_revision
module load python/3.11.0 >/dev/null 2>&1
source $D/venv/bin/activate
export METEOR_DATA=$R/code_snapshot/data
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $R/code_snapshot/panel108_gram.tsv)
GCA=$(echo "$LINE" | cut -f1); GRAM=$(echo "$LINE" | cut -f2)
echo "START $GCA $GRAM $(date)"
python3 -u $D/scripts/extract_meteor108_drafts.py --gca "$GCA" --gram "$GRAM" --outdir $D/per_genome
echo "END $GCA rc=$? $(date)"
