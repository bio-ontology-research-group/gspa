#!/bin/bash
#SBATCH --job-name=p07_merge
#SBATCH --partition=batch
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/encrypted/e3008/gspa/meteor108_repro_2026-10-20/logs/merge_%j.out
# P07 merge of the per-genome extractions into meteor108_drafts_compact.pkl.gz (GSPA code only).
set -eo pipefail
W=${W:-/encrypted/e3008/gspa/meteor108_repro_2026-10-20}
VENV=${VENV:-/ibex/scratch/projects/c2014/hohndor/meteor108_repro_2026-10-19/venv}
module load python/3.11.0 >/dev/null 2>&1
source $VENV/bin/activate
python3 -u $W/scripts/merge_meteor108_drafts.py --work-dir $W --data-dir $W/data --panel $W/scripts/panel108_gram.tsv
