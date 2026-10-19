#!/bin/bash
#SBATCH --job-name=p07_merge
#SBATCH --partition=batch
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/hohndor/meteor108_repro_2026-10-19/logs/merge_%j.out
set -eo pipefail
D=/ibex/scratch/projects/c2014/hohndor/meteor108_repro_2026-10-19
module load python/3.11.0 >/dev/null 2>&1
source $D/venv/bin/activate
export METEOR_DATA=/ibex/scratch/projects/c2014/kexin/funcarve/psb_revision/code_snapshot/data
python3 -u $D/scripts/merge_meteor108_drafts.py
