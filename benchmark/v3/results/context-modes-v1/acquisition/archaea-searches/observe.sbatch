#!/bin/bash
#SBATCH --job-name=gspa-archaea-pathways
#SBATCH --partition=debug
#SBATCH --exclude=node001,node005,node008,unimatrix
#SBATCH --array=0-9%4
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#SBATCH --output=/data/hohndor/gspa-v3/cohort-pathway-search-archaea-v2/logs/search-%A_%a.log
set -euo pipefail
export TMPDIR="$SLURM_TMPDIR"
cd /data/hohndor/gspa-v3/cohort-pathway-search-archaea-v2
/usr/bin/python3 observe_gapsmith_cohort.py --manifest manifest.json --manifest-sha256 00f67553b10f89fc288e28eb9f113e3e569163813d433728cf2b4297b21916f5 --index "$SLURM_ARRAY_TASK_ID"
