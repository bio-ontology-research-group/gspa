#!/bin/bash
#SBATCH --job-name=gspa-cohort-pathways
#SBATCH --partition=debug
#SBATCH --exclude=node001,node005,node008,unimatrix
#SBATCH --array=0-17%4
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#SBATCH --output=/data/hohndor/gspa-v3/cohort-pathway-search-v1/logs/search-%A_%a.log
set -euo pipefail
export TMPDIR="$SLURM_TMPDIR"
cd /data/hohndor/gspa-v3/cohort-pathway-search-v1
/usr/bin/python3 observe_gapsmith_cohort.py --manifest manifest.json --manifest-sha256 b9491693f7b34b29abe99c7f077c44f2b70b3003db258816d487d3badf2e8fa9 --index "$SLURM_ARRAY_TASK_ID"
