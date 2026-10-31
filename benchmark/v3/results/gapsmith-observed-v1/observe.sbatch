#!/bin/bash
#SBATCH --job-name=gspa-gapsmith-observed
#SBATCH --partition=debug
#SBATCH --nodelist=node004
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#SBATCH --output=/data/hohndor/gspa-v3/gapsmith-observed-v1/observe-%j.log
set -euo pipefail
export TMPDIR="$SLURM_TMPDIR"
cd /data/hohndor/gspa-v3/gapsmith-observed-v1
/usr/bin/python3 observe_gapsmith.py \
  --root /data/hohndor/gspa-v3/gapsmith-observed-v1 \
  --reference /data/hohndor/gspa/bin/gapsmith/data_merged \
  --binary /data/hohndor/gspa/bin/gapsmith/gapsmith-linux \
  --diamond /storage/miniforge3/envs/metagenomics/bin/diamond \
  --proteome-sha256 72cac5531a5b22ddba171c0bd3c79da968a9435980b60dd4aa73431d6392b476
