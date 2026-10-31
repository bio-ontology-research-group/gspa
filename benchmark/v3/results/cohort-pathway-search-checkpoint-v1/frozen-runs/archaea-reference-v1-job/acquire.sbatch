#!/bin/bash
#SBATCH --job-name=gspa-archaea-reference
#SBATCH --partition=debug
#SBATCH --exclude=node001,node005,node008,unimatrix
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/data/hohndor/gspa-v3/archaea-reference-v1-job/acquire-%j.log
set -euo pipefail
export TMPDIR="$SLURM_TMPDIR"
cd /data/hohndor/gspa-v3/archaea-reference-v1-job
sha256sum -c SHA256SUMS
/usr/bin/python3 acquire_gapsmith_archaea.py --root /data/hohndor/gspa-v3/archaea-reference-v1 --tables /data/hohndor/gspa/bin/gapsmith/data_merged
sha256sum -c SHA256SUMS
