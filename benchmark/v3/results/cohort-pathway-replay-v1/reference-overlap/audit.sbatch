#!/bin/bash
#SBATCH --job-name=gspa-reference-overlap
#SBATCH --partition=debug
#SBATCH --exclude=node001,node005,node008,unimatrix
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/data/hohndor/gspa-v3/cohort-reported-overlap-v1/audit-%j.log
set -euo pipefail
export TMPDIR="$SLURM_TMPDIR"
cd /data/hohndor/gspa-v3/cohort-reported-overlap-v1
sha256sum -c SHA256SUMS
/usr/bin/python3 verify_reported_reference_overlap.py --plan plan.json --plan-sha256 35d2b0faa992530ab7cc5d20c12871a59909e20caea7a1c5742324489e202c8d --out "$SLURM_TMPDIR/overlap.json"
sha256sum -c SHA256SUMS
cp "$SLURM_TMPDIR/overlap.json" overlap.json
