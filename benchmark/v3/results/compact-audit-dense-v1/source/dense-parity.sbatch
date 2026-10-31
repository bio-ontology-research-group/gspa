#!/bin/bash
#SBATCH --job-name=gspa-compact-dense
#SBATCH --partition=debug
#SBATCH --nodelist=node002
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=08:00:00
#SBATCH --output=/data/hohndor/gspa-v3/compact-audit-v1/slurm/dense-%j.log
set -euo pipefail
root=/data/hohndor/gspa-v3/compact-audit-v1
export TMPDIR=/scratch/$SLURM_JOB_ID
mkdir -p "$TMPDIR/source"
export JAVA_TOOL_OPTIONS="-Djava.io.tmpdir=$TMPDIR"
unset JDK_JAVA_OPTIONS _JAVA_OPTIONS
tar -xf "$root/source-v2.tar.gz" -C "$TMPDIR/source"
source=$TMPDIR/source/home/leechuck/Public/software/gspa
/data/hohndor/gspa-v3/price-production-v1/venv/bin/python "$root/dense_parity.py" "$source" "$root/jobs/10275/gspa-3.0.0.jar" "$root/dense/job-$SLURM_JOB_ID" "$TMPDIR"
