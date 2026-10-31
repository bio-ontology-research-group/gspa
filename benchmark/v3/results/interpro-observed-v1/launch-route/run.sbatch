#!/bin/bash
#SBATCH --job-name=gspa-interpro-route
#SBATCH --partition=debug
#SBATCH --nodelist=node003
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --output=/data/hohndor/gspa-v3/interpro-route-v1/run-%j.log
set -euo pipefail
export TMPDIR="${SLURM_TMPDIR:-/scratch/$SLURM_JOB_ID}"
export JAVA_TOOL_OPTIONS="-Djava.io.tmpdir=$TMPDIR"
export PATH="/usr/lib/jvm/java-11-openjdk-amd64/bin:$PATH"
mkdir -p "$TMPDIR"
cd /data/hohndor/gspa-v3/interpro-route-v1
sha256sum -c launch.sha256
python3 check_interpro_route.py --jar gspa.jar --fixture /data/hohndor/gspa-v3/interpro-observed-v1/observed --tool /data/hohndor/envs/interproscan-5.72-103.0/interproscan.sh --java /usr/lib/jvm/java-21-openjdk-amd64/bin/java --out replay
