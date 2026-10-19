#!/usr/bin/env bash
# Run a P10 step on node007 against the /scratch/hohndor_p10 mirror of /data (see node007.sh).
#   sbatch -w node007 -c <cpus> --mem=<mem> -o /scratch/hohndor_p10/logs/%x_%j.out node007_job.sbatch <step> [args]
# steps: runlist <runlist.tsv> <n_par>   |   py <python args...>   |   sh <script> [args]
#SBATCH -p debug
#SBATCH -t 24:00:00
set -uo pipefail
S7=/scratch/hohndor_p10
export JAVA=$S7/data/hohndor/envs/jdk-21.0.10+7/bin/java
export JAR=$S7/data/hohndor/gspa/bin/gspa-unified.jar
export REF=$S7/data/hohndor/gspa/reference
export P10_PIPE=$S7/data/hohndor/gspa/proteomes/pillar4_repro_2026-10-19/code
export XMX=${XMX:-6g}
echo "host=$(hostname) job=${SLURM_JOB_ID:-none} step=$1 start=$(date -Is)"
step=$1; shift
case "$step" in
  runlist) grep -v '^run_id' "$1" | xargs -P "$2" -L 1 "$P10_PIPE/integrate_one.sh" ;;
  py) python3 "$@" ;;
  sh) bash "$@" ;;
esac
echo "end=$(date -Is)"
