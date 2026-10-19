#!/usr/bin/env bash
# Run every line of a run-list TSV (run_id taxon profile claims operons gaps outdir)
# with integrate_one.sh, N_PAR at a time. Usage: sbatch run_list.sbatch <runlist.tsv> [n_par]
#SBATCH -p debug
#SBATCH --exclude=unimatrix,node007
#SBATCH -t 12:00:00
set -uo pipefail
LIST=$1; NPAR=${2:-2}
export P10_PIPE=${P10_PIPE:-/data/hohndor/gspa/proteomes/pillar4_repro_2026-10-19/code}
echo "host=$(hostname) job=${SLURM_JOB_ID:-none} list=$LIST npar=$NPAR start=$(date -Is)"
grep -v '^run_id' "$LIST" | xargs -P "$NPAR" -L 1 "$P10_PIPE/integrate_one.sh"
echo "end=$(date -Is)"
