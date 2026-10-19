#!/usr/bin/env bash
# P03 reference panels (Fig. 2d, Suppl. Table S2, former Table 2 external block).
# Runs on unimatrix01 via Slurm (node007, where /data is seen as /srv/data),
# then aggregates locally into results/*.tsv.
#   ./run.sh            sync code, submit both jobs, wait, fetch, aggregate
#   ./run.sh aggregate  only fetch + aggregate (remote outputs already exist)
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REMOTE=unimatrix01
R=/data/hohndor/gspa/proteomes/bench10/repro_2026-10-19
FETCH="${P03_FETCH:-$HERE/fetched}"
mode="${1:-all}"

if [[ "$mode" == all ]]; then
  ssh $REMOTE "mkdir -p $R/code ~/p03logs /data/hohndor/gspa/proteomes/bench9_unified/repro_2026-10-19/panel"
  scp -q "$HERE"/p03_*.py "$HERE"/p03_panel13.sbatch "$HERE"/p03_panel21.sbatch $REMOTE:$R/code/
  j13=$(ssh $REMOTE "cd $R/code && sbatch --parsable p03_panel13.sbatch")
  j21=$(ssh $REMOTE "cd $R/code && sbatch --parsable p03_panel21.sbatch")
  echo "submitted panel13=$j13 panel21=$j21"
  ssh $REMOTE "while squeue -h -j $j13,$j21 2>/dev/null | grep -q .; do sleep 60; done"
  ssh $REMOTE "grep -h '^STEP\|rc=' ~/p03logs/p03_panel13_$j13.out ~/p03logs/p03_panel21_$j21.out"
fi

mkdir -p "$FETCH"
scp -q $REMOTE:$R/eval/panel13_summary.json $REMOTE:$R/eval/panel13_per_genome.tsv \
       $REMOTE:$R/audit/panel13_audit_t05.tsv $REMOTE:$R/audit/panel13_audit_fmaxt.tsv $REMOTE:$R/perturbation/panel13_perturbation.tsv \
       $REMOTE:$R/panel21/eval/panel21_summary.json $REMOTE:$R/panel21/eval/panel21_per_genome.tsv "$FETCH/"
python3 "$HERE/p03_aggregate.py"
