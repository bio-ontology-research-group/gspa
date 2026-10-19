#!/usr/bin/env bash
# node007 has no /data mount (2026-10-19); only its local /scratch and NFS /home.
# This helper mirrors /data paths under /scratch/hohndor_p10 on node007 by piping tar
# streams through srun (run it on the unimatrix01 login node; it only moves files).
#   bash node007.sh push <abs path> [<abs path> ...]   copy /data/... -> node007:/scratch/hohndor_p10/data/...
#   bash node007.sh pull <abs /data path>              copy node007 mirror of that path back to /data
#   bash node007.sh run  <cmd ...>                     run a short command on node007 (e.g. ls, tail)
set -euo pipefail
S7=/scratch/hohndor_p10
SR=(srun -p debug -w node007 -n 1 -c 1 --mem=1G -t 01:00:00 --quiet --chdir=/tmp)
mode=$1; shift
case "$mode" in
  push) tar -C / -cf - "${@#/}" | "${SR[@]}" bash -c "mkdir -p $S7 && tar -C $S7 -xf -" ;;
  pull) for p in "$@"; do mkdir -p "$p"; "${SR[@]}" tar -C "$S7$p" -cf - . | tar -C "$p" -xf -; done ;;
  run)  "${SR[@]}" bash -c "$*" ;;
  *) echo "usage: node007.sh push|pull|run ..."; exit 1 ;;
esac
