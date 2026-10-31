#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=gspa-wp5-cobra-shard
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --no-requeue
#SBATCH --exclude=nid00022
#SBATCH --output=/scratch/project/k10018/gspa-v3/out/wp5-recovery-v1/slurm/%j.log
# WP5 COBRA replay fanned out: shards FIRST..FIRST+COUNT-1 of SHARDS (check_wp5_cobra.py --shard i/SHARDS, from the
# same frozen source.tar.gz and container recipe as wp5_recovery.sbatch, deterministic arms only, rest sample 0) run
# concurrently on one exclusive node. Reports go to oracle/cobra-deterministic.shards/shard-<i>-of-<SHARDS>.json;
# merge_wp5_cobra.py rebuilds the unsharded report from them. The Python env and JDK come from the shared read-only
# package tree used by arm-pack.sbatch (scratch inode quota), the source snapshot is extracted per job.
set -euo pipefail
: "${SHARDS:?}" "${FIRST:?}" "${COUNT:?}" "${ARM_STAGE:=deterministic}"
module load singularity
image=/project/k10018/sw/images/debian-trixie
data_root=/project/k10018/gspa-v3
logical_data=/data/hohndor/gspa-v3
input_root=$data_root/wp5-recovery-v1
out_root=/scratch/project/k10018/gspa-v3/out/wp5-recovery-v1
root=$logical_data/wp5-recovery-v1-out
work=/scratch/project/k10018/gspa-v3/work/$SLURM_JOB_ID
jobout=$out_root/jobs/$SLURM_JOB_ID
logical_source=/home/leechuck/tmp/arm-tmp/wp5-source
fsroot=$(ls -d /scratch/project/k10018/gspa-v3/shared/arm-fs-*/ | head -1); fsroot=${fsroot%/}
[ -f "$fsroot/.complete" ] || { echo "shared package tree missing: $fsroot" >&2; exit 3; }
shard_dir=oracle/cobra-$ARM_STAGE.shards
mkdir -p "$work/home-tmp" "$jobout" "$out_root/$shard_dir"
finish() {
 code=$?; trap - EXIT; set +e
 printf '%s\n' "$code" > "$jobout/exit-code.txt"; date --iso-8601=seconds > "$jobout/finished-at.txt"
 [ "${KEEP_WORK:-0}" = 1 ] || rm -rf "$work"
 exit "$code"
}
trap 'exit 143' TERM INT
trap finish EXIT
date --iso-8601=seconds > "$jobout/started-at.txt"; hostname > "$jobout/hostname.txt"; cp "$0" "$jobout/executed.sbatch"
echo "shards $FIRST..$((FIRST+COUNT-1)) of $SHARDS fsroot=$fsroot" > "$jobout/shards.txt"
( cd "$input_root" && sha256sum -c SHA256SUMS ) > "$jobout/input-verification.log"
mkdir -p "$work/home-tmp/wp5-source"
tar -xzf "$input_root/source.tar.gz" -C "$work/home-tmp/wp5-source"
( cd "$work/home-tmp/wp5-source" && sha256sum -c "$input_root/source.sha256" > "$jobout/source-check.log" )
container=(singularity exec --no-home --no-mount tmp --pwd "$root"
 --bind "$data_root:$logical_data:ro" --bind "$out_root:$root"
 --bind "$fsroot/home/leechuck:/home/leechuck:ro" --bind "$work/home-tmp:/home/leechuck/tmp/arm-tmp"
 --env TMPDIR=/home/leechuck/tmp/arm-tmp --env LC_ALL=C.UTF-8 --env PYTHONDONTWRITEBYTECODE=1 --env XDG_CACHE_HOME=/home/leechuck/tmp/arm-tmp/.cache --env XDG_CONFIG_HOME=/home/leechuck/tmp/arm-tmp/.config --env XDG_DATA_HOME=/home/leechuck/tmp/arm-tmp/.local/share --env MPLCONFIGDIR=/home/leechuck/tmp/arm-tmp/.mpl
 --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1
 --env PYTHONPATH="$logical_source/benchmark/v3"
 "$image")
printf '%s\n' "${container[@]}" > "$jobout/container-command.txt"
python=/home/leechuck/tmp/gspa-v3-testenv/bin/python
src=$logical_source/benchmark/v3
"${container[@]}" "$python" -c 'import sys,cobra; print(sys.version); print(cobra.__version__)' > "$jobout/python-version.txt" 2>&1
pids=() status=0
for i in $(seq "$FIRST" $((FIRST + COUNT - 1))); do
  [ "$i" -lt "$SHARDS" ] || break
  out=$shard_dir/shard-$i-of-$SHARDS.json
  [ -s "$out_root/$out" ] && { echo "shard $i exists, skipped"; continue; }
  ( c=0; /usr/bin/time -v -o "$jobout/shard$i.resources.txt" "${container[@]}" "$python" "$src/check_wp5_cobra.py" \
      --public "$root/bundle/public" --runs "$root/runs" --out "$root/$out.partial" --arm-stages "$ARM_STAGE" \
      --shard "$i/$SHARDS" > "$jobout/shard$i.log" 2>&1 || c=$?
    [ -s "$out_root/$out.partial" ] && mv "$out_root/$out.partial" "$out_root/$out"; exit $c ) &
  pids+=($!)
done
# check_wp5_cobra.py exits 1 when an oracle disagrees: that is a result, kept in the shard report, not a job failure.
for pid in "${pids[@]}"; do wait "$pid" || { c=$?; [ "$c" -eq 1 ] || status=1; }; done
missing=0; for i in $(seq "$FIRST" $((FIRST + COUNT - 1))); do [ "$i" -lt "$SHARDS" ] || break; [ -s "$out_root/$shard_dir/shard-$i-of-$SHARDS.json" ] || missing=$((missing+1)); done
echo "missing=$missing status=$status"; [ "$missing" -eq 0 ] && [ "$status" -eq 0 ]
