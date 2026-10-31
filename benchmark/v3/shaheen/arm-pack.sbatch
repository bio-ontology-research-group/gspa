#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=gspa-arm-pack
#SBATCH --nodes=1
#SBATCH --time=08:00:00
#SBATCH --no-requeue
#SBATCH --output=/scratch/project/k10018/gspa-v3/out/arm-logs/pack-%j.log
# Shaheen III: run many arms of an arms.tsv concurrently on one exclusive workq node (the shared QoS caps a job
# at 8 CPUs = 8 G, below the 12 G JVM heap of an arm). Same per-arm recipe as arm-array.sbatch (unimatrix jobs
# 10666 / 10686): the recorded command JSON runs unchanged inside the container with its own package extraction.
# Environment: BUNDLE, ARMS (arms.tsv relative to the bundle root), ARM_COLS (2 = ctx,kind; 3 = regime,ctx,kind),
# TASKS (space-separated 0-based row indices, default: all rows), PARALLEL (concurrent arms, default 12), KEEP_WORK,
# EXTRA_EXCLUDES (extra rsync --exclude options for the one-time bundle copy, e.g. large unrelated plan trees).
set -euo pipefail
: "${BUNDLE:?}" "${ARMS:?}" "${ARM_COLS:?}"
: "${NAME_COLS:=}" "${MANIFEST_COL:=}" "${COMMAND_COL:=}" "${VERIFIED:=0}" "${INFERENCE_COL:=}"
module load singularity
export script_self=$(readlink -f "$0")   # $0 is not the script inside xargs/bash -c workers
export image=/project/k10018/sw/images/debian-trixie data_root=/project/k10018/gspa-v3 logical_data=/data/hohndor/gspa-v3
export out_root=/scratch/project/k10018/gspa-v3/out/$BUNDLE root=$logical_data/$BUNDLE BUNDLE ARMS ARM_COLS NAME_COLS MANIFEST_COL COMMAND_COL VERIFIED INFERENCE_COL
export workbase=/scratch/project/k10018/gspa-v3/work/$SLURM_JOB_ID KEEP_WORK=${KEEP_WORK:-0}
# One read-only package tree shared by every arm job (scratch inode quota: 1,024,000 files per user; one
# extraction is ~36,000 files, so 26 concurrent per-job extractions exhausted it on 2026-10-29). Built once under a
# lock and identified by the sha256 of its three archives; jobs only read it.
pkg_sha=$(cat "$data_root/observed-gap-context-modes-v1-relocation/package.tar.zst" "$data_root/context-resume-v1/supplement.tar" \
  "$data_root/source-context-full-v1/parent-supplement.tar" | sha256sum | cut -c1-16)
export fsroot=/scratch/project/k10018/gspa-v3/shared/arm-fs-$pkg_sha
mkdir -p "$out_root/runs" "$workbase/tmp" /scratch/project/k10018/gspa-v3/shared
# One package extraction per node job (the scratch inode quota is 1 M files per user; a per-arm copy of the
# 28 GB / ~90 k-file package exhausted it on 2026-10-29). Arms mount it read-only; each arm gets its own
# writable /home/leechuck/tmp (JVM and Python temp, TMPDIR) on top of it.
if [ ! -f "$fsroot/.complete" ]; then
  if mkdir "$fsroot.lock" 2>/dev/null; then
    rm -rf "$fsroot.partial"; mkdir -p "$fsroot.partial"
    zstd -dc "$data_root/observed-gap-context-modes-v1-relocation/package.tar.zst" | tar -xf - -C "$fsroot.partial"
    tar -xf "$data_root/context-resume-v1/supplement.tar" -C "$fsroot.partial"
    tar -xf "$data_root/source-context-full-v1/parent-supplement.tar" -C "$fsroot.partial"
    mkdir -p "$fsroot.partial/home/leechuck/tmp/arm-tmp"
    mv "$fsroot.partial" "$fsroot" && touch "$fsroot/.complete"; rmdir "$fsroot.lock"
  else
    for i in $(seq 1 180); do [ -f "$fsroot/.complete" ] && break; sleep 10; done
    [ -f "$fsroot/.complete" ] || { echo "shared package tree did not appear: $fsroot" >&2; exit 3; }
  fi
fi
rsync -a --exclude runs/ --exclude slurm/ --exclude jobs/ --exclude 'runs-*/' --exclude 'mmseqs/tmp/' ${EXTRA_EXCLUDES:-} "$data_root/$BUNDLE/" "$out_root/"
nrows=$(( $(wc -l < "$out_root/$ARMS") - 1 ))
tasks=${TASKS:-$(seq -s ' ' 0 $((nrows-1)))}
echo "job $SLURM_JOB_ID on $(hostname): bundle=$BUNDLE arms=$ARMS rows=$nrows tasks=[$tasks] parallel=${PARALLEL:-12}"
run_arm() {
  local i=$1 row name manifest command arm work
  row=$(awk -v i=$((i+1)) -F'\t' 'NR>1{n++} n==i{print; exit}' "$out_root/$ARMS")
  [ -n "$row" ] || { echo "no arm row for task $i" >&2; return 2; }
  # Column layout: default = ARM_COLS leading name columns, then manifest, then command (Task L / role-genome
  # tables); the Task R schema (block,target,regime,ctx,kind,variant,…,manifest,command,…) is selected with
  # NAME_COLS (comma list, e.g. 3,4,5), MANIFEST_COL and COMMAND_COL (1-based).
  name=$(cut -f"${NAME_COLS:-1-$ARM_COLS}" <<<"$row" | tr '\t' '/'); manifest=$(cut -f"${MANIFEST_COL:-$((ARM_COLS+1))}" <<<"$row"); command=$(cut -f"${COMMAND_COL:-$((ARM_COLS+2))}" <<<"$row")
  arm=$out_root/runs/$name; work=$workbase/arm-$i; mkdir -p "$arm" "$work/tmp" "$work/home-tmp"
  echo "task $i arm=$name manifest=$manifest command=$command host=$(hostname) job=$SLURM_JOB_ID" | tee "$arm/task.txt"
  date --iso-8601=seconds > "$arm/started-at.txt"; cp "$script_self" "$arm/executed.sbatch"
  local code=0
  (
    set -euo pipefail
    export TMPDIR=$work/tmp
    container=(singularity exec --no-home --no-mount tmp --pwd /home/leechuck
     --bind /project/k10018/sw/envs/gspa-tools:/project/k10018/sw/envs/gspa-tools:ro
     --bind "$data_root:$logical_data:ro" --bind "$out_root:$root"
     --bind "$fsroot/home/leechuck:/home/leechuck:ro" --bind "$work/home-tmp:/home/leechuck/tmp/arm-tmp"
     --bind "$fsroot/usr/lib/jvm/java-21-openjdk-amd64:/usr/lib/jvm/java-21-openjdk-amd64:ro"
     --env TMPDIR=/home/leechuck/tmp/arm-tmp --env JAVA_TOOL_OPTIONS=-Djava.io.tmpdir=/home/leechuck/tmp/arm-tmp --env PYTHONDONTWRITEBYTECODE=1
     --env PATH=/usr/lib/jvm/java-21-openjdk-amd64/bin:/usr/bin:/bin:/project/k10018/sw/envs/gspa-tools/bin --env LC_ALL=C.UTF-8 --env XDG_CONFIG_HOME=/home/leechuck/.config --env XDG_CACHE_HOME=/home/leechuck/.cache --env XDG_DATA_HOME=/home/leechuck/.local/share
     --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1
     --bind /project/k10018/code/gspa-v3:/opt/gspa:ro --env PYTHONPATH=/opt/gspa/benchmark/v3:/opt/gspa/vendor --env IJSON_BACKEND=python
     "$image")
    printf '%s\n' "${container[@]}" > "$arm/container-command.txt"
    "${container[@]}" java -version > "$arm/java-version.txt" 2>&1
    mapfile -t cmd < <(python3 -c 'import json,sys; [print(t) for t in json.load(open(sys.argv[1]))]' "${command/#$logical_data/$data_root}")
    printf '%s\n' "${cmd[@]}" > "$arm/command.txt"
    if [ "${VERIFIED:-0}" = 1 ]; then
      # Verified-batch contract (receipt verified_run.json in the inference dir; earlier outputs go to history/).
      : "${INFERENCE_COL:?VERIFIED=1 needs INFERENCE_COL}"
      inference=$(cut -f"$INFERENCE_COL" <<<"$row")
      /usr/bin/time -v -o "$arm/resources.txt" "${container[@]}" /home/leechuck/tmp/gspa-v3-testenv/bin/python \
        /opt/gspa/benchmark/v3/run_verified_arm.py --command-json "$command" --manifest "$manifest" \
        --out-dir "$inference" > "$arm/integrate.log" 2>&1
    else
      /usr/bin/time -v -o "$arm/resources.txt" "${container[@]}" "${cmd[@]}" > "$arm/integrate.log" 2>&1
    fi
  ) || code=$?
  printf '%s\n' "$code" > "$arm/exit-code.txt"; date --iso-8601=seconds > "$arm/finished-at.txt"
  [ "$code" -eq 0 ] && [ "$KEEP_WORK" != 1 ] && rm -rf "$work"
  echo "task $i arm=$name exit=$code"
  return 0
}
export -f run_arm
printf '%s\n' $tasks | xargs -P "${PARALLEL:-12}" -I{} bash -c 'run_arm {}'
fails=0; for i in $tasks; do row=$(awk -v i=$((i+1)) -F'\t' 'NR>1{n++} n==i{print; exit}' "$out_root/$ARMS"); name=$(cut -f"${NAME_COLS:-1-$ARM_COLS}" <<<"$row" | tr '\t' '/'); c=$(cat "$out_root/runs/$name/exit-code.txt" 2>/dev/null || echo missing); [ "$c" = 0 ] || { echo "FAILED task $i $name exit=$c"; fails=$((fails+1)); }; done
# Per-job work (arm scratch) is removed whatever the outcome unless KEEP_WORK=1; failed arms keep their logs in runs/.
[ "$KEEP_WORK" = 1 ] || rm -rf "$workbase"
echo "done, failures=$fails"; [ "$fails" -eq 0 ]
