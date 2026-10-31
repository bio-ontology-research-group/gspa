#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=shared
#SBATCH --job-name=gspa-wp5-recovery
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --no-requeue
#SBATCH --output=/scratch/project/k10018/gspa-v3/out/wp5-recovery-v1/slurm/%j.log
# WP5 masked assignment recovery on Shaheen III (design tasks/v3-wp5-masked-recovery-design-2026-10-29.md §8-§9).
# Prepared by Task N; NOT submitted. Container conventions as rolefence-arm.sbatch: the unimatrix data root
# /data/hohndor/gspa-v3 is the read-only mirror /project/k10018/gspa-v3 bound at its original path; the writable
# output root is bound at the logical path of the run; the Python test env and the JDK come from the relocated
# package (observed-gap-context-modes-v1-relocation), rsync/zstd/DIAMOND from the gspa-tools conda env.
#
# Inputs (read-only mirror, staged beforehand under /project/k10018/gspa-v3/wp5-recovery-v1/):
#   source.tar.gz + source.sha256          source snapshot containing the WP5 scripts and scripts/*experiments.py
#   inputs/...                             the same layout as benchmark/v3/wp5_recovery.sbatch (Keio cases,
#                                          nutrient and sulfonate inputs/observations, iML1515/iYO844, NCBI CDS
#                                          proteomes Ec.faa/Bs.faa, uniprot_sprot.dat.gz, gspa-3.0.0.jar)
#   SHA256SUMS                             hashes of source.tar.gz and every inputs/ file
# Outputs: /scratch/project/k10018/gspa-v3/out/wp5-recovery-v1/{bundle,runs,oracle,jobs/<jobid>,slurm}.
#
# Stages (STAGES, default below; STAGES=none stops after the WP5 unit tests): prepare reference query diamond
# claims construct audit protocol exhaustive deterministic cobra. There is no Qwen stage here: the Qwen endpoint
# (GSPA_LLM_URL, unimatrix01) is not reachable from Shaheen compute nodes, so the qwen and transcripts stages run
# on unimatrix (benchmark/v3/wp5_recovery.sbatch) against the same frozen bundle/protocol, and the COBRA replay
# here covers the deterministic arms only (--arm-stages deterministic). The audit is a gate: a failed masking audit
# stops the job before any LP. Evaluation (join/calibrate/evaluate/report) is not part of the job.
set -euo pipefail
module load singularity
image=/project/k10018/sw/images/debian-trixie
data_root=/project/k10018/gspa-v3
logical_data=/data/hohndor/gspa-v3
input_root=$data_root/wp5-recovery-v1
out_root=/scratch/project/k10018/gspa-v3/out/wp5-recovery-v1
root=$logical_data/wp5-recovery-v1-out
work=/scratch/project/k10018/gspa-v3/work/$SLURM_JOB_ID/wp5
jobout=$out_root/jobs/$SLURM_JOB_ID
logical_source=/home/leechuck/tmp/wp5-source
stages=${STAGES:-prepare reference query diamond claims construct audit protocol exhaustive deterministic cobra}
cobra_rest=${COBRA_REST:-0}
mkdir -p "$work/fs" "$jobout" "$out_root/slurm"
export TMPDIR=/scratch/project/k10018/gspa-v3/work/$SLURM_JOB_ID/tmp; mkdir -p "$TMPDIR"
finish() {
 code=$?; trap - EXIT; set +e
 printf '%s\n' "$code" > "$jobout/exit-code.txt"; date --iso-8601=seconds > "$jobout/finished-at.txt"
 if [ "$code" -eq 0 ] && [ "${KEEP_WORK:-0}" != 1 ]; then rm -rf "/scratch/project/k10018/gspa-v3/work/$SLURM_JOB_ID"; fi
 exit "$code"
}
trap 'exit 143' TERM INT   # scancel: let the EXIT trap record 143, not the last command's status
trap finish EXIT
date --iso-8601=seconds > "$jobout/started-at.txt"; hostname > "$jobout/hostname.txt"; cp "$0" "$jobout/executed.sbatch"
singularity --version > "$jobout/singularity-version.txt"
( cd "$input_root" && sha256sum -c SHA256SUMS ) > "$jobout/input-verification.log"
# Python test env (numpy, scipy, python-libsbml, cobra/GLPK) and JDK 21 from the relocated package, as rolefence-arm.
zstd -dc "$data_root/observed-gap-context-modes-v1-relocation/package.tar.zst" | tar -xf - -C "$work/fs"
mkdir -p "$work/fs$logical_source" "$work/fs/home/leechuck/tmp"
tar -xzf "$input_root/source.tar.gz" -C "$work/fs$logical_source"
( cd "$work/fs$logical_source" && sha256sum -c "$input_root/source.sha256" > "$jobout/source-check.log" )
container=(singularity exec --no-home --no-mount tmp --pwd "$root"
 --bind /project/k10018/sw/envs/gspa-tools:/project/k10018/sw/envs/gspa-tools:ro
 --bind "$data_root:$logical_data:ro" --bind "$out_root:$root"
 --bind "$work/fs/home/leechuck:/home/leechuck"
 --bind "$work/fs/usr/lib/jvm/java-21-openjdk-amd64:/usr/lib/jvm/java-21-openjdk-amd64"
 --env TMPDIR=/home/leechuck/tmp --env JAVA_TOOL_OPTIONS=-Djava.io.tmpdir=/home/leechuck/tmp
 --env PATH=/usr/lib/jvm/java-21-openjdk-amd64/bin:/usr/bin:/bin:/project/k10018/sw/envs/gspa-tools/bin --env LC_ALL=C.UTF-8 --env HOME=/home/leechuck --env XDG_CONFIG_HOME=/home/leechuck/.config --env XDG_CACHE_HOME=/home/leechuck/.cache --env XDG_DATA_HOME=/home/leechuck/.local/share
 --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1
 --env PYTHONPATH="$logical_source/benchmark/v3"
 "$image")
printf '%s\n' "${container[@]}" > "$jobout/container-command.txt"
python=/home/leechuck/tmp/gspa-v3-testenv/bin/python
src=$logical_source/benchmark/v3
inputs=$logical_data/wp5-recovery-v1/inputs
"${container[@]}" java -version > "$jobout/java-version.txt" 2>&1
"${container[@]}" "$python" -c 'import sys,numpy,scipy,libsbml,cobra; print(sys.version); print(numpy.__version__, scipy.__version__, libsbml.getLibSBMLDottedVersion(), cobra.__version__)' > "$jobout/python-version.txt" 2>&1
"${container[@]}" "$python" -m unittest discover -s "$logical_source/benchmark/v3/tests" -p 'test_wp5_*.py' > "$jobout/tests.log" 2>&1
[ "$stages" = none ] && { echo smoke-ok > "$jobout/smoke.txt"; exit 0; }

run() {  # run <stage> <args...>: one containerized stage with resources and a log
  local stage=$1; shift
  date --iso-8601=seconds > "$jobout/$stage.started"
  /usr/bin/time -v -o "$jobout/$stage.resources.txt" "${container[@]}" "$@" > "$jobout/$stage.log" 2>&1
  date --iso-8601=seconds > "$jobout/$stage.finished"
}
shards() {  # shards <stage>: SLURM_CPUS_PER_TASK concurrent run_wp5_recovery.py shards
  local stage=$1 n=$SLURM_CPUS_PER_TASK pids=() status=0
  date --iso-8601=seconds > "$jobout/$stage.started"
  for i in $(seq 0 $((n - 1))); do
    "${container[@]}" "$python" "$src/run_wp5_recovery.py" "$stage" --public "$root/bundle/public" --out "$root/runs" \
      --shard "$i/$n" > "$jobout/$stage.shard$i.log" 2>&1 &
    pids+=($!)
  done
  for pid in "${pids[@]}"; do wait "$pid" || status=1; done
  date --iso-8601=seconds > "$jobout/$stage.finished"
  return $status
}
for stage in $stages; do
  case "$stage" in
  prepare)
    [ -e "$out_root/bundle" ] && { echo "bundle exists in $out_root; refusing to overwrite" >&2; exit 2; }
    run prepare "$python" "$src/prepare_wp5_cases.py" --out "$root/bundle" \
      --keio-cases "$inputs/keio/cases_manuscript_regulatory.jsonl" --keio-model "$inputs/models/iML1515.xml" \
      --ec-model "$inputs/models/iML1515.xml" --bs-model "$inputs/models/iYO844.xml" \
      --nutrient-inputs "$inputs/nutrient/inputs.jsonl.gz" --nutrient-observations "$inputs/nutrient/observations.jsonl.gz" \
      --sulfonate-inputs "$inputs/sulfonate/inputs.jsonl" --sulfonate-observations "$inputs/sulfonate/observations.jsonl" ;;
  reference)
    run reference "$python" "$src/build_wp5_candidates.py" reference --bundle "$root/bundle" \
      --swissprot "$inputs/swissprot/uniprot_sprot.dat.gz" ;;
  query)
    run query "$python" "$src/build_wp5_candidates.py" query --bundle "$root/bundle" \
      --proteome "Ec=$inputs/proteomes/Ec.faa" --proteome "Bs=$inputs/proteomes/Bs.faa" ;;
  diamond)
    "${container[@]}" sh -c 'command -v diamond' > "$jobout/diamond-path.txt" || { echo "diamond not in the container PATH" >&2; exit 2; }
    run diamond "$python" "$src/build_wp5_candidates.py" diamond --bundle "$root/bundle" --threads "$SLURM_CPUS_PER_TASK" \
      --diamond-bin diamond ;;
  claims)
    run claims "$python" "$src/build_wp5_candidates.py" claims --bundle "$root/bundle" ;;
  construct)
    run construct "$python" "$src/build_wp5_candidates.py" construct --bundle "$root/bundle" --jar "$inputs/gspa-3.0.0.jar" \
      --java java --python "$python" --sidecar "$logical_source/scripts/metabolic_experiments.py" ;;
  audit)
    run audit "$python" "$src/audit_wp5_masking.py" all --bundle "$root/bundle" --jar "$inputs/gspa-3.0.0.jar" --java java \
      --python "$python" --sidecar "$logical_source/scripts/metabolic_experiments.py" --scratch /home/leechuck/tmp/wp5-audit ;;
  protocol)
    run protocol "$python" "$src/run_wp5_recovery.py" protocol --public "$root/bundle/public" --out "$root/runs" ;;
  exhaustive|deterministic)
    shards "$stage" ;;
  cobra)
    mkdir -p "$out_root/oracle"
    run cobra "$python" "$src/check_wp5_cobra.py" --public "$root/bundle/public" --runs "$root/runs" \
      --out "$root/oracle/cobra-deterministic.json" --arm-stages deterministic --rest-sample "$cobra_rest" ;;
  qwen|transcripts)
    echo "stage $stage needs the unimatrix01 Qwen endpoint; run it with benchmark/v3/wp5_recovery.sbatch" >&2; exit 2 ;;
  *)
    echo "unknown stage: $stage" >&2; exit 2 ;;
  esac
done
echo done > "$jobout/done.txt"
