#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=gspa-resume
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --no-requeue
#SBATCH --output=/scratch/project/k10018/gspa-v3/out/source-context-resume-v2/slurm/%j.log
# Resume of the 10406 production replay on Shaheen III, in two phases driven by MODE:
#   MODE=arm-pack  one exclusive node; every arm of $ARMS_FILE ("source/context/kind" per line) runs as a concurrent
#                  process (PARALLEL at a time): run_source_context_remote.py resume --replay-arm <arm>. All share one
#                  restored results tree (restored once, guarded by restore.lock / restored.json; Task M2) and each
#                  writes only its own arm. The shared QoS caps a job at 8 G, below an arm's JVM heap, hence workq.
#   MODE=final     one process, after every arm exited 0: plain resume → reuses every audited arm, runs the analysis
#                  stages, writes replay/verification.json and the analysis.
# Container, package and binds are the unimatrix recipe (source_context_comparison.sbatch, resume branch) with the
# Shaheen host mapping of rolefence-arm.sbatch; seal.py runs once per bundle inside the container, where its /data
# paths resolve through the read-only mirror bind. TEST_GATES (JSON pattern→count) pins the test counts.
set -euo pipefail
: "${MODE:?arm-pack or final}"
module load singularity
export script_self=$(readlink -f "$0")   # $0 is not the script inside xargs/bash -c workers
export image=/project/k10018/sw/images/debian-trixie data_root=/project/k10018/gspa-v3 logical_data=/data/hohndor/gspa-v3
export bundle=${BUNDLE:-source-context-resume-v2}
export out_root=/scratch/project/k10018/gspa-v3/out/$bundle root=$logical_data/$bundle
export shared=$out_root/results workbase=/scratch/project/k10018/gspa-v3/work/$SLURM_JOB_ID MODE KEEP_WORK=${KEEP_WORK:-0}
export logical_source=/home/leechuck/tmp/source-context-current logical_results=/home/leechuck/tmp/source-context-results
export jar=$logical_data/ec-prefix-identity-v3/jobs/10402/gspa-3.0.0.jar
export probe=$logical_data/multisource-comparison-probe-v4/jobs/10403/results/verification.json
export training=$logical_data/source-context-full-v1/jobs/10406/checkpoint/training
export replay=$logical_data/source-context-full-v1/jobs/10406/results/replay
export TEST_GATES=${TEST_GATES:-{}}
mkdir -p "$out_root/slurm" "$out_root/jobs" "$shared" "$workbase/fs" "$workbase/tmp"
rsync -a --ignore-existing --exclude jobs/ --exclude slurm/ --exclude results/ "$data_root/$bundle/" "$out_root/"
# One package + source extraction per node job (scratch inode quota: 1 M files per user); arms mount it
# read-only with their own writable /home/leechuck/tmp on top.
zstd -dc "$data_root/observed-gap-context-modes-v1-relocation/package.tar.zst" | tar -xf - -C "$workbase/fs"
tar -xf "$data_root/context-resume-v1/supplement.tar" -C "$workbase/fs"
tar -xf "$out_root/parent-supplement.tar" -C "$workbase/fs"
mkdir -p "$workbase/fs$logical_source" "$workbase/fs/home/leechuck/tmp/arm-tmp"
tar -xzf "$out_root/source.tar.gz" -C "$workbase/fs$logical_source"
( cd "$workbase/fs$logical_source" && sha256sum -c "$out_root/source.sha256" > "$out_root/jobs/$SLURM_JOB_ID-source-check.log" )
run_one() {   # $1 = task label, $2 = arm id or "final"
  local task=$1 arm=$2 work=$workbase/$1 jobout=$out_root/jobs/$SLURM_JOB_ID-$1 code=0
  mkdir -p "$work/tmp" "$work/home-tmp/jobout" "$jobout"
  date --iso-8601=seconds > "$jobout/started-at.txt"; hostname > "$jobout/hostname.txt"; cp "$script_self" "$jobout/executed.sbatch"
  echo "$arm" > "$jobout/arm.txt"
  (
    set -euo pipefail; export TMPDIR=$work/tmp
    singularity --version > "$jobout/singularity-version.txt"
    cp "$out_root/jobs/$SLURM_JOB_ID-source-check.log" "$jobout/source-check.log"
    container=(singularity exec --no-home --no-mount tmp --pwd "$logical_source"
     --bind /project/k10018/sw/envs/gspa-tools:/project/k10018/sw/envs/gspa-tools:ro
     --bind "$data_root:$logical_data:ro" --bind "$out_root:$root"
     --bind "$workbase/fs/home/leechuck:/home/leechuck:ro" --bind "$work/home-tmp:/home/leechuck/tmp/arm-tmp" --bind "$jobout:/home/leechuck/tmp/arm-tmp/jobout" --bind "$shared:$logical_results"
     --bind "$workbase/fs/usr/lib/jvm/java-21-openjdk-amd64:/usr/lib/jvm/java-21-openjdk-amd64:ro"
     --env TMPDIR=/home/leechuck/tmp/arm-tmp --env JAVA_TOOL_OPTIONS=-Djava.io.tmpdir=/home/leechuck/tmp/arm-tmp --env PYTHONDONTWRITEBYTECODE=1
     --env PATH=/usr/lib/jvm/java-21-openjdk-amd64/bin:/usr/bin:/bin:/project/k10018/sw/envs/gspa-tools/bin --env LC_ALL=C.UTF-8 --env XDG_CONFIG_HOME=/home/leechuck/.config --env XDG_CACHE_HOME=/home/leechuck/.cache --env XDG_DATA_HOME=/home/leechuck/.local/share
     --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1
     --env PYTHONPATH="$logical_source/vendor" --env IJSON_BACKEND=python --env TEST_GATES="$TEST_GATES"
     "$image")
    printf '%s\n' "${container[@]}" > "$jobout/container-command.txt"
    python=/home/leechuck/tmp/gspa-v3-testenv/bin/python
    "${container[@]}" java -version > "$jobout/java-version.txt" 2>&1
    "${container[@]}" rsync --version > "$jobout/rsync-version.txt" 2>&1
    # Seal once per bundle (seal.py creates its products exclusively); other processes wait for the sealer.
    if [ ! -f "$out_root/packages.sha256" ]; then
      if mkdir "$out_root/seal.lock" 2>/dev/null; then
        "${container[@]}" "$python" "$root/seal.py" "$root" > "$jobout/seal.log" 2>&1 || { rm -rf "$out_root/seal.lock"; exit 3; }
        rm -rf "$out_root/seal.lock"
      else
        for i in $(seq 1 180); do [ -f "$out_root/packages.sha256" ] && break; sleep 10; done
        [ -f "$out_root/packages.sha256" ] || { echo "seal did not appear" >&2; exit 3; }
      fi
    fi
    cp "$out_root/packages.sha256" "$out_root/source.sha256" "$jobout/"
    "${container[@]}" sha256sum -c "$root/packages.sha256" > "$jobout/package-check.log"
    "${container[@]}" "$python" - <<'PY' > "$jobout/tests.log" 2>&1
import unittest, json, os
expected = json.loads(os.environ.get('TEST_GATES', '{}'))
suite = unittest.TestSuite()
for pattern in ('test_source_arm_context*.py', 'test_verified_batch.py', 'test_streaming_ledger.py', 'test_ec_reference_panel.py',
                'test_observed_gap_context_analysis.py', 'test_source_context_remote.py', 'test_cohort_sequence_inputs.py'):
    tests = unittest.defaultTestLoader.discover('benchmark/v3/tests', pattern=pattern)
    if pattern in expected:
        assert tests.countTestCases() == expected[pattern], (pattern, tests.countTestCases(), expected[pattern])
    suite.addTests(tests)
result = unittest.TextTestRunner(verbosity=2).run(suite)
assert result.wasSuccessful() and not result.skipped
PY
    args=(resume --out "$logical_results" --logs /home/leechuck/tmp/arm-tmp/jobout/execution --bindings "$root/packages.sha256" --jar "$jar"
     --clean-root "$logical_data/clean-cohort-v3/collection" --clean-audit "$logical_data/clean-cohort-v3/jobs/10368/audit.json"
     --probe "$probe" --resume-training "$training" --resume-replay "$replay"
     --resume-changed-source benchmark/v3/check_source_arm_context_models.py)   # the only patched *probed* source (10403 probe); run_source_context_remote.py is not probed
    [ "$arm" = final ] || args+=(--replay-arm "$arm")
    # set -e is inert inside a subshell that is part of a || list: propagate the runner's status explicitly.
    /usr/bin/time -v -o "$jobout/resources.txt" "${container[@]}" "$python" "$logical_source/benchmark/v3/run_source_context_remote.py" "${args[@]}" > "$jobout/execution.log" 2>&1 || exit $?
    if [ "$arm" = final ]; then rsync -a --exclude 'training/' "$shared/" "$jobout/results/" || exit $?; fi
  ) || code=$?
  printf '%s\n' "$code" > "$jobout/exit-code.txt"; date --iso-8601=seconds > "$jobout/finished-at.txt"
  [ "$code" -eq 0 ] && [ "$KEEP_WORK" != 1 ] && rm -rf "$work"
  echo "task $task arm=$arm exit=$code"
}
export -f run_one
if [ "$MODE" = final ]; then
  run_one final final; exit "$(cat "$out_root/jobs/$SLURM_JOB_ID-final/exit-code.txt")"
fi
: "${ARMS_FILE:?file with one source/context/kind per line}"
mapfile -t arms < "$ARMS_FILE"
echo "job $SLURM_JOB_ID on $(hostname): ${#arms[@]} arms, parallel ${PARALLEL:-12}"
for i in "${!arms[@]}"; do printf '%s\t%s\n' "arm$i" "${arms[$i]}"; done | xargs -P "${PARALLEL:-12}" -L1 bash -c 'run_one "$0" "$1"'
fails=0; for i in "${!arms[@]}"; do c=$(cat "$out_root/jobs/$SLURM_JOB_ID-arm$i/exit-code.txt" 2>/dev/null || echo missing); [ "$c" = 0 ] || { echo "FAILED arm$i ${arms[$i]} exit=$c"; fails=$((fails+1)); }; done
[ "$KEEP_WORK" = 1 ] || rm -rf "$workbase/fs" "$workbase/tmp"; rmdir "$workbase" 2>/dev/null || true
echo "done, failures=$fails"; [ "$fails" -eq 0 ]
