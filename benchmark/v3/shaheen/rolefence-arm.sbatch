#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=shared
#SBATCH --job-name=gspa-rolefence-arm
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=20G
#SBATCH --time=12:00:00
#SBATCH --no-requeue
#SBATCH --output=/scratch/project/k10018/gspa-v3/out/rolefence-arm-v1/slurm/%j.log
# WP2e on Shaheen III: the unimatrix recipe of rolefence-arm-v1/rolefence.sbatch (job 10673), unchanged inside the
# container. Only the host side differs: the unimatrix data root /data/hohndor/gspa-v3 is the read-only mirror
# /project/k10018/gspa-v3, bound at its original path, and the bundle root (training/, jobs/) is a writable copy
# on /scratch bound at the original path too. Compute nodes see /project read-only and have no $HOME.
# Environment: STAGES="prepare export train evaluate" (set STAGES=none for a smoke test that stops after the tests),
# KEEP_WORK=1 keeps the extracted package after a successful run.
set -euo pipefail
module load singularity
# Container root: a plain Debian 13 sandbox (the unimatrix recipe bound the Debian host's /usr,/lib,/bin,/etc into an
# empty root; Shaheen compute nodes run SLES without /usr/lib/jvm, so the JDK bind point cannot be created there).
# Java and Python still come from the package; rsync from the gspa-tools conda env.
image=/project/k10018/sw/images/debian-trixie
data_root=/project/k10018/gspa-v3
out_root=/scratch/project/k10018/gspa-v3/out/rolefence-arm-v1
logical_data=/data/hohndor/gspa-v3
root=$logical_data/rolefence-arm-v1
base=$logical_data/observed-gap-context-modes-v1-relocation
supplement=$logical_data/context-resume-v1
parent=$logical_data/source-context-full-v1
work=/scratch/project/k10018/gspa-v3/work/$SLURM_JOB_ID/rolefence
jobout=$out_root/jobs/$SLURM_JOB_ID
logical_source=/home/leechuck/tmp/source-context-current
stages=${STAGES:-prepare export train evaluate}
mkdir -p "$work/fs" "$jobout" "$out_root/slurm"
export TMPDIR=/scratch/project/k10018/gspa-v3/work/$SLURM_JOB_ID/tmp; mkdir -p "$TMPDIR"
finish() {
 code=$?; trap - EXIT; set +e
 printf '%s\n' "$code" > "$jobout/exit-code.txt"; date --iso-8601=seconds > "$jobout/finished-at.txt"
 if [ "$code" -eq 0 ] && [ "${KEEP_WORK:-0}" != 1 ]; then rm -rf "/scratch/project/k10018/gspa-v3/work/$SLURM_JOB_ID"; fi
 exit "$code"
}
trap finish EXIT
date --iso-8601=seconds > "$jobout/started-at.txt"; hostname > "$jobout/hostname.txt"; cp "$0" "$jobout/executed.sbatch"
singularity --version > "$jobout/singularity-version.txt"
# Writable bundle root: the immutable inputs copied once from the mirror (never training/, jobs/, slurm/).
rsync -a --exclude training/ --exclude jobs/ --exclude slurm/ "$data_root/rolefence-arm-v1/" "$out_root/"
( cd "$out_root" && sha256sum source.tar.gz claims/*.jsonl claims/claims_report.json > "$jobout/inputs.sha256" )
( cd "$data_root/rolefence-arm-v1" && sha256sum source.tar.gz claims/*.jsonl claims/claims_report.json | sed "s#  #  $out_root/#" | (cd / && sha256sum -c --quiet) ) > "$jobout/inputs-match-mirror.log" 2>&1
zstd -dc "$data_root/observed-gap-context-modes-v1-relocation/package.tar.zst" | tar -xf - -C "$work/fs"
tar -xf "$data_root/context-resume-v1/supplement.tar" -C "$work/fs"
tar -xf "$data_root/source-context-full-v1/parent-supplement.tar" -C "$work/fs"
mkdir -p "$work/fs$logical_source" "$work/fs/home/leechuck/tmp"
tar -xzf "$out_root/source.tar.gz" -C "$work/fs$logical_source"
( cd "$work/fs$logical_source" && sha256sum -c "$out_root/source.sha256" > "$jobout/source-check.log" )
container=(singularity exec --no-home --no-mount tmp --pwd "$logical_source"
 --bind /project/k10018/sw/envs/gspa-tools:/project/k10018/sw/envs/gspa-tools:ro
 --bind "$data_root:$logical_data:ro" --bind "$out_root:$root"
 --bind "$work/fs/home/leechuck:/home/leechuck"
 --bind "$work/fs/usr/lib/jvm/java-21-openjdk-amd64:/usr/lib/jvm/java-21-openjdk-amd64"
 --env TMPDIR=/home/leechuck/tmp --env JAVA_TOOL_OPTIONS=-Djava.io.tmpdir=/home/leechuck/tmp
 --env PATH=/usr/lib/jvm/java-21-openjdk-amd64/bin:/usr/bin:/bin:/project/k10018/sw/envs/gspa-tools/bin --env LC_ALL=C.UTF-8 --env HOME=/home/leechuck --env XDG_CONFIG_HOME=/home/leechuck/.config --env XDG_CACHE_HOME=/home/leechuck/.cache --env XDG_DATA_HOME=/home/leechuck/.local/share
 --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1
 --env PYTHONPATH="$logical_source/vendor:$logical_source/benchmark/v3" --env IJSON_BACKEND=python
 "$image")
printf '%s\n' "${container[@]}" > "$jobout/container-command.txt"
python=/home/leechuck/tmp/gspa-v3-testenv/bin/python
"${container[@]}" java -version > "$jobout/java-version.txt" 2>&1
"${container[@]}" "$python" -c 'import sys,numpy,ijson; print(sys.version); print(numpy.__version__,ijson.__version__,ijson.backend)' > "$jobout/python-version.txt" 2>&1
"${container[@]}" "$python" - <<'PY' > "$jobout/tests.log" 2>&1
import unittest
suite = unittest.TestSuite()
for pattern, expected in (('test_run_rolefence_arm.py',6), ('test_source_arm_context_comparison.py',16),
                         ('test_verified_batch.py',21), ('test_cohort_sequence_inputs.py',22)):
    tests = unittest.defaultTestLoader.discover('benchmark/v3/tests', pattern=pattern)
    assert tests.countTestCases() == expected, (pattern, tests.countTestCases(), expected)
    suite.addTests(tests)
result = unittest.TextTestRunner(verbosity=2).run(suite)
assert result.wasSuccessful() and not result.skipped
PY
[ "$stages" = none ] && { echo smoke-ok > "$jobout/smoke.txt"; exit 0; }
out=$root/training
for stage in $stages; do
  args=(--out "$out")
  if [ "$stage" = prepare ]; then
    args+=(--training /home/leechuck/.cache/gspa-v3-data/observed-gap-models-v1
      --context /home/leechuck/.cache/gspa-v3-data/cohort-context-modes-v2-inputs
      --jar "$logical_data/ec-prefix-identity-v3/jobs/10402/gspa-3.0.0.jar"
      --claims-dir "$root/claims" --claims-report "$root/claims/claims_report.json")
  fi
  date --iso-8601=seconds > "$jobout/$stage.started"
  /usr/bin/time -v -o "$jobout/$stage.resources.txt" "${container[@]}" "$python" "$logical_source/benchmark/v3/run_rolefence_arm.py" "$stage" "${args[@]}" > "$jobout/$stage.log" 2>&1
  date --iso-8601=seconds > "$jobout/$stage.finished"
done
echo done
