#!/bin/bash
#SBATCH --job-name=gspa-sulfonate
#SBATCH --partition=debug
#SBATCH --cpus-per-task=2
#SBATCH --nodelist=node002
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/data/hohndor/gspa-v3/sulfonate-diagnostic-v1/slurm/%j.log
set -euo pipefail
package=/data/hohndor/gspa-v3/empty-input-v1
current=/data/hohndor/gspa-v3/sulfonate-diagnostic-v1
work=/scratch/$SLURM_JOB_ID/gspa-sulfonate
jobout=$current/jobs/$SLURM_JOB_ID
mkdir -p "$work/fs" "$jobout"
export TMPDIR=/scratch/$SLURM_JOB_ID/tmp
mkdir -p "$TMPDIR"
finish() {
 code=$?
 trap - EXIT
 printf '%s\n' "$code" > "$jobout/exit-code.txt"
 if [ -d "$work/fs/home/leechuck/tmp/gspa-sulfonate-biology-v1" ]; then
  rsync -a "$work/fs/home/leechuck/tmp/gspa-sulfonate-biology-v1/" "$jobout/results/" || true
 fi
 exit "$code"
}
trap finish EXIT
hostname > "$jobout/hostname.txt"
cp "$current/run-v2.sbatch" "$jobout/"
sha256sum "$package/before.tar.zst" "$package/python-runtime.tar.zst" "$current/source-v2.tar" > "$jobout/packages.sha256"
zstd -dc "$package/before.tar.zst" | tar -xf - -C "$work/fs"
zstd -dc "$package/python-runtime.tar.zst" | tar -xf - -C "$work/fs"
tar -xf "$current/source-v2.tar" -C "$work/fs"
container=(singularity exec --no-home --no-mount tmp --pwd /home/leechuck/Public/software/gspa
 --bind /usr,/lib,/lib64,/bin,/etc
 --bind "$work/fs/home/leechuck:/home/leechuck"
 --env TMPDIR=/home/leechuck/tmp --env OPENBLAS_NUM_THREADS=1 --env OMP_NUM_THREADS=1
 --env PATH=/home/leechuck/tmp/gspa-v3-testenv/bin:/usr/bin:/bin
 /data/hohndor/gspa-v3/observed-gap-context-modes-v1-relocation/container-root)
base=/home/leechuck/tmp/gspa-sulfonate-biology-v1
runner=benchmark/v3/run_sulfonate_benchmark.py
"${container[@]}" python -m unittest discover -s benchmark/v3/tests -p 'test_sulfonate*.py' -v > "$jobout/tests.log" 2>&1
"${container[@]}" python "$runner" prepare --observations-sha256 d4e786d1260ce21ffcbc07ae9c14c282dcf4547be5c7686a3893418cb43c1f35 --inputs "$base/derived-v2/inputs.jsonl" --model "$base/iML1515.xml" --out "$base/plan" > "$jobout/prepare.log" 2>&1
set +e
"${container[@]}" python "$runner" run --plan "$base/plan" --model "$base/iML1515.xml" --out "$base/run" > "$jobout/run.log" 2>&1
run_code=$?
# Retain evaluation even if execution or the oracle reports an incomplete result.
"${container[@]}" python "$runner" oracle --plan "$base/plan" --run "$base/run" --model "$base/iML1515.xml" --out "$base/oracle.json" > "$jobout/oracle.log" 2>&1
oracle_code=$?
set -e
"${container[@]}" python "$runner" evaluate --plan "$base/plan" --run "$base/run" --observations "$base/derived-v2/observations.jsonl" --out "$base/evaluation.json" > "$jobout/evaluate.log" 2>&1
"${container[@]}" python benchmark/v3/check_sulfonate_medium.py --plan "$base/plan" --model "$base/iML1515.xml" --out "$base/medium-controls.json" > "$jobout/controls.log" 2>&1
if [ "$run_code" -ne 0 ]; then exit "$run_code"; fi
exit "$oracle_code"
