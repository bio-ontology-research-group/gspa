#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=gspa-wp4-search
#SBATCH --nodes=1
#SBATCH --time=08:00:00
#SBATCH --no-requeue
#SBATCH --output=/scratch/project/k10018/gspa-v3/out/wp4-search-v1/slurm/%j.log
# WP4 design §6: fenced gapsmith searches of one (set, superkingdom) manifest written by
# benchmark/v3/prepare_wp4_search_manifests.py, on one exclusive workq node.
#
# observe_gapsmith_cohort.py runs ONE genome per invocation (--index, requires SLURM_JOB_ID) and stages the whole
# reference (~46 k files) into $SLURM_TMPDIR/gapsmith-observation, which it neither shares nor deletes. So this
# job runs PARALLEL invocations at a time (xargs), each with its own SLURM_TMPDIR, SLURM_CPUS_PER_TASK = manifest
# threads (4), and deletes that staging directory as soon as the genome finishes (success or failure): at most
# PARALLEL reference copies exist at once. gapsmith and DIAMOND run natively (no container); Python is the gspa-v3
# conda env; the code is the read-only worktree slice. No package copies per worker.
#
# Environment:
#   MANIFEST         staged manifest on /project, e.g. /project/k10018/gspa-v3/wp4-search-v1/peers/archaea/manifest.json
#                    (its directory holds manifest.json, staging-hashes.json and proteomes/)
#   MANIFEST_SHA256  sha256 of that manifest (receipt.json "manifests" of the preparer)
#   PARALLEL         concurrent genomes (default 8; 4 threads each)
#   TASKS            space-separated 0-based genome indices (default: every genome without an observation dir)
#   STAGE_ROOT       parent of the job-private staging dir <STAGE_ROOT>/$SLURM_JOB_ID-wp4-search (default
#                    /scratch/project/k10018/gspa-v3/work; node-local memory such as /dev/shm if the reference fits)
#   KEEP_WORK=1      keep the staging dirs of failed genomes for inspection
# The search root is $OUT = dirname(output_root) = /scratch/project/k10018/gspa-v3/out/wp4-search-v1/<set>/<sk>;
# it receives manifest.json, staging-hashes.json, proteomes/ (copied once from /project), observations/<genome>/
# and jobs/$SLURM_JOB_ID/ (per-genome logs and exit codes). A genome whose observations/<genome> exists is never
# rerun here (the runner refuses to overwrite): a retry uses a new manifest and search root.
set -euo pipefail
: "${MANIFEST:?}" "${MANIFEST_SHA256:?}"
source /project/k10018/sw/eq-env.sh
export PY=/project/k10018/sw/envs/gspa-v3/bin/python
SRC=/project/k10018/code/gspa-v3
BASE=/scratch/project/k10018/gspa-v3/out/wp4-search-v1
export PARALLEL=${PARALLEL:-8} KEEP_WORK=${KEEP_WORK:-0}
export STAGE_ROOT=${STAGE_ROOT:-/scratch/project/k10018/gspa-v3/work}/$SLURM_JOB_ID-wp4-search   # job-private
mkdir -p "$BASE/slurm"
[ "$(sha256sum "$MANIFEST" | cut -d' ' -f1)" = "$MANIFEST_SHA256" ] || { echo "manifest sha256 differs: $MANIFEST" >&2; exit 2; }
BUNDLE=$(dirname "$MANIFEST")
# One python pass over the manifest: output root, threads, runner path, genome ids.
eval "$("$PY" - "$MANIFEST" "$BASE" <<'PY'
import json, shlex, sys
from pathlib import Path
m = json.load(open(sys.argv[1]))
assert m["schema"] == "gspa.cohort-pathway-search/1", "unsupported manifest"
out = Path(m["output_root"]).parent
assert out.parent.parent == Path(sys.argv[2]) and m["output_root"] == str(out / "observations"), \
    "output_root must be " + sys.argv[2] + "/<set>/<superkingdom>/observations"
assert (out.parent.name, out.name) == (m["set"], m["superkingdom"].lower()), "output root does not name set/superkingdom"
runner = [p for p in m["sources"] if p.endswith("/observe_gapsmith_cohort.py")]
assert len(runner) == 1, "manifest sources must name exactly one observe_gapsmith_cohort.py"
print("OUT=" + shlex.quote(str(out)))
print("THREADS=" + str(int(m["threads"])))
print("RUNNER=" + shlex.quote(runner[0]))
print("GENOMES=(" + " ".join(shlex.quote(g["genome_id"]) for g in m["genomes"]) + ")")
PY
)"
export OUT THREADS RUNNER
JOB=$OUT/jobs/$SLURM_JOB_ID
mkdir -p "$JOB"
date --iso-8601=seconds > "$JOB/started-at.txt"; hostname > "$JOB/hostname.txt"; cp "$(readlink -f "$0")" "$JOB/executed.sbatch"
printf '%s\n' "MANIFEST=$MANIFEST" "MANIFEST_SHA256=$MANIFEST_SHA256" "PARALLEL=$PARALLEL" "THREADS=$THREADS" \
  "STAGE_ROOT=$STAGE_ROOT" "TASKS=${TASKS:-}" > "$JOB/environment.txt"
# Stage the search root once (manifest, staging hashes, proteomes); a re-run must find identical bytes.
if [ -e "$OUT/manifest.json" ]; then
  cmp -s "$BUNDLE/manifest.json" "$OUT/manifest.json" || { echo "$OUT holds a different manifest" >&2; exit 2; }
else
  rsync -a "$BUNDLE/manifest.json" "$BUNDLE/staging-hashes.json" "$BUNDLE/proteomes" "$OUT/"
fi
# Fail before any search: staged bytes, frozen sources (binary, DIAMOND, runners) and runner path resolution
# (observe_gapsmith.py records Path(__file__).resolve(), which the audit compares with the manifest sources).
"$PY" - "$OUT" > "$JOB/preflight.json" <<'PY'
import hashlib, json, os, sys
from pathlib import Path
root = Path(sys.argv[1])
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()
staged = json.load(open(root / "staging-hashes.json"))
bad = [p for p, h in staged.items() if sha(root / p) != h]
m = json.load(open(root / "manifest.json"))
bad += ["source:" + p for p, h in m["sources"].items() if sha(p) != h]
bad += ["unresolved:" + p for p in m["sources"] if p.endswith(".py") and os.path.realpath(p) != p]
bad += ["proteome:" + g["genome_id"] for g in m["genomes"] if sha(g["proteome"]) != g["sha256"]]
bad += ["reference:" + m["reference"] + "/seq/" + m["superkingdom"] + "/" + c
        for c in ("rev", "unrev", "rxn") if not (Path(m["reference"]) / "seq" / m["superkingdom"] / c).is_dir()]
print(json.dumps({"staged_files": len(staged), "sources": m["sources"], "genomes": len(m["genomes"]),
                  "failures": bad}, indent=2))
sys.exit(1 if bad else 0)
PY
"$("$PY" -c 'import json,sys; print(json.load(open(sys.argv[1]))["diamond"])' "$OUT/manifest.json")" version > "$JOB/diamond-version.txt"
[ ! -e "$STAGE_ROOT" ] || { echo "staging dir exists: $STAGE_ROOT" >&2; exit 2; }
mkdir -p "$OUT/observations" "$STAGE_ROOT"
if [ -n "${TASKS:-}" ]; then
  tasks=$TASKS
else
  tasks=""
  for i in "${!GENOMES[@]}"; do [ -e "$OUT/observations/${GENOMES[$i]}" ] || tasks="$tasks $i"; done
fi
printf '%s\n' "${GENOMES[@]}" > "$JOB/genomes.txt"
echo "job $SLURM_JOB_ID on $(hostname): manifest=$MANIFEST genomes=${#GENOMES[@]} tasks=[${tasks# }] parallel=$PARALLEL"
run_genome() {
  local i=$1 gid stage code=0
  gid=$(sed -n "$((i+1))p" "$JOB/genomes.txt")
  [ -n "$gid" ] || { echo "no genome for index $i" >&2; return 0; }
  stage=$STAGE_ROOT/g$i; mkdir -p "$stage/tmp"
  (
    export SLURM_TMPDIR=$stage SLURM_CPUS_PER_TASK=$THREADS TMPDIR=$stage/tmp OMP_NUM_THREADS=1
    cd "$OUT"
    "$PY" "$RUNNER" --manifest "$OUT/manifest.json" --manifest-sha256 "$MANIFEST_SHA256" --index "$i"
  ) > "$JOB/$gid.log" 2>&1 || code=$?
  printf '%s\n' "$code" > "$JOB/$gid.exit"
  if [ "$code" -eq 0 ] || [ "$KEEP_WORK" != 1 ]; then rm -rf "$stage"; fi
  echo "index $i genome=$gid exit=$code"
  return 0
}
export -f run_genome
export JOB MANIFEST_SHA256
printf '%s\n' $tasks | xargs -r -P "$PARALLEL" -I{} bash -c 'run_genome {}'
fails=0
for i in $tasks; do
  gid=${GENOMES[$i]}; c=$(cat "$JOB/$gid.exit" 2>/dev/null || echo missing)
  printf '%s\t%s\t%s\n' "$i" "$gid" "$c" >> "$JOB/exit-codes.tsv"
  [ "$c" = 0 ] || { echo "FAILED index $i $gid exit=$c"; fails=$((fails+1)); }
done
[ "$KEEP_WORK" = 1 ] || rm -rf "$STAGE_ROOT"
date --iso-8601=seconds > "$JOB/finished-at.txt"; printf '%s\n' "$fails" > "$JOB/failures.txt"
echo "done, failures=$fails"; [ "$fails" -eq 0 ]
