#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=gspa-routes-cohort18
#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --no-requeue
#SBATCH --exclude=nid00022,nid00497,nid01033
#SBATCH --output=/scratch/project/k10018/gspa-v3/out/routes-closure-v1/slurm/cohort18-%j.log
#
# Matched 18-genome cohort comparison: CLI batch (`gspa integrate --batch-manifest`) versus Nextflow
# (`gspa-nf/cohort.nf`), same release-candidate JAR, same inputs. REPO_ROOT must be passed explicitly.
set -euo pipefail
: "${REPO_ROOT:?pass REPO_ROOT explicitly (never derived from the spooled \$0)}"
JOB=${JOB:-/scratch/project/k10018/gspa-v3/out/routes-closure-v1/cohort18}
JAR=${JAR:-/project/k10018/gspa-v3/release-candidate-3/gspa-3.0.0.jar}
EXPECT_JAR_SHA256=${EXPECT_JAR_SHA256:-8391d73fc1f6826492918a6116d40de4ba8147607ab406a895c29c32e067e80d}
FSROOT=${FSROOT:-/scratch/project/k10018/gspa-v3/shared/arm-fs-d04c9e11ee64db47}
JAVA21=$FSROOT/usr/lib/jvm/java-21-openjdk-amd64/bin
THREADS=${THREADS:-16}
unset JAVA_TOOL_OPTIONS JDK_JAVA_OPTIONS _JAVA_OPTIONS PYTHONPATH PYTHONHOME
mkdir -p "$JOB" "$JOB/tmp" "$JOB/direct"
export TMPDIR="$JOB/tmp"
echo "$(sha256sum "$JAR")" | tee "$JOB/jar.sha256"
[ "$(cut -d' ' -f1 "$JOB/jar.sha256")" = "$EXPECT_JAR_SHA256" ] || { echo "JAR hash mismatch" >&2; exit 2; }
"$JAVA21/java" -version 2>&1 | tee "$JOB/java21.version"

python3 "$REPO_ROOT/benchmark/v3/routes_closure/build_cohort_inputs.py" \
  --cohort /project/k10018/gspa-v3/training-cohort \
  --reference /project/k10018/gspa-v3/price-production-v1/reference \
  --gapsmith /project/k10018/gspa-v3/gapsmith-reference-bacteria-v1/reference \
  --model "$JOB/model_logistic.json" --out "$JOB" --threads "$THREADS"
sha256sum "$JOB/model_logistic.json" "$JOB"/manifest-*.tsv "$JOB/resources.json" "$JOB/direct.json" "$JOB/runtime.config" "$JOB/overlay.json" > "$JOB/inputs.sha256"
python3 - "$JOB/manifest-direct.tsv" "$JOB/resources.json" >> "$JOB/inputs.sha256" <<'PY'
import csv, hashlib, json, sys
def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as fh:
        for c in iter(lambda: fh.read(1 << 22), b''): h.update(c)
    return h.hexdigest()
for r in csv.DictReader(open(sys.argv[1]), delimiter='\t'):
    for k in ('claims_path', 'genome_layout_path'): print(sha(r[k]), r[k])
for k, v in json.load(open(sys.argv[2])).items():
    for p in (v if isinstance(v, list) else [v]):
        import os
        if os.path.isfile(p): print(sha(p), p)
PY

echo "== direct CLI batch: $(date)"
( export JAVA_TOOL_OPTIONS="-Xmx64g -Djava.io.tmpdir=$TMPDIR"
  /usr/bin/time -v "$JAVA21/java" -jar "$JAR" integrate --workflow-config "$JOB/direct.json" \
    --batch-manifest "$JOB/manifest-direct.tsv" --batch-threads "$THREADS" --batch-stats-out "$JOB/direct/batch_stats.tsv" \
    > "$JOB/direct.log" 2> "$JOB/direct.time" ) && echo "direct exit 0" | tee "$JOB/direct.exit" || { echo "direct exit $?" | tee "$JOB/direct.exit"; }
echo "== direct done: $(date)"

echo "== nextflow: $(date)"
module load nextflow/25.10.4
export NXF_OFFLINE=true NXF_ANSI_LOG=false NXF_HOME="$JOB/.nextflow"
export PATH="$JAVA21:$PATH"   # the process script's `java -jar` must be Java 21; the launcher keeps its own JAVA_CMD
nextflow -version | tee "$JOB/nextflow.version"
mkdir -p "$JOB/nf"; cp "$REPO_ROOT/gspa-nf/cohort.nf" "$JOB/nf/"; mkdir -p "$JOB/nf/modules"; cp "$REPO_ROOT/gspa-nf/modules/cohort.nf" "$JOB/nf/modules/"
( cd "$JOB" && export JAVA_TOOL_OPTIONS="-Xmx64g"
  /usr/bin/time -v nextflow -C "$JOB/runtime.config" run "$JOB/nf/cohort.nf" -work-dir "$JOB/work" -with-trace "$JOB/trace.tsv" \
    --cohort_manifest "$JOB/manifest-nf.tsv" --cohort_resources "$JOB/resources.json" --gspa_jar "$JAR" --outdir "$JOB/published" \
    > "$JOB/nextflow.log" 2> "$JOB/nextflow.time" ) && echo "nextflow exit 0" | tee "$JOB/nextflow.exit" || { echo "nextflow exit $?" | tee "$JOB/nextflow.exit"; }
echo "== nextflow done: $(date)"
cp "$JOB/.nextflow.log" "$JOB/nextflow.engine.log" 2>/dev/null || true

echo "== compare"
python3 "$REPO_ROOT/benchmark/v3/routes_closure/compare_cohort_routes.py" --direct "$JOB/direct" \
  --nextflow "$JOB/published/cohort" --manifest "$JOB/manifest-direct.tsv" --out "$JOB/comparison.json" | tee "$JOB/comparison.summary"
echo "== done: $(date)"
