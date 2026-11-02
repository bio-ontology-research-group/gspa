#!/bin/bash
# Three cohort routes on matched inputs (plan 2026-11-01, check A6): CLI batch, Nextflow cohort.nf, web cohort job.
# Environment (all explicit): JOB (output dir), INPUTS (claims/, layout/, cohort_metadata.json, cohort_orthogroups.tsv),
# GENOMES ("GCF_x=taxid ..."), JAR, EXPECT_JAR_SHA256, JAVA_BIN (dir with java 21), NEXTFLOW (launcher), PYTHON,
# REPO_ROOT, GO_OWL, EC2GO, REACTION_GRAPH, REACTION_ALIASES, MODEL, GO_OBO, THREADS, ROUTES (default "direct nextflow web").
set -euo pipefail
: "${JOB:?}" "${INPUTS:?}" "${GENOMES:?}" "${JAR:?}" "${EXPECT_JAR_SHA256:?}" "${JAVA_BIN:?}" "${PYTHON:?}" "${REPO_ROOT:?}"
: "${GO_OWL:?}" "${EC2GO:?}" "${REACTION_GRAPH:?}" "${REACTION_ALIASES:?}" "${MODEL:?}" "${GO_OBO:?}"
THREADS=${THREADS:-8}; ROUTES=${ROUTES:-"direct nextflow web"}; NEXTFLOW=${NEXTFLOW:-nextflow}
mkdir -p "$JOB/tmp" "$JOB/direct" "$JOB/web"
export TMPDIR="$JOB/tmp"
export PATH="$JAVA_BIN:$PATH"
unset JDK_JAVA_OPTIONS _JAVA_OPTIONS
echo "$(sha256sum "$JAR")" | tee "$JOB/jar.sha256"
[ "$(cut -d' ' -f1 "$JOB/jar.sha256")" = "$EXPECT_JAR_SHA256" ] || { echo "JAR hash mismatch" >&2; exit 2; }
"$JAVA_BIN/java" -version 2>&1 | tee "$JOB/java.version"
"$PYTHON" --version | tee "$JOB/python.version"
IDS=$(for g in $GENOMES; do echo "${g%%=*}"; done | sort | tr '\n' ' ')

"$PYTHON" "$REPO_ROOT/benchmark/v3/routes_closure/build_learned_cohort_inputs.py" \
  --claims-dir "$INPUTS/claims" --layout-dir "$INPUTS/layout" --cohort-metadata "$INPUTS/cohort_metadata.json" \
  --cohort-orthogroups "$INPUTS/cohort_orthogroups.tsv" ${LAYOUT_METADATA:+--layout-metadata "$LAYOUT_METADATA"} --genomes $GENOMES --model "$MODEL" --go-owl "$GO_OWL" --ec2go "$EC2GO" \
  --reaction-graph "$REACTION_GRAPH" --reaction-aliases "$REACTION_ALIASES" --out "$JOB" --threads "$THREADS" | tee "$JOB/build.log"

# A launcher for the web worker: the same JAR, Java 21, job-local temp.
cat > "$JOB/gspa-cli" <<EOF
#!/bin/bash
export TMPDIR="$JOB/tmp"
export JAVA_TOOL_OPTIONS="-Xmx24g -Djava.io.tmpdir=$JOB/tmp"
exec "$JAVA_BIN/java" -jar "$JAR" "\$@"
EOF
chmod +x "$JOB/gspa-cli"

if [[ " $ROUTES " == *" direct "* ]]; then
  echo "== direct CLI batch: $(date)"
  ( export JAVA_TOOL_OPTIONS="-Xmx24g -Djava.io.tmpdir=$TMPDIR"
    /usr/bin/time -v "$JAVA_BIN/java" -jar "$JAR" integrate --workflow-config "$JOB/direct.json" \
      --batch-manifest "$JOB/manifest-direct.tsv" --batch-threads "$THREADS" --batch-stats-out "$JOB/direct/batch_stats.tsv" \
      > "$JOB/direct.log" 2> "$JOB/direct.time" ) && echo "direct exit 0" | tee "$JOB/direct.exit" || echo "direct exit $?" | tee "$JOB/direct.exit"
  echo "== direct done: $(date)"
fi

if [[ " $ROUTES " == *" nextflow "* ]]; then
  echo "== nextflow: $(date)"
  export NXF_OFFLINE=true NXF_ANSI_LOG=false NXF_HOME="$JOB/.nextflow"
  "$NEXTFLOW" -version | tee "$JOB/nextflow.version" || true
  mkdir -p "$JOB/nf/modules"; cp "$REPO_ROOT/gspa-nf/cohort.nf" "$JOB/nf/"; cp "$REPO_ROOT/gspa-nf/modules/cohort.nf" "$JOB/nf/modules/"
  ( cd "$JOB" && export JAVA_TOOL_OPTIONS="-Xmx24g"
    /usr/bin/time -v "$NEXTFLOW" -C "$JOB/runtime.config" run "$JOB/nf/cohort.nf" -work-dir "$JOB/work" -with-trace "$JOB/trace.tsv" \
      --cohort_manifest "$JOB/manifest-nf.tsv" --cohort_resources "$JOB/resources.json" --gspa_jar "$JAR" --outdir "$JOB/published" \
      > "$JOB/nextflow.log" 2> "$JOB/nextflow.time" ) && echo "nextflow exit 0" | tee "$JOB/nextflow.exit" || echo "nextflow exit $?" | tee "$JOB/nextflow.exit"
  cp "$JOB/.nextflow.log" "$JOB/nextflow.engine.log" 2>/dev/null || true
  echo "== nextflow done: $(date)"
fi

if [[ " $ROUTES " == *" web "* ]]; then
  echo "== web cohort job: $(date)"
  ( export PYTHONPATH="$REPO_ROOT" GSPA_VIZ_PYTHON="$PYTHON"
    /usr/bin/time -v "$PYTHON" "$REPO_ROOT/benchmark/v3/routes_closure/web_cohort_job.py" --bundle "$JOB/bundle.zip" --jobs "$JOB/web/jobs" \
      --resources "$JOB/resources-web.json" --gspa "$JOB/gspa-cli" --threads "$THREADS" --go-obo "$GO_OBO" --feature-export \
      --job-id-file "$JOB/web/job.id" > "$JOB/web.log" 2> "$JOB/web.time" ) && echo "web exit 0" | tee "$JOB/web.exit" || echo "web exit $?" | tee "$JOB/web.exit"
  echo "== web done: $(date)"
fi

echo "== compare: $(date)"
ARGS=(--genomes $IDS --out "$JOB/comparison.json")
[[ " $ROUTES " == *" direct "* ]] && ARGS+=(--route "direct=$JOB/direct={g}/integrated.tsv" --features "direct=$JOB/direct/features")
[[ " $ROUTES " == *" nextflow "* ]] && ARGS+=(--route "nextflow=$JOB/published/cohort={g}/{g}_integrated.tsv" --features "nextflow=$JOB/published/cohort/features")
if [[ " $ROUTES " == *" web "* ]] && [ -s "$JOB/web/job.id" ]; then
  JID=$(cat "$JOB/web/job.id"); ARGS+=(--route "web=$JOB/web/jobs/$JID/gspa_out/cohort={g}/integrated.tsv" --features "web=$JOB/web/jobs/$JID/gspa_out/cohort/features")
fi
"$PYTHON" "$REPO_ROOT/benchmark/v3/routes_closure/compare_cohort_outputs.py" "${ARGS[@]}" | tee "$JOB/comparison.summary" || true
echo "== done: $(date)"
