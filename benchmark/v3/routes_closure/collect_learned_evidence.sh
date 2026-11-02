#!/bin/bash
# Copy the small evidence of the learned-output / web-cohort closure (plan 2026-11-01) into the repository and write
# SHA256SUMS; files over 5 MB are listed with their hashes in LARGE_FILES_REMOVED.tsv and stay under the work root.
# Usage: collect_learned_evidence.sh <work root: ~/tmp/routes-learned> <repo root>
set -euo pipefail
WORK=${1:?work root}; REPO=${2:?repo root}
OUT=$REPO/benchmark/v3/results/routes-learned-cohort-v1
mkdir -p "$OUT"/{plan,service,jvm,a1,a2,a3,a4,a5,a6-ws,a6-shaheen,a7}
cp "$WORK/service-tests.log" "$OUT/service/"
cp "$WORK/jvm/"* "$OUT/jvm/"
cp "$WORK/engine/jar.sha256" "$OUT/jvm/engine-jar.sha256"
# A1: predictor regeneration with the contract
cp "$WORK/a1/report.json" "$WORK/a1.log" "$OUT/a1/"
for g in "$WORK"/a1/*/; do [ -d "$g/search" ] || continue; n=$(basename "$g"); mkdir -p "$OUT/a1/$n"; cp "$g/search/diamond.log" "$OUT/a1/$n/"; done
# A2-A5 single-genome routes
cp "$WORK/single/results.json" "$OUT/a2/results-A2-A5.json"
cp "$WORK/single.log" "$WORK/single-a5.log" "$OUT/a2/" 2>/dev/null || true
for g in "$WORK"/single/a2/*/; do n=$(basename "$g"); mkdir -p "$OUT/a2/$n"; for f in workflow.json run.log integrated.tsv.stages.json integrated.tsv.learned.json; do [ -f "$g/$f" ] && cp "$g/$f" "$OUT/a2/$n/"; done; done
cp "$WORK/single/a3/workflow.json" "$WORK/single/a3/manifest.tsv" "$WORK/single/a3/run.log" "$WORK/single/a3/batch_stats.tsv" "$OUT/a3/" 2>/dev/null || true
cp "$WORK/single/a3/"*.stages.json "$WORK/single/a3/"*.learned.json "$OUT/a3/" 2>/dev/null || true
for g in "$WORK"/single/a4/*/; do n=$(basename "$g"); mkdir -p "$OUT/a4/$n"; cp "$g/integration-config.json" "$g/config.json" "$g/run.log" "$OUT/a4/$n/" 2>/dev/null || true; cp "$g"/annot/*_integrated.tsv.stages.json "$g"/annot/*_integrated.tsv.learned.json "$g"/annot/*_predictor_status.json "$OUT/a4/$n/" 2>/dev/null || true; done
for g in "$WORK"/single/a5/*/; do n=$(basename "$g"); mkdir -p "$OUT/a5/$n"; cp "$g/run.log" "$OUT/a5/$n/" 2>/dev/null || true
  jid=$(cat "$g/job.id" 2>/dev/null || true); [ -n "$jid" ] && for f in status.json manifest.json integration.workflow.json config.json run.log gspa_out/integrated.tsv.stages.json gspa_out/integrated.tsv.learned.json gspa_out/analyses.json; do [ -f "$g/jobs/$jid/$f" ] && cp "$g/jobs/$jid/$f" "$OUT/a5/$n/job-$(basename $f)"; done; done
[ -d "$WORK/single/a5-failed-attempt-no-fastapi" ] && { mkdir -p "$OUT/a5/failed-attempt-no-fastapi"; cp "$WORK"/single/a5-failed-attempt-no-fastapi/*/run.log "$OUT/a5/failed-attempt-no-fastapi/" 2>/dev/null || true; }
# A6: three cohort routes, ws run (informational + web via the real endpoint) and Shaheen run
collect_cohort() { src=$1; dst=$2; mkdir -p "$dst"
  for f in build.log jar.sha256 java.version python.version nextflow.version direct.json resources.json resources-web.json runtime.config manifest-direct.tsv manifest-nf.tsv inputs.sha256 references.sha256 comparison.json comparison.summary direct.log direct.time direct.exit nextflow.log nextflow.time nextflow.exit web.log web.time web.exit trace.tsv; do [ -f "$src/$f" ] && cp "$src/$f" "$dst/"; done
  [ -f "$src/direct/batch_stats.tsv" ] && cp "$src/direct/batch_stats.tsv" "$dst/"
  mkdir -p "$dst/inputs"; cp "$src/inputs/cohort_metadata.json" "$dst/inputs/" 2>/dev/null || true
  for g in "$src"/direct/*/; do n=$(basename "$g"); [ "$n" = features ] && continue; mkdir -p "$dst/direct/$n"; cp "$g"/integrated.tsv.stages.json "$g"/integrated.tsv.learned.json "$g"/integrated.tsv.cohort.json "$dst/direct/$n/" 2>/dev/null || true; done
  jid=$(cat "$src/web/job.id" 2>/dev/null || true)
  if [ -n "$jid" ]; then mkdir -p "$dst/web"; for f in status.json manifest.json route.json cohort.workflow.json batch.tsv run.log; do [ -f "$src/web/jobs/$jid/$f" ] && cp "$src/web/jobs/$jid/$f" "$dst/web/"; done
    for g in "$src"/web/jobs/$jid/gspa_out/cohort/*/; do n=$(basename "$g"); [ "$n" = features ] && continue; mkdir -p "$dst/web/$n"; cp "$g"/integrated.tsv.stages.json "$g"/integrated.tsv.learned.json "$dst/web/$n/" 2>/dev/null || true; done; fi
  for g in "$src"/published/cohort/*/; do n=$(basename "$g"); [ "$n" = features ] && continue; mkdir -p "$dst/nextflow/$n"; cp "$g"/*_integrated.tsv.stages.json "$g"/*_integrated.tsv.learned.json "$dst/nextflow/$n/" 2>/dev/null || true; done
  cp "$src/published/cohort/batch.tsv" "$src/published/cohort/workflow.json" "$dst/nextflow/" 2>/dev/null || true
}
collect_cohort "$WORK/cohort3-ws" "$OUT/a6-ws"
[ -d "$WORK/cohort3-ws-attempt1-metadata" ] && { mkdir -p "$OUT/a6-ws/failed-attempt-metadata"; cp "$WORK/cohort3-ws-attempt1-metadata/comparison.json" "$WORK/cohort3-ws-attempt1-metadata/direct.exit" "$WORK/cohort3-ws-attempt1-metadata/nextflow.exit" "$WORK/cohort3-ws-attempt1-metadata/web.exit" "$OUT/a6-ws/failed-attempt-metadata/" 2>/dev/null || true; grep -h "Batch not run" "$WORK/cohort3-ws-attempt1-metadata/direct.log" | head -2 > "$OUT/a6-ws/failed-attempt-metadata/direct-error.txt" 2>/dev/null || true; }
[ -d "$WORK/shaheen/cohort3" ] && collect_cohort "$WORK/shaheen/cohort3" "$OUT/a6-shaheen"
mkdir -p "$OUT/a6-shaheen/slurm"; cp "$WORK"/shaheen/slurm/*.log "$OUT/a6-shaheen/slurm/" 2>/dev/null || true
cp "$WORK/a6-cross-machine.json" "$OUT/a6-shaheen/" 2>/dev/null || true
# A7 reporting
for b in batch-a2-failed-attempt-expectation-format batch-a2 batch-web; do [ -d "$WORK/a7/$b" ] || continue; mkdir -p "$OUT/a7/$b"; cp "$WORK/a7/$b/verification.json" "$OUT/a7/$b/" 2>/dev/null || true
  for c in "$WORK/a7/$b"/*/; do [ -d "$c" ] || continue; n=$(basename "$c"); mkdir -p "$OUT/a7/$b/$n"; for f in command.json render.log portal-summary.json portal-text.txt browser-ledger-text.txt results-browser.png portal.png status.json; do [ -f "$c/$f" ] && cp "$c/$f" "$OUT/a7/$b/$n/"; done; done; done
cp "$WORK/a7/"cases-*.json "$WORK/a7/"*.log "$OUT/a7/" 2>/dev/null || true
cd "$OUT"
: > LARGE_FILES_REMOVED.tsv
find . -type f -size +5M | while read -r f; do echo -e "$(sha256sum "$f" | cut -d' ' -f1)\t$f\t$(stat -c %s "$f")\tkept under $WORK" >> LARGE_FILES_REMOVED.tsv; rm -f "$f"; done
find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
du -sh "$OUT"; wc -l SHA256SUMS LARGE_FILES_REMOVED.tsv
