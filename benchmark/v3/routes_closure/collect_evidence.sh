#!/bin/bash
# Copy the small evidence of the routes/reporting closure into the repository and write SHA256SUMS.
# Usage: collect_evidence.sh <work root: ~/tmp/routes-closure> <repo root>
set -euo pipefail
WORK=${1:?work root}; REPO=${2:?repo root}
OUT=$REPO/benchmark/v3/results/routes-reporting-closure-v1
mkdir -p "$OUT"/{plan,service,jvm,r2,r2e,r3,r4,p10}
cp "$REPO/tasks/v3-routes-reporting-closure-2026-11-01.md" "$OUT/plan/"
sha256sum "$REPO/tasks/v3-routes-reporting-closure-2026-11-01.md" > "$OUT/plan/plan.sha256"
cp "$WORK/service-tests.log" "$OUT/service/"
cp "$WORK/jvm-gate.log" "$WORK/jvm-gate-rerun.log" "$WORK/shadowjar.log" "$WORK/spock-focused.log" "$OUT/jvm/" 2>/dev/null || true
cp "$WORK/jars/jar-diff-rc3-rc4.json" "$WORK/jars/gspa-3.0.0-rc4.sha256" "$OUT/jvm/"
python3 - "$REPO" "$OUT/jvm/test-results-summary.json" <<'PY'
import sys, glob, json, xml.etree.ElementTree as ET
repo, out = sys.argv[1], sys.argv[2]
t = f = e = s = 0; mods = {}
for x in glob.glob(repo + '/gspa-*/build/test-results/test/*.xml'):
    r = ET.parse(x).getroot(); m = x.split('/')[-5]
    mods.setdefault(m, [0, 0, 0, 0])
    for i, k in enumerate(('tests', 'failures', 'errors', 'skipped')): mods[m][i] += int(r.get(k) or 0)
json.dump({'per_module': mods, 'total': [sum(v[i] for v in mods.values()) for i in range(4)], 'order': ['tests', 'failures', 'errors', 'skipped']}, open(out, 'w'), indent=1)
PY
# R2: single-genome routes
for r in routeA routeA-rc3-divergence routeA-noprov routeB routeB-rawclaims; do
  mkdir -p "$OUT/r2/$r"
  for f in run.log config.json integration-config.json workflow.json; do [ -f "$WORK/r2/$r/$f" ] && cp "$WORK/r2/$r/$f" "$OUT/r2/$r/"; done
  [ -d "$WORK/r2/$r/annot" ] && for f in protein_integrated.tsv.stages.json protein_predictor_status.json protein_timing.tsv; do [ -f "$WORK/r2/$r/annot/$f" ] && cp "$WORK/r2/$r/annot/$f" "$OUT/r2/$r/"; done
  [ -f "$WORK/r2/$r/integrated.tsv.stages.json" ] && cp "$WORK/r2/$r/integrated.tsv.stages.json" "$OUT/r2/$r/"
done
cp "$WORK/r2/env.sh" "$WORK/r2/routeC.log" "$WORK/r2/routeC-rc4.log" "$OUT/r2/" 2>/dev/null || true
for tag in rc3 rc4; do
  jid=$(cat "$WORK/r2/routeC$([ $tag = rc4 ] && echo -rc4).job" 2>/dev/null || true)
  [ -n "$jid" ] && mkdir -p "$OUT/r2/web-$tag" && for f in status.json manifest.json integration.workflow.json config.json run.log gspa_out/integrated.tsv.stages.json gspa_out/analyses.json; do
    [ -f "$WORK/r2/jobs-$tag/$jid/$f" ] && cp "$WORK/r2/jobs-$tag/$jid/$f" "$OUT/r2/web-$tag/$(basename $f)"; done
done
cp "$WORK/r2/comparison.json" "$WORK/r2/comparison-rawclaims.json" "$OUT/r2/" 2>/dev/null || true
# R2e: pooled four-source claims, portal path versus plain CLI flags
for a in portal cli; do mkdir -p "$OUT/r2e/$a"; for f in run.log command.json integration.workflow.json integrated.tsv.stages.json actions.tsv; do [ -f "$WORK/r2e/$a/$f" ] && cp "$WORK/r2e/$a/$f" "$OUT/r2e/$a/"; done; done
cp "$WORK/r2e/comparison.json" "$OUT/r2e/" 2>/dev/null || true
# R3: 18-genome cohort, CLI batch versus Nextflow on Shaheen
for f in comparison.json comparison.summary direct.json resources.json overlay.json runtime.config manifest-direct.tsv manifest-nf.tsv inputs.sha256 jar.sha256 java21.version nextflow.version trace.tsv direct.time nextflow.time direct.exit nextflow.exit nextflow.log; do
  [ -f "$WORK/cohort18/$f" ] && cp "$WORK/cohort18/$f" "$OUT/r3/"; done
[ -f "$WORK/cohort18/direct/batch_stats.tsv" ] && cp "$WORK/cohort18/direct/batch_stats.tsv" "$OUT/r3/"
mkdir -p "$OUT/r3/slurm"; cp "$WORK/cohort18/slurm/"*.log "$OUT/r3/slurm/" 2>/dev/null || true
cp "$WORK/cohort18/direct/GCA_000007005.1.integrated.tsv.stages.json" "$WORK/cohort18/direct/GCA_000007005.1.integrated.tsv.cohort.json" "$OUT/r3/" 2>/dev/null || true
# R4: states through the routes
cp "$WORK/r4/results.json" "$WORK/r4/agreement.json" "$WORK/r4/agreement-combined.json" "$WORK/r4/run.log" "$OUT/r4/" 2>/dev/null || true
mkdir -p "$OUT/r4/nextflow-rerun"; cp "$WORK/r4-nf/results.json" "$WORK/r4-nf/agreement.json" "$WORK/r4-nf/run.log" "$OUT/r4/nextflow-rerun/" 2>/dev/null || true
for s in requested missing failed; do d="$WORK/r4-nf/$s/nextflow"; [ -d "$d" ] || continue; mkdir -p "$OUT/r4/nextflow-rerun/$s"; for f in nextflow.log direct.log resources.json runtime.config manifest.tsv trace.tsv direct.json; do [ -f "$d/$f" ] && cp "$d/$f" "$OUT/r4/nextflow-rerun/$s/"; done; cp "$d/direct/"*.stages.json "$OUT/r4/nextflow-rerun/$s/" 2>/dev/null || true; for l in "$d"/work/*/*/results/*/*_integrated.tsv.stages.json "$d"/published/cohort/*/*_integrated.tsv.stages.json; do [ -f "$l" ] && cp "$l" "$OUT/r4/nextflow-rerun/$s/nf-$(basename $l)"; done; done
for s in requested missing failed no-work; do for r in cli annotate web nextflow; do
  d="$WORK/r4/$s/$r"; [ -d "$d" ] || continue; mkdir -p "$OUT/r4/$s/$r"
  for f in run.log nextflow.log direct.log workflow.json integration-config.json resources.json runtime.config manifest.tsv trace.tsv integrated.tsv.stages.json; do [ -f "$d/$f" ] && cp "$d/$f" "$OUT/r4/$s/$r/"; done
  [ -d "$d/annot" ] && cp "$d/annot/"*_integrated.tsv.stages.json "$OUT/r4/$s/$r/" 2>/dev/null || true
  jid=$(cat "$d/job.id" 2>/dev/null || true); [ -n "$jid" ] && for f in status.json run.log gspa_out/integrated.tsv.stages.json; do [ -f "$d/jobs/$jid/$f" ] && cp "$d/jobs/$jid/$f" "$OUT/r4/$s/$r/$(basename $f)"; done
  cp "$d/direct/"*.stages.json "$OUT/r4/$s/$r/" 2>/dev/null || true
done; done
# P10: reporting replays
for b in batch1-failed-attempt batch1-rc3 batch1 batch2 batch3 batch3-failed-web; do
  [ -d "$WORK/p10/$b" ] || continue; mkdir -p "$OUT/p10/$b"
  cp "$WORK/p10/$b/verification.json" "$OUT/p10/$b/" 2>/dev/null || true
  for c in "$WORK/p10/$b"/*/; do [ -d "$c" ] || continue; n=$(basename "$c"); mkdir -p "$OUT/p10/$b/$n"
    for f in command.json render.log portal-summary.json portal-text.txt browser-ledger-text.txt results-browser.png portal.png status.json; do [ -f "$c/$f" ] && cp "$c/$f" "$OUT/p10/$b/$n/"; done; done
done
cp "$WORK/p10/"cases-*.json "$WORK/p10/"*.log "$OUT/p10/" 2>/dev/null || true
cd "$OUT" && find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
du -sh "$OUT"; wc -l SHA256SUMS
