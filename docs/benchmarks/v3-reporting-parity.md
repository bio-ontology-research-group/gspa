# Production report and portal parity

The reporting audit found three defects using existing production artifacts:

- The portal recognized only learned and legacy scoring. Actual mixed-scoring
  runs appeared as unknown, and their learned reports were marked untrusted.
- A valid model abstention also appeared untrusted. The model had rejected all
  input profiles and retained baseline scores, but its report was useful evidence
  explaining that decision. The results browser already supported this case.
- Nested stage metadata became `[object Object]` in the results browser.
  Assembly uncertainty, hypothesis comparisons and structured assumptions were
  present in the exported JSON but unreadable in the stage table.

The portal now recognizes mixed scoring and displays the bound model identity
and unsupported-candidate count. A skipped learned stage can expose an abstention
report only when the scoring remains legacy, zero candidates were scored, the
model hashes match, and every reported candidate is unsupported with no score or
margin. Failed stages, mismatched hashes and reports containing scores retain
their previous rejection behavior. Missing reports remain unavailable.

The results browser now offers expandable, escaped JSON for structured stage
notes. It preserves arrays, empty arrays, objects, nulls and their nested values.
Scalar notes and stage statuses remain visible. This makes recorded assumptions
and unknown outcomes inspectable without changing any prediction.

## Production replay

Seven stored production-output sets were rendered using the rebuilt
`gspa visualize` command, then read through the actual FastAPI job-summary and
artifact endpoints. Chromium exercised both the results browser and the portal
page. HTTP requests were routed to a local TestClient; no public server or
prediction worker was started.

| Input | Runs | Candidate rows per run | Unsupported per run | Expected scoring |
|---|---:|---:|---:|---|
| Homology, logistic and GBDT | 2 | 26,328 | 7,630 | Mixed |
| Observed gapsmith, logistic and GBDT | 2 | 3,469 | 3,469 | Legacy; learned model abstained |
| Absence mechanism: unknown metadata, insufficient budget, disabled baseline | 3 | 1 | Not applicable | Legacy |

All seven final replays pass. Stage names, ordering, statuses and reasons match
the source ledgers. Browser counts reconcile with the TSV exports and all net
changed probabilities match. Every learned candidate identity, probability,
score status and supporting-source list matches its JSON report. There are no
unjoined learned candidates. The portal serves the exact learned-report bytes
and exposes its download link for all four valid reports. No page errors occur.

The protein-only replay retains all 2,661 supplied proteins, including the
1,859 without candidate annotations; its portal reports 802 annotated proteins.
The first replay adapter incorrectly used only annotated proteins for the total.
That harness defect, an intermediate compact-key assertion error, and their
logs remain in the archive. `replay-v4` is the final evidence. Inference files
were not changed by these adapter corrections.

The absence cases use a synthetic mechanism model. Their `UNKNOWN` values refer
to assembly/function-search metadata, not a measured biological absence or a
solver failure. Structured comparison records, including these unknown values,
are checked against the rendered detail text. An additional browser regression
covers an inconclusive outcome and HTML-like source text without executing it.

## Tests and packaging

Three new service regressions cover mixed scoring, valid complete abstention,
and rejection of invalid abstention reports. A new live-browser regression
covers nested notes, unknowns, empty arrays and escaping. The two original portal
failures and the original browser failure are recorded before the fixes.

The final gates pass 46 service tests and 73 visualization tests, including
headless Chromium. The initial full service gate lacked FastAPI; the repository's
pinned service requirements were installed in the local test environment, with
package lists retained before and after. The completed gate then passed.

The rebuilt JAR has SHA-256
`b410dadec00ead589f1f2fd0280ced8eea460755c202c41784f4de3df8a012d5`.
An entry-by-entry comparison with the preceding `a625288a…` JAR confirms that
only `visualize/make_viz.py` changed. Every other entry, including inference
classes, is byte-identical. The seven replay commands exercise that packaged
visualization resource; Java inference was not rerun or altered.

```bash
TMPDIR="$HOME/tmp" python -m unittest discover -s service -p 'test_*.py'
TMPDIR="$HOME/tmp" GSPA_BROWSER_TESTS=1 python -m unittest \
  scripts.tests.test_visualize_browser \
  scripts.tests.test_visualize_ledger \
  scripts.tests.test_visualize_learned
JAVA_TOOL_OPTIONS="-Djava.io.tmpdir=$HOME/tmp" ./gradlew --offline :gspa-cli:shadowJar
```

The evidence archive is `benchmark/v3/results/reporting-parity-v1/`. It includes
input/output snapshots, source hashes, exact CLI commands, old/new portal
summaries, screenshots, initial failures, package versions, final tests and JAR
comparison. Its storage index records uncompressed hashes; `SHA256SUMS` verifies
stored files, including compressed payloads.

This checks reporting correctness on stored inference, not biological accuracy.
Missing batch provenance remains empty; no provenance-completeness claim is
made. Protein-only inputs have no genomic coordinates. Broader route and
production-state coverage, independent biological validation, and the other
original acceptance requirements remain open. No whole item is signed off.
