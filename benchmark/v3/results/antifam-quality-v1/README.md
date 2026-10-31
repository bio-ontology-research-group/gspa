# AntiFam quality preservation and scoring eligibility

This archive closes the specific quality-evidence loss in the working-tree review.
It does not complete items 2–4 or establish biological AntiFam accuracy. Read the
[implementation and interpretation](../../../../docs/benchmarks/v3-antifam-quality.md).

`cli/verification.json` binds seven actual CLI executions, their inputs and outputs.
The genuine installed HMMER 3.4 searches an explicitly synthetic HMM and sequences.
Both old predictors report SUCCESS, but the old engine drops AntiFam from raw
claims and from final annotations when integration replaces functional evidence.
The repaired engine retains both records, exports quality with no calibrated
functional probability, preserves regional metadata, and leaves the functional
prediction TSV byte-identical. Known-coverage annotation and direct integration
quality reports agree. No-hit, failed and unknown scan states remain distinct.

The original observed fixture is copied under `observed-synthetic-fixture/`, with
its own tool/input/output provenance. It is not a real AntiFam family. The GA
failure is genuine: this synthetic HMM contains no gathering threshold.

`logs/` retains the three failing-before tests, intermediate fixture/checker
failures, passing final gates and engine identities. `tests/` contains final XML
for the affected regression classes. `initial-order-check-failure/` preserves the
first CLI replay: its checker incorrectly assumed concurrent predictors returned
records in a fixed order. The corrected checker compares all records as a
multiset; nothing was dropped to obtain parity.

Final gates: 1,110 JVM tests, no failures/errors, two optional real-data demo skips;
552 benchmark tests with 12 discovery skips, followed by all nine live Nextflow
cases passing against the final JAR (three optional benchmark skips remain);
41 trainer and nine input-support Python tests pass. Twenty-eight JVM and four
Python cases were added in this repair. The final core suite is unchanged from its
preceding passing invocation; the last build reruns the affected predictor/CLI
suites after adding the coverage-binding regression.

Current models cannot interpret quality observations. Flagged candidates are
unsupported by learned scoring; existing heuristic evidence remains. Training
preparation and the trainer reject flagged rows. No models were refitted: the
numeric feature contract and no-quality export fixtures remain unchanged.

Before JAR SHA-256:
`94db839ff7c8eaf1bcb03b9d28d710d2cdbe45342f5941dfc520fe3809b0469c`.
Final JAR SHA-256:
`cb74fc910024c3d199a38e94a1a5d24466cb9740fc11fd88ff72071b67579122`.
Frozen JARs are in `~/.cache/gspa-v3-data/engines/<sha>.jar`. The source snapshots
are under `code/`; `source-files.json` hashes their original repository files.
This is a dirty working-tree repair based on release HEAD `498560f`, not a commit
or release. Preserve the old engine and artifacts; older engines may drop new
quality evidence and should not be used to replay the new raw contract.

To repeat from the repository root, use `benchmark/v3/check_antifam_quality.py`
with the five paths shown in `cli/runs.json` and `cli/verification.json`; `--out`
must be a fresh directory. Set `TMPDIR=$HOME/tmp` and
`JAVA_TOOL_OPTIONS=-Djava.io.tmpdir=$HOME/tmp`. The replay rechecks all input hashes
before accepting its receipt. The complete archive inventory is `sha256.json`,
excluding only that inventory itself.

Remote jobs were inspected but not restarted or changed. At the final observation,
10080 was running at 1:02:40, 10108 pending on its dependency, and 10109 running at
45:35 with 462/851 collection files rehashed. These are unfinished checks, not
capacity-comparison or collection-audit success claims.
