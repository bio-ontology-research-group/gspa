# Genomic TSV and geNomad native-output contracts

Status: before/after regressions, the full JVM gate, sixteen production CLI
fixture runs, and independent retained-output checks pass. No whole audit item
is accepted; these are contract checks, not predictor accuracy measurements.

The genomic TSV parser previously discarded its first line without validating
the header, accepted extra or missing fields, accepted infinite/out-of-range
scores, coerced unknown region types, and silently dropped malformed metadata.
A failed rerun could leave the preceding full-region records in genome metadata.
The parser now requires the declared six-column header and rows, finite scores
in `[0,1]`, supported genomic region types, and unambiguous metadata keys. It
validates rows before applying the score threshold. Each scan clears its prior
full-region records before attempting execution. Existing coordinate/coverage
checks still report partial results rather than inventing CDS membership.

The Python geNomad runner had an additional upstream failure boundary: it could
publish a header-only success when neither native summary existed, publish
only one summary, filter malformed rows into zero hits, or leave a partial or
stale output after failure. It now requires both native summary files, validates
their required columns, scores, lengths and provirus coordinates, and publishes
the six-column table atomically only after both summaries pass. A nonzero native
exit fails the sidecar. Valid empty native summaries remain successful zero-hit
results. Provirus names retain the original contig ID, including embedded pipes.

The native file requirements and coordinate conventions were checked against
the [pinned geNomad summary writer](https://github.com/apcamargo/genomad/blob/1886c7fa85053defba8751601acf81e2e7336b25/genomad/modules/summary.py)
and [upstream output documentation](https://github.com/apcamargo/genomad#understanding-the-outputs).
The writer opens both summaries and writes their headers before iterating over
predictions. Its source and commit metadata are retained with the checkpoint.
No model inference was inferred from those sources or replaced by a fixture.

## Verification

Unimatrix job 10173 reproduces thirteen JVM failures and fifteen failing Python
subcases. The other JVM cases retain valid zero/one scores, complete zero-hit
tables and the existing NaN rejection. Job 10174 passes all six Python test
methods, all 33 focused JVM checks, and the forced full gate: **1,384 tests,
zero failures/errors, six optional skips**, across 147 suites.

The CLI replay executes the actual JVM, actual old/new Python sidecars and a
controlled native executable. It does not run geNomad's model. All commands,
fixtures, source versions and engine hashes are registered before execution.

| Native fixture | Original outcome | Repaired outcome |
|---|---|---|
| Three valid regions, including a zero score | SUCCESS; three regions | Same; GFF byte-identical |
| Both valid headers, no rows | SUCCESS_ZERO | Same; GFF byte-identical |
| Both summaries missing | SUCCESS_ZERO | FAILED |
| Virus summary only | SUCCESS; incomplete result | FAILED |
| NaN score | SUCCESS_ZERO | FAILED |
| Infinite score | FAILED at JVM numeric parsing | FAILED at native parsing |
| Invalid header | SUCCESS_ZERO | FAILED |
| Native process exits 7 | FAILED | FAILED; GFF byte-identical |

The first replay in job 10174 incorrectly expected the original engine to accept
the sidecar's lowercase `inf`. Java rejects that spelling, although the original
JVM parser accepts `Infinity` when supplied directly. This fixture expectation
was corrected, not counted as a production defect. Its failed run is retained.
Job 10175 then completed all sixteen runs with exit zero using the same engines.
The packaged resource zipapp contains the exact repaired sidecar bytes.

An independent checker re-evaluates all saved predictor outcomes, region counts,
coordinates, scores, empty functional claims and integration-ledger presence.
It rehashes 240 CLI files, verifies current sources/tests against their tested
overlays, and checks the complete test reports and engine bindings. The final
engine SHA-256 is
`2fa391755c6b11451d175e437939e85a5a19291f6a68caefff5c3d86a042cd4c`.
All execution, builds and tests ran on Unimatrix node005.

## Remaining scope

These fixtures do not establish real geNomad execution or biological accuracy.
The JVM contract is shared by geNomad, CheckV and PhiSpy, but their native sidecar
parsers have not all been repaired or certified. PhiSpy still conflates a missing
coordinate file with zero hits; CheckV's missing completeness is not an observed
zero probability. VirSorter2/VIBRANT native parsing, complete manifest coverage,
configuration validation and broader input-ID handling also require review.
Scores retain their tool-specific semantics: a nominal call or completeness
estimate is not a calibrated probability of viral origin. A zero classification
after native filtering and `minScore` is not proof of biological absence.

Remote records: `/data/hohndor/gspa-v3/genomic-contract-v1/`.
Compressed evidence: `benchmark/v3/results/genomic-contract-v1/`.
The preceding [GFF export proof](v3-genomic-region-output.md) remains unchanged.
