# Genomic-region export and empty-input CLI checks

Status: repaired regression/full JVM gates and all twelve production CLI replays
pass, with independent retained-output verification.
No whole audit item is accepted.

Nucleotide predictors retain their original regions in genome metadata and make
separate copies for overlapping CDS. The GFF writer previously exported only
features and residue regions attached to proteins. Consequently, original
genomic regions disappeared from GFF, including successful predictions on
contigs with no CDS. The earlier empty-input dispatch repair alone did not fix
this output loss.

The writer now exports each original genomic record once as a `region`, before
the optional FASTA section. Coordinates remain 1-based and inclusive. Source,
raw score (including zero), annotation type, region type, evidence, metadata and
provenance are retained. Attributes are percent-escaped; structured metadata and
provenance use JSON within their escaped attributes. Generated IDs avoid existing
feature and residue-subfeature IDs. These rows have unknown strand and no
`Parent` or `Derives_from`: location does not establish a protein function.
Per-CDS copies do not create extra genomic rows. Disabling annotation export
continues to suppress these annotations.

The overlap regression also exposed incorrect protein matching: Groovy's `==`
on `Feature` uses its coordinate-based `compareTo`, so colocated CDS with
different IDs could receive the first protein's annotations. The writer now
matches the actual feature object or protein ID. Distinct colocated CDS retain
their own GO terms and residue subfeatures.

Five focused regressions check empty-proteome export and evidence, overlapping
CDS copies and ID collisions, annotation opt-out compatibility, and FASTA
ordering, plus distinct annotations on colocated CDS. Initially, three fail at the missing-region assertions and
the opt-out check passes (Unimatrix job 10168). Job 10169 built the intermediate
dispatch-repaired engine successfully, but its collection command expected an
incorrect `*-all.jar` filename and stopped before repaired tests. That attempt
is retained; job 10170 uses the configured `gspa-3.0.0.jar` filename. Its four
existing writer tests and three new tests passed; the overlap case exposed
the coordinate-equality defect above. Job 10171 repeats the original-writer
comparison with all five new tests before validating both repairs.

Job 10171 reproduces four intended failures against the original writer, then
passes all nine focused tests and the forced full gate: **1,369 tests, zero
failures or errors, six optional skips**, across 146 suites. Its first CLI
replay stops on a fixture error: it expected the genome-level operon predictor
to be skipped. That predictor correctly returns `SUCCESS_ZERO`. The corrected
fixture enables an available protein-level adapter on the no-CDS cases to test
`SKIPPED`, while keeping the operon expectation at `SUCCESS_ZERO`. Replay job
10172 uses the same built engines; no production source changes or repeat full
test gate are needed for this fixture correction.

The production replay registers twelve runs before execution: original engine,
dispatch-repaired engine, and final writer, each on no-CDS hits, overlapping CDS,
header-only zero results and missing output. A controlled Python sidecar emits
the genomic adapter's TSV contract. This exercises the packaged CLI, input
loading, dispatch, status handling, integration and GFF output. It is explicitly
not a geNomad execution or biological accuracy benchmark.

The completed checks verify one original region per hit, no invented protein
assignments, retained zero score and escaped metadata, an integration ledger
on empty inputs, `SKIPPED` protein predictors, and distinct successful-zero and
failed genomic scans. They also verify unchanged raw claims and integrated
prediction tables between the dispatch-only and writer repairs, and byte-identical
GFF for zero and failed scans.

Job 10172 completed with exit zero. An independent standard-library check
recomputed all twelve runs' region counts, coordinates, scores, metadata,
predictor states and empty-input ledger presence, and rehashed 177 retained
CLI files. It also verifies all source/test overlay bytes against the current
files, both rebuilt engine hashes, and all 146 final test reports. The original
engine omits both predictor execution and integration on no-CDS inputs; the
dispatch-only engine runs them but loses the genomic GFF rows; the final engine
preserves the rows. The controlled fixture's zero/missing-output distinction
does not resolve the real sidecar's native-summary ambiguity below.

The final engine SHA-256 is
`cc664c5bba2384c2f47fe202599c1e7a917bb387c1d56d482283b28e60b442cd`.
Its predecessor with only the dispatch repair is
`c543abac1d66e10bc542942822b1e74131d66d6f752bb831b6b8def2c2c2c776`.
The original comparison engine is the immutable `b23a7f6d…` used by the ongoing
context experiment. That experiment remains unchanged.

The subsequent [genomic contract repair](v3-genomic-contract.md) addresses the
adapter's header/nonfinite-score checks and geNomad's false successful-zero
result when native summaries are missing. Real-tool output/coverage and
biological validation remain open. The
gene-caller interface still conflates failure with genuine zero genes. These
limits prevent acceptance of the broader predictor and reporting requirements.

Execution package and job records:
`/data/hohndor/gspa-v3/genomic-output-v1/`.

The compressed evidence checkpoint is
`benchmark/v3/results/genomic-region-output-v1/`. It preserves the original
failures, intermediate defect, fixture/publication failures, exact tested
overlays, complete test reports, CLI outputs and independent verification.
Large runtime and engine binaries remain hash-bound on Unimatrix.
