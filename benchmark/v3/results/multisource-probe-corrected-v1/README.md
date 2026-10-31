# Corrected mixed-source real-genome probe

Job 10403 completes with exit zero on genome `GCA_000007005.1`. DIAMOND,
native CLEAN and their combined claims are exported under sequence-only,
UNIVERSAL and reaction-local-context-disabled controls: nine exports in total.
The combined arm contains 29,512 unique candidates, including 511 shared
protein–EC assignments retaining both evidence sources and their structured
predictor contracts. Original sequence claim files are unchanged.

Independent job 10413 completes with exit zero. Its stdlib reader imports no
GSPA ingestion or benchmark audit modules. It verifies all 145 published file
hashes and recomputes source identities, source contracts and score features
across 178,605 feature rows. All sequence-only context values are zero. The
independent receipt binds the original probe receipt by SHA-256.

The earlier independent attempt, 10411, passed the 145 output hashes and two
feature exports before failing on an incorrect assumed path for the disabled
reaction-local-context export. Its failure and reader are retained under
`independent-audit/failed-10411/`; the corrected reader and complete successful
receipt are retained separately.

This archive contains exact probe source, launch scripts, logs and receipts.
Full raw claims and feature outputs remain in
`/data/hohndor/gspa-v3/multisource-comparison-probe-v4/jobs/10403/results` on
Unimatrix. Independent audit outputs remain in
`/data/hohndor/gspa-v3/source-probe-audit-v2/jobs/10413`. The full published
outputs were rehashed there by job 10413; they are not duplicated in this
repository archive. `archive.sha256` covers the files retained here.

Both historical learned-output numerical checks also pass. They do not fit
new models. This one-genome result verifies integration and export behavior;
it does not establish improved function prediction, recurrence across genomes,
or acceptance of any whole requirement. The full source/context comparison is
tracked in [the comparison record](../../../../docs/benchmarks/v3-source-context-comparison.md).
