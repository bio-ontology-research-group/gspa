# Measured nutrient-environment data and model preflight

The [protocol and limitations](../../../../docs/benchmarks/v3-nutrient-environments.md)
describe this raw-data audit. It is preparation for empirical validation, not an
accuracy result or completion of an audit item.

The pinned author repository and version-3 preprint supply 16,608 wells and
2,464 assay summaries; 2,240 belong to primary environments. The fixed input
catalogue retains outcomes and missingness separately. All author means match
raw reconstruction. All 616 supplemented-export SEMs repeat the mean; every
discrepancy is recorded and uncertainty is recomputed from the raw replicates.
Source files remain unchanged.

- `data/raw/`: thirteen pinned author files, verified against Git blob hashes.
- `data/paper-v3.*`: retrieved primary methods, with version and identity checks.
- `data/derived-v3/`: current, source-bound inputs, measurements, observations
  and complete source/output hash report.
- `data/derived-v1/`, `data/derived-v2/`: historical intermediate outputs. Their
  reports retain implementation hashes, but those earlier generator source
  snapshots are unavailable; they are not the current reproducible release.
- `data/models/`: both supplied models and preflight results. All twelve gene
  symbols map; iML1515 serA affects three reactions and leuB has an isozyme bypass.
- `source/`: current preparation, model registration, tests and documentation.
- `verification/`: exact commands, initial acquisition metadata failure and
  final passing logs. Thirteen focused tests and 626 benchmark tests pass;
  the latter has three optional skips and ten executed live Nextflow cases.

No model outcomes were used to choose cases. Model GPRs here support known-genotype
phenotype testing, not independent gene-function truth. Flux assumptions, matched
simulations, actual sequence-loss cases and gene-attribution evaluation remain open.

`storage-index.json` records original paths, sizes and SHA-256 values. Large files
use deterministic gzip, and every payload was decompressed and rehashed.
`archive-verification.json` records counts and sizes. `SHA256SUMS` covers every
stored file except itself.
