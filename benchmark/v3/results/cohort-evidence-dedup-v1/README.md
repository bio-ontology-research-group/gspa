# Heterogeneous-strain evidence deduplication

See [the protocol](../../../../docs/benchmarks/v3-cohort-evidence-identity.md).

The cohort engine groups complete, identifier-independent observed inputs only
within each existing strain cluster. Distinct observations, independent clusters,
all output candidates and raw descriptive counts remain present. Canonical
permutation pools remove incidental record-order effects. Feature schema remains
version 3; cohort engine provenance advances to `cohort-v3`.

`replay/` preserves the original failure and corrected output: support stays 0.5
with 0, 1, 3 or 5 copies, with zero duplicate-induced change in every model column.
`source/` contains the production/test changes and hash-checked Python sources
bound by the refit receipt. `tests/` preserves the four failing pre-fix cases,
the intermediate test-fixture error and passing focused/full JVM gates.

The fresh export preserves 484,202 identities and all non-cohort values; only
5,043 permutation-excess and 893 permutation-p values change on the 18-genome
cohort. All 4,262 supervision records/splits are preserved. Both frozen models
and the table evaluation are included. Null layers remain untrained; these models
predict reference-annotation support within their declared coverage.

Use the pinned engine in `engine.json` with these models. The prior engine and
models remain in the contract-v3 archive. All six production arms complete and
independently verify, covering 118 genome outputs. `all-genome-duplicates.json`
checks all 18 original genomes: 346,269 supported numeric probabilities per model
are exactly invariant; 137,933 unsupported candidates remain unsupported, and no
candidate is missing. Frozen policies and all-case reports retain the negative
calibration-transfer and neighborhood-benefit results.

Final gates: 1,082 JVM tests (six optional skips), 551 benchmark tests (three
optional skips, including nine executed Nextflow cases). `sha256.json` inventories
the archive. No whole recommendation or independent biological accuracy claim
is made by this slice.
