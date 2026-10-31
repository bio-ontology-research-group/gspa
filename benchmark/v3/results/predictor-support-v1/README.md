# Predictor identity and score-range support

Eight frozen development models retain exactly the previous numerical scoring
parameters, calibration and thresholds while requiring verified predictor
identities and observed training/calibration score ranges. Four matched exports
preserve 484,202 candidates each and 4,262 supervision records per refit. Forty
old/new/baseline CLI runs and 34,096 JVM/Python comparisons validate eligibility
and probability parity. These are coverage mechanisms, not biological accuracy
or calibration guarantees. No whole audit requirement is accepted.

See [the report](../../../../docs/benchmarks/v3-predictor-support.md).
`index.json` maps original files to content-addressed gzip objects with original
and stored hashes. `bound-input:` entries use the sources executed during the
migration; `current-source:` entries include subsequent default-path updates.
Restore bound sources with their original relative paths in an isolated checkout
to reproduce historical receipts. Initial verifier failures are retained.
Plain `models/` files are byte-identical convenience copies for current runners.
`SHA256SUMS` covers all stored files except itself. Every decompressed payload
has also been independently hashed.

The final full JVM suite has 1,307 tests, zero failures and six optional skips.
`final-full-jvm:` entries preserve its complete XML after all Python checks;
Python visualization tests had replaced the earlier core XML with one focused
writer test. The initial 434-case XML inventory remains documented in
`archive-verification-initial.json`; the primary verification report records the
final complete suite. The benchmark gate has 667 tests (13 optional skips), all
10 optional live Nextflow tests pass separately, and 286 Python script tests pass
(16 optional skips).
