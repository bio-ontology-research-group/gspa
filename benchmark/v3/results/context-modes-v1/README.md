# Context modes and complete pathway acquisition

This archive is verified: 1,281 payloads, 16,270,693,645 original bytes and
959,252,822 stored payload bytes. See the
[protocol and limitations](../../../../docs/benchmarks/v3-context-modes.md).
No whole audit item is accepted.

- `acquisition/`: ten successful archaeal retry searches, reference acquisition
  and scheduler receipts, hashes and the combined 18-genome audit. The first
  ten failed searches remain in the earlier acquisition archive and audit history.
- `inputs/initial`, `exports/initial`: preserved preparer/loader mismatch and
  failed first matched arm; the successful default compatibility control remains.
- `inputs/final`, `exports/final`: corrected, hash-bound inputs and four complete
  eighteen-genome exports. Default mode exactly reproduces the preceding engine;
  all modes retain candidate/source identities, non-cohort features and functional
  prediction bytes.
- `source/`: source and tests used for these runs, including the model coverage
  and training changes. Receipts bind execution-time sources; later checkout
  edits do not change the archived evidence.
- `verification/`: initial failures, final gates, graph-loader ID parity,
  scheduler evidence, model scoring parity and archive scripts.

The final gates pass 1,245 JVM tests with two optional skips, 610 benchmark
tests with three optional skips (ten live Nextflow tests), 41 trainer tests
and 13 model-input-support tests. All four production export arms retain
484,202 candidates. Six models preserve all 4,262 supervision records;
25,572 actual JVM/Python row-model checks agree to at most 6.11e−16.
On the common 460 supported test rows, all paired genome-level Brier intervals
cross zero. Null outcomes remain untrained. No mode is promoted on these results.

The numerical feature contract remains v3. Existing universal-only models
abstain on new modes. Matched-model comparisons use the existing supervision
and splits; the target remains reference support on a previously examined cohort
with confirmed sequence overlap. These results do not establish unseen-family
accuracy, calibrated biological inference, null learning or flux feasibility.

`storage-index.json` maps payloads to original paths, sizes and SHA-256 values.
Files above 512 KiB use deterministic gzip; validation decompresses and rehashes
each payload. `SHA256SUMS` covers every stored file except itself.
The production engine remains in the immutable local cache as
`7f9d1d477098087469f90564d79a923ee28e85ba0a010c33840f2db44131671b.jar`.
