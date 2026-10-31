# Observed-gap cohort failure checkpoint

The first verifier expected unchanged functional scores. Both homology production
batches completed, but that assertion failed because full-workflow gap inputs also
activate the gap-filling prior. This archive preserves that failure and its source.
A complete separate diagnostic finds 202 changes among 484,202 predictions:
198 increases, four decreases and 197 upward crossings of 0.5. Sequence likelihoods
and candidate identities are unchanged. These are heuristic changes, not verified
corrections. The corrected four-arm experiment is registered and running; it is
not accepted as complete at this checkpoint.

Twelve focused tests and the 708-test benchmark gate pass (13 optional skips).
See [the report](../../../../docs/benchmarks/v3-observed-gap-cohort.md).
`index.json` maps inputs and outputs to gzip objects with original/stored hashes;
every decompressed payload is checked. Bound historical source snapshots remain
separate from corrected sources. `SHA256SUMS` covers every file except itself.
