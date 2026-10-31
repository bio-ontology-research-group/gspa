# Completed community oracle replay

Job 10158 completed with exit code 0. All 6,600 numerical comparisons and 2,640
GPR comparisons agree. All 6,204 preserved records from the interrupted run match
exactly. Production predictions and registered tolerances are unchanged.

The collection audit covers 13,200 of 39,600 comparisons, with 13,199 agreements
and the original shard-0 solver unknown. Four shards remain outstanding at this
checkpoint; full numerical and empirical acceptance remain open.
See [the report](../../../../docs/benchmarks/v3-nutrient-community-audit.md).

`index.json` maps original paths and bound inputs to gzip objects with original
and stored hashes. All decompressed objects are checked. `SHA256SUMS` covers every
archive file except itself. This preserves the terminal replay, interruption
history, audit, inputs and scheduler record; the prior archive remains unchanged.
