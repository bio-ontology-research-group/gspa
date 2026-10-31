# Verified observed-gap export and model training

All 72 genome/arm executions preserve candidates, provenance and sequence evidence
while measuring the effect of 521 supplied observed gaps. There are 202 homology
and 106 family-excluded score changes; these are not biologically verified corrections.
All four models are frozen before test evaluation; all 13,412 Python/JVM table
comparisons agree. Null remains untrained. The benchmark gate passes 716 tests
with 13 optional skips. Full model production verification and matched-score
analysis are pending and are not certified by this archive.

See [the report](../../../../docs/benchmarks/v3-observed-gap-cohort.md). The earlier
failed-verifier archive remains unchanged. This archive retains all bound inputs,
complete verified exports, fitted models, test evaluation, reader parity, source
snapshots and logs. `index.json` maps original bytes to gzip objects with original
and stored hashes; every decompressed payload is verified. `SHA256SUMS` covers
all files except itself. No biological or whole-item acceptance is claimed.
