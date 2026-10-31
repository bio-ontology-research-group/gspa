# Complete observed-gap production and matched reference comparison

All 72 genome/model runs and 1,513,136 candidate/model reader comparisons pass.
The matched analysis retains all 608 homology and 416 family-excluded test rows,
including unsupported fallback and taxon enforcement. Calibration-selected and
fixed 0.5 policies, all correction/harm counts and paired genome sensitivity are
retained. Homology rules have higher F1 than either learned arm; both learned
Brier resampling ranges include zero. Family-excluded improvements do not prove
a gap-context benefit: that test table has no gap-linked examples. No model is
promoted and no whole audit item is accepted. Reference-unconfirmed rows are not
experimentally proven negatives; all four null models remain untrained.

See [the report](../../../../docs/benchmarks/v3-observed-gap-cohort.md). This archive
binds the earlier complete export/training archive index and preserves all bound
inputs, production outputs, analysis policies, scored rows, code and test logs.
`index.json` maps original bytes to gzip objects with original/stored hashes;
every decompressed payload is checked. `SHA256SUMS` covers all files except itself.
