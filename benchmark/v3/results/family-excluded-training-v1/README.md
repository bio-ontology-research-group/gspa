# Family-excluded reference-support training

Two models use the verified family-excluded predictor identity and production
features for all 18 genomes / 272,366 candidates. Training/calibration/test splits
and selection rules were registered before the run. Both models were frozen before
test evaluation. Test precision remains only 25–29% on supported reference-support
rows; these are not calibrated biological function probabilities. The null layer
remains untrained and no whole audit requirement is accepted.

The complete production replay checks 544,732 candidate/model pairs against the
original Python model reader, retains unsupported fallbacks and respects taxon
constraints on projection. The failed first checker and its source are preserved.
See [the report](../../../../docs/benchmarks/v3-family-excluded-training.md).

`index.json` maps original paths to gzip objects with stored and expanded hashes.
Bound historical source snapshots are retained separately from current benchmark
defaults. Plain frozen models are in `models/default/`. `SHA256SUMS` covers every
archive file except itself.
