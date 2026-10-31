# Predictor contract reconstruction

A fresh DIAMOND execution reproduces 164,787 alignment rows and all 756,568
historical claims across 18 genomes and two transfer modes. Enriched copies add
tool/reference/parameter/adapter/score-transform identities. Removing that one
metadata field restores every original byte. Seven focused tests and the full
benchmark gate (659 tests, 13 optional skips) pass.

This prepares verified data for the learned-model version/distribution gate.
Metadata propagation, eligibility, model refitting and biological validation
remain open. No model or whole audit item is accepted here.

`index.json` maps original source paths to gzip-compressed content-addressed
objects and records original/stored hashes. Duplicate inputs/reconstructions
share objects. Source, reference/search inputs, the executable, original and new
claim files, run receipts, test logs and documentation are retained. Historical
inputs remain untouched. `SHA256SUMS` covers every stored file except itself.
