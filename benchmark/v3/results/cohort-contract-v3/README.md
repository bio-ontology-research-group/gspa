# Contract-v3 migration evidence

See [the protocol](../../../../docs/benchmarks/v3-cohort-contract.md).

The migration drops eight descriptive cohort columns, weights two coverage
fractions, rejects old artifacts, and refits both models with unchanged supervision.
All 484,202 exported identities and non-cohort feature values are preserved; all
4,262 labelled rows retain their evidence, splits and locus groups. Frozen test
results remain reference-support measurements, not biological accuracy evidence.

Both models complete 59 genome/arm outputs each. Independent stage/hash checks
pass for all 118 outputs. `all-genome-duplicates.json` verifies the 18-genome exact-
copy control at tolerance 1e-12; `mixed-cluster-failure/` preserves a broader failing
case. General duplication invariance remains unachieved.

- `models/`: frozen logistic and tree artifacts; both null layers are untrained.
- `inputs.json`, `export.json`, `training.json`, `evaluation.json`: migration
  identities, preservation checks and frozen test evaluation.
- `controls/`, `*-policy.json`, `*-report.json`: execution receipts, logs and
  separately calibrated whole-genome policies/results.
- `source/`, `tests/`: original and current source/fixture snapshots, failing
  pre-fix evidence, passing gates and intermediate fixture-migration failures.
- `check_all_duplicates.py`: receipt-verified comparisons over original genomes.
- `sha256.json`: hashes of the archive; large exports/outputs remain at recorded
  cache paths. Historical inputs are preserved, not overwritten or recertified.

Final gates: 1,055 JVM tests (two skips), 543 benchmark tests (three skips, nine
live Nextflow cases), 40 trainer tests. No whole acceptance item is signed off.
