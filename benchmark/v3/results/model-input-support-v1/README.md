# Model input coverage evidence

This archive verifies source identity/coverage enforcement and report behavior.
It does not establish biological accuracy, valid predictor transfer, a trained
null layer, or completion of an entire acceptance item.

- `replay/`: original 18-genome feature replay, training manifest and input receipt.
- `models/`: frozen logistic/tree artifacts and supported-only test evaluation.
- `cli/`: six matched production executions against the final rebuilt engine.
- `browser/`: real mixed-report Chromium check, screenshot and reproducible adapter.
- `tests/`: JVM, Python, browser and live Nextflow test logs; the deliberately
  failing no-projection regression is retained beside final passing checks.
- `validation.json`: cache locations, final engine identity and rechecked files.
- `checksums.sha256`: hashes of every other archived file.

The source-coverage replay used the engine preserved under its hash in
`validation.json`. The final CLI/browser checks use the later engine containing
the report-only mixed-state fix. Absolute original paths stay in receipts so
these records are auditable; see the protocol for reproducing into fresh paths.
Large feature exports, prediction files and browser HTML remain in the cache.

Protocol: `docs/benchmarks/v3-model-input-support.md`.
