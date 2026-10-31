# EC candidate coverage v1

Read [the report](../../../../docs/benchmarks/v3-ec-candidate-coverage.md).

```sh
python3 benchmark/v3/results/ec-candidate-coverage-v1/verify_archive.py
```

This archive retains the 335-assertion, two-regime candidate trace and the complete
schema-2 production-report replay. `trace-inputs/` includes the original reference
annotation inventory, search hits, protein identity map, labels, splits and
reference-locus EC claims extracted from the complete raw claim files. The large
TSV sources are gzip-compressed; the verifier checks their uncompressed identities.
`trace.py` is the exact diagnostic executed against the original local cache and
preceding context-audit archive, not a portable reimplementation of the predictor.
`trace-report.json` binds all 41 original inputs used in that diagnostic.

`job-10224/` retains the terminal job, source/package hashes, nine focused remote
tests, old reports and the independent audit receipt, and all new report outputs.
Its registration explicitly records this as a replay of previously examined
predictions, without model execution or new independent validation. Only the
consumed prediction/table/reference files are restaged and rehashed; the upstream
complete-run audit is retained as provenance, not claimed rerun.

`source-v1.tar.gz` is the exact remote source overlay. `source/` retains final
production code and tests; production files equal the executed versions. The
cutoff test was strengthened after the remote replay and all nine focused tests
pass locally in `final-tests.log`. These are pure metric tests, with no models,
LPs or JVM execution on the laptop. `before-tests.log` retains the expected old
report failure. Both Sonnet reviews and their adjudications are included.

The portable verifier checks every retained file, reconstructs all 670 trace
rows from references/hits/extracted claims, and compares every new EC panel and
unchanged original summary/cutoff against the earlier independently verified
report. Raw full claim files, integrated TSVs and the execution runtime are not
copied here. Their recorded hashes and retained successful check receipts do not
replace those originals for a full source replay. No model, default, predictor
search parameter or historical result was changed.
