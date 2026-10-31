# dbCAN contract repair

See the [change and scope](../../../../docs/benchmarks/v3-dbcan-contract.md).
This is parser/adapter evidence, not a biological benchmark or a full dbCAN run.

The archive retains the upstream source and acquisition hashes, controlled
component inputs, formatter outputs, exact before/after adapter sources, regression
test, remote batch scripts and test results. `validation.json` records the final
counts and archived execution status. `SHA256SUMS` binds every archived file.

Job 10182 exits zero. The old adapter fails 12 of 15 focused cases; the repaired
adapter passes all 15. The full predictor gate reports 413 tests with four
optional skips and no failures or errors. Exact XML cases, including the original
failure messages, are preserved in `job-10182/`.

`reproduce_formatter.py` executes only the source block identified by the
source hash and line bounds in `formatter-inputs/formatter-extraction.json`.
The source was inspected before execution. To reproduce both outputs:

```sh
python reproduce_formatter.py --source upstream/run_dbcan_v4.py \
  --inputs formatter-inputs --out "$TMPDIR/dbcan-formatter-new"
```

The output directory must not already exist. No model execution or database
download occurs in this replay. The two generated files must equal the fixtures
byte-for-byte. Full JVM validation runs on unimatrix under Slurm; its frozen
source package and runtime hashes are retained with the job records.

Job 10181 failed on node007 with scheduler signal 53 before creating any script
output. Its scheduler record is preserved. The subsequent validation uses node002.
The failed attempt is not counted as either a regression failure or a passing test.
