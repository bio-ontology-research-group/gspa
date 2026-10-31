# Shared stage-ledger validation

Two regressions are reproduced before the repair. Final local and remote job
10244 gates contain 106 cases: 105 pass and one optional JVM smoke test skips.
Both new verification entry points match the historical real DPZ artifact hashes
(2.68GB ledger, 9,790,188 prediction rows); replay peak RSS is 24,400 KiB.
No inference, fitting or label evaluation occurs in this replay.

`failed-10241/` preserves the failed test-bundle attempt and its source archive.
`remote/source.tar.gz` contains the complete 133-file final source/test bundle;
`remote/job-10244/result.json` binds every source and parser file read during replay.
`prior-replay.json` is the independently retained preceding verification record.
The actual large outputs and installed native-library bytes remain remote.

Run `python3 benchmark/v3/results/generic-ledger-v1/verify_archive.py` to check this
archive's inventory, all bundled source bindings, measured replay consistency and
test results. This is not an independent reread of the remote bulk artifacts and
does not certify biological accuracy, the full benchmark collection or a release.

See [the benchmark record](../../../../docs/benchmarks/v3-generic-ledger.md).
