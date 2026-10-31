# Price production checkpoint v1

This is an **incomplete production collection**. Registration covers all 132
single-source parts and both prior settings; only the first completed cell is
retained here. The other cells, ensemble/cohort comparisons and biological scoring
remain open. Do not interpret registration or candidate preservation as accuracy.

* Final streaming replay 10226: existing 2.68GB ledger and 9,790,188-row prediction
  hashes match; 53.75 seconds, 23,424 KiB peak RSS, 24 focused tests pass.
* Broader local gate: 44 tests pass, with no JVM/model execution.
* Pilot 10228: both arms preserve all 3,481 candidates; 111.05 seconds total.
* Array 10229: indices 1–65, at most two concurrent tasks, started after pilot.
* Dense cell 10235: both arms for DPZ index 66, 64g heap/80G allocation;
  running at the updated checkpoint, with prior-version output comparison queued.
* Diagnostic 10230: largest CLEAN part, only stage0, 168g heap/192G allocation;
  running at this checkpoint. It cannot complete a two-arm cell.

`source/` contains the frozen executed runner/test/wrapper files. `replay-10225/`
also preserves the earlier parser revision. Registration hashes cover the remote
JAR, references, publications, layouts, Python/Java binaries and parser/native
library. Bulk matrix claims, reference files and full dense outputs remain remote;
they are not bundled or reverified by the portable archive verifier.

Run `python3 benchmark/v3/results/price-production-checkpoint-v1/verify_archive.py`
for local inventory, source-registration, retained replay and complete pilot
receipt checks. This does not certify the rest of the production collection.

See [the benchmark record](../../../../docs/benchmarks/v3-price-production-collection.md)
for interpretation, constraints and remaining work.
