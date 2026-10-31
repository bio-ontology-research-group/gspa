# Price full-matrix work: partial evidence snapshot

This is an in-progress acquisition/capacity record, not completion of the full
benchmark and not biological validation. See `export-snapshot.json` for the capture
time, published-job counts and explicit limitations. Inspect current Slurm state
before relying on the job IDs in `slurm-state.txt`.

`plan/` freezes all 132 matrices (1,526,919,081 complete/preliminary nonzero claims).
`repaired-parts/` records two real failures repaired by reproducing the historical
compact statistics projection; raw claims remain unchanged. The original failure
receipts and diagnostic remain under `superseded/`. Published claims and other
large inputs stay at the absolute paths recorded in receipts.

`projection-reference/` contains the historical mapping, vocabulary and extraction
code used solely to verify compact target vectors. The mapping pickle hash matches
the existing PSB source snapshot. This is not a current-nomenclature certification.

`capacity-topec/` records actual baseline/priors production inference with complete
candidate retention and resource measurements on one sparse matrix. It does not
prove dense/cohort capacity. The larger EnzBERT measurement is a separate job.
`tests/` retains the failing projection regression and final passing tests.

Protocol: `docs/benchmarks/v3-price-streaming-export.md`.
