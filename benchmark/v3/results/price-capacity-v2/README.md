# Price capacity follow-up: partial evidence

EnzBERT GCF_000236665.1 contains 349,400 nonzero claims across 3,494 matrix rows.
Here stage 0 is the production no-priors combiner, not the raw-predictor maximum
score used as an accuracy baseline. Both complete production outputs preserve the full candidate universe and are
byte-identical to each other. This measures capacity/retention, not accuracy or
an improvement from priors. The independent bounded-memory key comparison
rechecks every candidate and 33 recorded input/output/code hashes.

`enzbert/` contains execution receipts, commands, resources and logs. Large TSVs
and stage ledgers remain at the recorded remote paths; their hashes were checked,
but those bulk files are deliberately not copied into this archive. The original
JAR and input/code identities remain frozen in the receipts.

The candidate comparer now externally sorts all identity keys with a 64 MiB
sort buffer on caller-provided scratch and compares every key. It rejects
missing, extra, duplicate and altered keys, without storing the candidate set
or output rows in Python memory. This does not establish bounded memory for
production inference or the separate ledger validation in verified_batch.py.

`launch/` freezes the scripts and JAR identity used by Slurm job 10037 for the
9,790,188-candidate DeepProZyme GCF_000236665.1 baseline. The job is ongoing at
this snapshot; check Slurm/logs rather than treating its launch as success.
The heap limit is 64 GiB within an 80 GiB allocation. Only stage 0 is requested.

The 70/132 export snapshot is receipt-only progress: published bulk bytes were
not rehashed for that count. All completed parts were independently checked by
the publishing workers. Acquisition, dense inference, cohort comparisons and
biological validation remain separate unfinished gates.

Tests: eight focused regressions pass; benchmark discovery reports 523 tests
with 11 skips (eight opt-in live Nextflow tests plus three other optional checks).
The earlier eight live Nextflow results remain applicable to the unchanged route
and production JAR. No JVM implementation changed in this follow-up.
