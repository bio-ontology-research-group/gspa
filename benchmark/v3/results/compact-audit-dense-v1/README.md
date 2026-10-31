# Full DPZ compact-audit parity

Job 10279 completed successfully. Both full 9,790,188-candidate arms match original
job 10235 in every TSV byte and every ledger line except exactly two stage-duration
lines. All candidate identities were compared. `source/dense_parity.py` records
the complete remote check; production files remain at the registered remote paths.

Run `python3 benchmark/v3/results/compact-audit-dense-v1/verify_archive.py` to verify
this portable receipt/resource archive. It checks hash chains and the common TSV
digest, but does not reread the large remote TSV/ledger files.

Peak JVM RSS falls 15.4% without priors and 4.2% with priors in these single paired
observations. Default-prior wall time increases from 37:56.95 to 40:30.11. These are
not repeated performance experiments or a general speedup. All twelve whole
acceptance requirements remain open. This success released the dependent full
production schedule; the full comparison is still running.
