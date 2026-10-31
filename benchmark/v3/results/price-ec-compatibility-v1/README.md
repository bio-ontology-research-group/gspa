# Current-JAR compatibility on one complete Price cell

Job 10405 runs the tested EC-identity/provenance JAR from 10402 against the full
TopEC GCF_000236665.1 cell (3,481 candidates), with no priors and with default
priors. Both production arms complete successfully and preserve the full input
candidate universe. The job then exits 1 because its comparison incorrectly
requires byte-identical TSVs across the intended identifier-format change.
That failed comparison and exit status are retained, not rewritten as success.

The independent `verify_archive.py` checks every retained TSV field and full
stage ledger, plus the original artifact hashes and 34 matching input hashes per
arm. All 3,481 identifiers change from bare EC numbers to `EC:`-prefixed numbers
in each arm. No other TSV field changes, including probabilities, source counts
and ordering. Ledgers agree after explicitly normalizing only EC spelling and
elapsed time. The two different JAR hashes and relocated manifest hashes remain
explicit; receipt contents are unchanged.

`before/` holds job 10292's complete small outputs; `after/` holds 10405's outputs.
The latter receipts retain their original scratch paths. The portable verifier
maps those paths to the archived files and checks their original hashes. It does
not claim to rehash external raw/reference inputs locally.

Run `python3 benchmark/v3/results/price-ec-compatibility-v1/verify_archive.py`.
This is one-cell compatibility evidence, not a complete Price comparison or a
biological accuracy result. The ongoing 132-cell collection still uses the earlier
JAR. No whole acceptance item is closed.
