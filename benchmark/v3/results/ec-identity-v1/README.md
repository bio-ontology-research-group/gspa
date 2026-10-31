# EC identifier audit

The pinned ENZYME release of 02-Sep-2026 contains 8,472 entries. All entries and
transfer paths are retained. The old compact-extraction projection has 1,309
replacements: 1,202 resolve unambiguously to an active entry, 101 pick a branch
of a split, and six lead to deleted entries.

No labels were read and no production predictions, historical exports or metrics
were altered. The legacy projection remains necessary to reproduce old extraction
statistics; it is not a safe unconditional equivalence map for new scoring.

The compressed enzyme.dat preserves the official source, including its CC BY 4.0
notice and attribution. The source receipt records retrieval, release bytes and
SHA-256. The audit contains every entry, resolution and legacy disposition.
The source folder contains the audit code and eight unit tests. The verifier separately
reconstructs all 8,472 entries and 1,309 dispositions, including the actual R/S
acetoin split and distinct NADH/F420 reductases, without importing the audit code.

Run: python3 benchmark/v3/results/ec-identity-v1/verify_archive.py

This verifies the retained source audit, not a biological prediction benchmark.
See [the record](../../../../docs/benchmarks/v3-ec-identity.md).
