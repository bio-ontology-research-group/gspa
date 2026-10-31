# Largest Price cell and queued collection finalizer

Job 10295 completes both no-prior and default-prior integration on the complete
CLEAN GCA_004402975.1 matrix: 45,290,880 candidates per arm, with exact candidate
universe verification and successful checked publication. This checkpoint retains
the small completion, production, stage and publication receipts. All retained
published files match their recorded hashes. Large TSVs, stage ledgers and raw
inputs remain remote; they are not independently rehashed by this small archive.

The collection had 71 completed cells out of 132 at this checkpoint. The remaining
large-cell wave has started. Finalizer 10410 waits for successful arrays 10297 and
10298, then invokes the original registered finalizer on all cells. It checks
exactly one completed attempt per cell, both arms, candidate counts, every output
and input hash, publication identities and inventories, and writes the final
collection receipt exclusively. It requests 2 CPUs/4 GiB on node004 for up to
three days; it performs verification, not model inference.

The first queued wrapper, 10409, was cancelled before execution after inspection
found its assumed Python-file count was wrong. The fresh v2 wrapper checks the
exact registered source-file set; its read-only preflight passes on all eight
registered Python files, the original plan and registration. Both wrapper versions
are retained. The actual finalization has not run.

The benchmark uses the earlier frozen JAR (SHA-256
`0606bec0bf009bc6b9c1301a28ae8d21998730c36138c24cb5720a44eb31a50a`). These results
are neither current-JAR full-system validation nor a biological accuracy result.
No whole acceptance requirement is closed. `archive.sha256` verifies only the
files stored here when run from this directory.
