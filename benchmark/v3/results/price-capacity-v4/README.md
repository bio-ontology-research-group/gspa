# Complete dense default-prior capacity

Job 10058 completed with exit 0 on 2026-10-27. The unchanged historical engine
retains all 9,790,188 candidates from the complete DeepProZyme GCF_000236665.1
matrix (3,494 rows). Java wall time is 2:44:34; peak RSS is 52,331,280 KiB
(49.91 GiB). End-to-end execution and candidate verification took 10,229.94 s.

`dpz/` contains the verification result, production receipt, commands, manifest,
batch statistics and measured resources. The production receipt hash is
`bea01577ef55c1ede45da30e6b2f70e4eb720260367510c60b1717777c41d3f2`.
The terminal Slurm state was captured before expiry. Bulk predictions and ledgers
remain at the receipt's remote paths.

The indexed-engine counterpart, job 10080, also completed with exit 0. Its
`indexed/` records retain all 9,790,188 candidates: 1:09:50 Java wall time,
63,324,216 KiB peak RSS (60.39 GiB), 4,550.47 s including verification.
Independent comparison job 10108 completed with exit 0. `comparison/` preserves
its receipt, log and terminal state; `comparison-launch/` preserves its source.
Every candidate identity is compared in both arms, the complete prediction TSVs
are byte-identical, and 37 bound input/output files are rehashed twice.

The TSV SHA-256 is
`64701c3dd55cd6dc753139df7b93eb10e49f1efad361ea943f19f42f7aa6136f`.
These measurements are capacity and equivalence, not accuracy or a controlled
speedup comparison. They do not validate larger CLEAN matrices or cohort inference
at that scale.
