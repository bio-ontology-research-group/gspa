# Essentiality indexing and key parsing cost

The live stack in `evidence/` samples frozen dense-prior job 10058 inside repeated
candidate-key parsing. The repaired implementation indexes eligible GO candidates
once per iteration, restores original ordering before ranking, and constructs
key arrays directly. It preserves scoring, floors and rescue choices.

The recorded old/new replay compares all 240 prior states from 80 synthetic
scenarios, including 90 states with rescues. Boosts, pinned floors and logs match
exactly. For an EC-only 100,000-candidate, 32-category microbenchmark, three
post-warm-up samples have medians 4.493 s before and 0.072 s after. This measures
only `beginIteration`, not application throughput or biological validity.

Two scan-count cases fail before repair (64,000 visits for 2,000 candidates);
the semantic-ordering cases already pass. The final build passes 1,010 JVM tests
(two skips) and 525 benchmark tests (three optional skips), including eight live
Nextflow cases. `code/` freezes the changed sources and semantic/performance
replay. `launch/` freezes the complete-matrix capacity scripts and code, excluding
the large JAR whose expected hash is recorded. Job status alone is not a pass.

`capacity/` records completed EnzBERT job 10079 and independent comparison 10082:
349,400 candidates retained in both builds; the entire prediction TSV is
byte-identical. Java takes 224.56 s, peak RSS 3,531,416 KiB (3.37 GiB), with an
8 GiB heap/12 GiB allocation. Total measurement including verification is 243.40 s.
The independent job rehashes all recorded inputs/outputs twice and compares all
candidate identities by external sorting. Bulk TSV and stage-ledger files stay
at their bound remote paths. The older run had a different heap/allocation and
concurrent load; this is not a controlled application-speedup claim.
Optimized dense DeepProZyme job 10080 remains dependent on original job 10058.

Pre-index JAR: ee72ddd579f606ded6f86a6a7b3af2c91877b6881ac315b5e86714000ae5b665.
Indexed JAR: 86805e4019bbe239b43b07664bd7a5d432a60b2a145907a08ddc807c5a6086cf.
See `docs/benchmarks/v3-essentiality-index.md` and `tasks/v3-acceptance.md`.

A later post-convergence stack sample from original job 10058 is preserved in
`evidence/post-convergence-10058.*`. It shows the same old array conversion in
`IntegrationState.splitFunctionKey`, now called by annotation synchronization.
The indexed engine already repairs this conversion; the sample does not establish
completion or end-to-end speedup. The live original job was not restarted.
