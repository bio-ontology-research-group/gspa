# Complete 18-genome pathway replay

All 1,001 payloads independently rehash successfully: 9,710,009,047 original bytes,
498,237,328 stored payload bytes. See the
[protocol and limitations](../../../../docs/benchmarks/v3-cohort-pathway-replay.md).

All 72 production CLI calls pass on the complete 18-genome search acquisition.
The independent table audit verifies 175,272 reaction rows, 19,424 pathway rows,
56,613 raw candidates, 497 eligible reaction groups and 521 gap records (464
with GO). Automatic and explicit gap routes match across 484,202 homology
candidates, 7,448 RLC rows and 373 suggestions; functional TSV bytes are unchanged.

This is a frozen-output replay of actual searches, not a new search, biological
accuracy estimate, trained-model comparison or test of the new absence analysis.
The confirmed reference overlap and weak-hit search reproducibility issue from
the preceding archive remain open. No whole audit item is accepted.

`cli/verification.json` binds 338 inputs and 992 outputs. `cli/` retains every
annotation and integration output; `source/` freezes the replay implementation
and protocol; `verification/` retains the complete run log, summary and archive
script. The engine is the prior context-mode engine:
`7f9d1d477098087469f90564d79a923ee28e85ba0a010c33840f2db44131671b`.

`storage-index.json` records original paths, sizes and hashes. Large files use
deterministic gzip; verification decompresses and rehashes every payload.
`SHA256SUMS` covers every stored file except itself. The archive-name field in
the original console metadata was corrected without changing any payload;
`archive-generation.log` and `archive-verification.json` preserve this correction.
