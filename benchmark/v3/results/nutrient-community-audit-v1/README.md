# Community audit checkpoint

This archive preserves two terminal shards and their provenance. The first has
6,600 independent comparisons: 6,599 agreements and one production solver unknown.
The interrupted second verifier has 6,204 preserved prefix rows but no completed
verification receipt. The full expected collection has 39,600 comparisons and
remains incomplete. Job 10158 is a separately recorded independent-verifier replay;
no partial prefix is reused as a new solve.

Sixteen audit regressions and the 683-test benchmark gate pass (13 optional skips).
Two completed numerical diagnostic jobs retain ten solver calls, the exact
production LP matrices and primal solutions. Alternative configurations solve the
problem at zero growth, while the default reproducibly returns unknown. The frozen
unknown is preserved, and no biological outcome is used for selection or tuning.
See [the report](../../../../docs/benchmarks/v3-nutrient-community-audit.md).

`index.json` maps original and bound-source paths to content-addressed gzip
objects. Original and stored hashes are provided; each decompressed payload is
verified. Initial auditor/test failures and earlier source snapshots are retained.
`SHA256SUMS` covers every stored file except itself. The checkpoint verifies the
archive and audit integrity, not full numerical or biological acceptance.
