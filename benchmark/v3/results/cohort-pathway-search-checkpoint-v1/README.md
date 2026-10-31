# Cohort pathway-search checkpoint — incomplete

This checkpoint contains **eight verified bacterial searches out of 18 fixed
cohort genomes**, plus all **ten failed archaeal preflight attempts**. It is not
a completed context evaluation or biological validation. See the
[protocol and remaining work](../../../../docs/benchmarks/v3-cohort-pathway-acquisition.md).

The completed searches preserve **147,818 reaction rows** and **15,384 pathway
rows**, with unchanged inputs and matching reference hashes across genomes.
The independent audit checks the source/tool/proteome bindings, output bytes,
recomputed table counts, taxonomy and shared reference identities. Requiring
complete acquisition correctly exits nonzero and lists the ten missing genomes.

The frozen cohort contains 49,463 proteins across eight bacterial and ten
archaeal genomes. No function labels were used to select searches. Its initial
manifest SHA-256 is
`b9491693f7b34b29abe99c7f077c44f2b70b3003db258816d487d3badf2e8fa9`.
Array 10112 has terminated; its per-task scheduler records are retained.

Archaeal failures exposed an incomplete installed reference and the binary's
requirement for a `Bacteria` directory marker. The new recorder rejects missing
reference categories and supplies only an empty marker for archaeal staging.
Acquisition job 10120 has downloaded the pinned database and is publishing its
46,477 extracted files. Its wall time was extended to three hours without a
restart. Retry array 10131 waits on `afterok:10120`; neither is certified complete
by this checkpoint. Recheck current scheduler and receipts before acting.

Eighteen added tests cover observation/array execution contracts, actual
preflight regressions, input/output/reference substitutions, failed-attempt
preservation and archive extraction. The full benchmark gate passes **591 tests**,
with twelve optional skips. No production JVM code changed in this work.

Contents:

- `frozen-runs/cohort-pathway-search-v1/`: exact staged proteomes, manifest,
  runtime source and all first-array observations. Copied tool executables are
  omitted; their exact original hashes are bound in the receipts.
- `frozen-runs/cohort-pathway-search-archaea-v2/`: staged retry source, manifest
  and job script; no completed retry outputs are claimed.
- `frozen-runs/archaea-reference-v1-job/`: checksum-bound acquisition source and
  job script; the large reference publication is still remote and in progress.
- `verification/`: failed-before/passed-after tests, full gate, audit reports,
  rejected submission, scheduler snapshots and time-limit update evidence.
- `source/`: current source/tests/docs. `inputs/` retains the original selection,
  assembly provenance and split records used to freeze cohort membership.

`storage-index.json` maps original payloads to stored paths and SHA-256 values.
Large files use deterministic gzip; receipts refer to uncompressed bytes.
`source-files.json` binds current source snapshots; `SHA256SUMS` covers every
stored file except itself. This is an immutable checkpoint; subsequent completed
acquisition belongs in a new archive.
