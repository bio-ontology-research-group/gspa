# Partial cohort pathway replay evidence, 2026-10-27

Status: **partial_verified, 8/18 genomes**. No whole audit item accepted.
See the [protocol, counts and limitations](../../../../docs/benchmarks/v3-cohort-pathway-replay.md).

- `cli/`: all 32 successful production CLI runs and their complete outputs,
  frozen configuration, 156 bound input hashes, 442 output hashes, and final
  receipt. Eight raw-table annotation replays preserve 42,795 candidates.
  Twenty-four matched integration runs preserve 264,687 homology candidates
  and verify automatic/explicit parity for 389 gap records and 282 suggestions.
  The receipt lists all ten missing genomes and does not certify full acquisition.
- `reference-overlap/`: frozen plan, scanner, scheduler script and completed
  job 10133 output. Full-sequence comparison confirms 395 proteins identical to
  named references, including 93 in test genomes. This is a lower bound and
  invalidates an unseen-sequence interpretation of these searches.
- `verification/`: full 600-test benchmark log (twelve optional skips), replay
  log, failed initial overlap job 10132, successful retry, archive script and
  repeat-search diagnostic. The repeat adds two weak EC candidates, without
  changing pathway predictions or the fifty production gap records.
- `source/`: replay/audit sources and nine added regression tests.
- `upstream-inspection/`: source relevant to a possible query-order explanation
  for repeated-search variation. This checkout is not proven to match the
  installed binary; it does not establish causation.

The replay uses frozen genuine search outputs, not fresh searches. Acquisition
evidence is in `../cohort-pathway-search-checkpoint-v1/`; the earlier single-genome
observation and gap replay remain preserved in their original archives. The
current production JAR is retained in the local engine cache, SHA-256
`1f0a7da89894896b05626d386879be494fd1e8d1ecbacba3b048e600086befd6`.
No JVM code changed; its prior 1,221-test JVM gate and nine live Nextflow checks
are not counted as new runs here.

`storage-index.json` binds **467 payloads**, original paths, sizes and
uncompressed SHA-256 values. Files larger than 512 KiB use deterministic gzip;
receipt hashes refer to decompressed bytes. All replay input/output hashes were
rechecked before archiving; every archived payload was independently decompressed
and rehashed. Original payloads total 5,369,228,411 bytes and stored payloads
280,914,039 bytes. `SHA256SUMS` covers every stored file except itself, including
the nested job checksum file. There are 470 stored files in total.

Outstanding: ten archaeal searches and full replay, search reproducibility,
reference/family independence, matched context and model comparisons, calibrated
locus/role/null prediction and biological validation. RLC suggestions remain
heuristic hypotheses. Successful execution does not satisfy these requirements.
