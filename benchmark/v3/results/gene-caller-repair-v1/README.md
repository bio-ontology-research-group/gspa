# Gene-caller repair evidence

Final source v4 passes 57 focused and 1,552 full JVM cases (six optional skips), job 10520. Job 10521 passes 23 actual annotation/import CLI invocations using native Pyrodigal 3.7.1 and official Prodigal 2.6.3. Normal whole-assembly gene exports are byte-identical before/after; true zero, native failure, unavailable execution, malformed/incomplete output and numeric-contig mapping are distinguished. Both imports retain 65 unscored memberships on the 4,359-protein Pyrodigal gene set (not the earlier 70 memberships on 4,915 supplied native antiSMASH CDS).

This is scoped contract and route verification, not gene/function prediction accuracy, learned confidence, or whole-item acceptance. Both single/meta modes use the same assembly, not a community benchmark. Controlled reordered/missing output tests are explicitly injected corruption.

Final JAR: `/data/hohndor/gspa-v3/gene-caller-status-v1/jobs/10520/gspa-3.0.0.jar`, SHA256 `1ae842e5905eecfff4ecbcb07bfa21d320ff9de0797273e4aeafd2abaccbdb5d`. The old executed JAR is 10481, SHA256 `069263f7b876222815928a637d02106b9c3b31a756edd1cf692432acc2e52303`. Frozen Price/source-context artifacts were not replaced.

The archive retains native observations, before failures, intermediate successes and final outputs, source/package manifests, four scoped source snapshots, reviews/adjudication and independent checks. 10513's compilation/wrapper failure and 10514's availability failure remain visible. 10514's final CLI XML is an old-phase result because Gradle stopped before rerunning it; it is not an after-pass. 10519 was cancelled while pending before execution, then replaced by final-source job 10521. The cancellation description is retained. No failed run is relabelled successful.

`retained-bindings.json` maps 826 original relative names to retained files and decompressed SHA256s. Larger files are gzip-compressed. Full source packages/JARs, Java, host Python, the native Prodigal executable and two wheels remain remotely hash-bound; the portable archive does not contain those binaries or complete runtime. Native source records and caller output are retained. Seven runtime binary bindings are explicitly enumerated by the verifier. Relevant sources are checked against all four full manifests; the full 843-file final snapshot is remote.

Run from any directory:

```bash
python3 benchmark/v3/results/gene-caller-repair-v1/verify_archive.py
```

The verifier checks the complete archive inventory, all retained bindings, source snapshots, failed/successful test gates, all output hashes, every native sequence/coordinate, unchanged gene exports, status/atomicity controls, numeric contigs and unscored import parity. It does not execute native tools or JVM tests.

Full report: `docs/benchmarks/v3-gene-caller-repair.md`. Remote root and all launch scripts: `/data/hohndor/gspa-v3/gene-caller-status-v1`; local preparation: `/home/leechuck/tmp/gspa-gene-caller-status-v1`.
