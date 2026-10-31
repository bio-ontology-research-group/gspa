# RLC mapping evidence, 2026-10-27

Verified engineering repair; no biological accuracy claim and no whole audit
item accepted. See [behavior, results and limits](../../../../docs/benchmarks/v3-rlc-mapping.md).

Standalone EC-to-GO input now reaches RLC without a pathway database. All mapping
alternatives survive shared loading; correlated labels do not independently add
diversity. The previous engine fails three new workflow tests, and an intermediate
fourth regression reproduces duplicate-label score inflation. Both the initial
and intermediate failures are retained under `verification/`.

- Full JVM gate: 1,221 tests, zero failures/errors, two optional skips.
- Benchmark discovery: 573 tests, twelve optional skips; nine live Nextflow
  cases pass separately, leaving three other optional checks unexecuted.
- `cli/rlc-mapping-v1/`: nine matched real-input CLI runs. The old standalone
  route yields zero RLC features; the repair restores 820 rows and 41 suggestions,
  byte-identical to the old pathway-backed route and invariant to mapping order.
  Empty/missing mapping, no-anchor and no-gap controls produce no suggestions.
- `cli/rlc-mapping-missing-evidence-v1/`: ten synthetic before/after CLI runs
  preserve missing-evidence behavior and explicit 0.5 evidence.
- `cli/rlc-mapping-gap-bridge-v1/`: ten real-output replay runs preserve all
  observed pathway rows, 3,469 raw candidates and 50 audited gap hypotheses,
  including annotation/direct/explicit route parity and failure controls.

The real mapping has 4,822 ECs and no multiple GO associations per EC. Its replay
proves standalone recovery and compatibility; the multiple-association checks
are synthetic. RLC scores remain heuristic, and the gaps remain hypotheses.

Engines are retained in the local immutable engine cache, not duplicated here:

- Before mapping repair: `fb7f750e451799a93446d6931a29817ec29adcd74ed00cf24828238fcdf3ee3b`.
- After: `1f0a7da89894896b05626d386879be494fd1e8d1ecbacba3b048e600086befd6`.
- Historical missing-evidence control: `cb74fc910024c3d199a38e94a1a5d24466cb9740fc11fd88ff72071b67579122`.

Each CLI receipt contains exact commands, Java/engine/input hashes and output
hashes. External searches are frozen-output replays of the verified gapsmith
observation, not new searches. `verification/jvm-summary.json` binds the final
engine to the full test counts; XML results are under `test-results/`.

`storage-index.json` identifies 504 uncompressed payloads, original paths,
SHA-256 values and stored paths. Files over 512 KiB use deterministic gzip;
receipt hashes refer to their uncompressed bytes. `source-files.json` records
current source snapshots. `SHA256SUMS` covers all stored files except itself.
Inputs and outputs were rehashed after every replay and again during archiving.
