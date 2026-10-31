# Community result ledger and summary repair

The old production verifier accepted missing solver ledgers, excessive call
counts and incomplete organism summaries. Twenty-three regression cases fail
before the repair. All 76 final focused cases and 1,356 JVM cases pass (six
optional skips). Fourteen old/new CLI runs prove all four malformed reports now
fail without publishing a community result; valid functional predictions,
candidate features and raw community results remain byte-identical. The initial
CLI checker's incorrect disabled-state expectation and baseline run are retained.

See [the report](../../../../docs/benchmarks/v3-community-result-contract.md).
This is synthetic real-SBML/LP contract validation, not biological cross-feeding
or gene-attribution evidence. No whole audit item is accepted. This archive
preserves both engines, before/after source, current module sources/resources,
commands, input/output hashes, complete test XML and failure/pass logs.
`index.json` maps original payloads to gzip objects; original and stored hashes
are recorded and every decompressed object is checked. `SHA256SUMS` covers every
file except itself. One unused first-attempt fixture was overwritten by the retry;
its expected hash and limitation are retained in unavailable-historical-inputs.json.
The first attempt stopped at the baseline before using that fixture. All final
run inputs and outputs are present and verified. Earlier archives and frozen numerical runs are unchanged.
