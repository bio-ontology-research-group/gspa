# Reporting parity: evidence and limitations

Three reporting defects are repaired: mixed scoring appeared unknown in the
portal; complete model abstention reports were marked untrusted; nested stage
assumptions and outcomes rendered as `[object Object]` in the results browser.

**Final evidence is `verification/replay-v4/`.** Seven stored production-output
sets pass through the packaged visualization CLI, real FastAPI endpoints and
Chromium. Sources, scores, candidate identities and all input proteins remain
intact. Two models each retain 26,328 mixed-scoring candidates (7,630 unsupported)
and 3,469 wholly unsupported candidates. Three synthetic absence-mechanism
outputs retain their uncertainty and budget/no-work states.

`verification/replay/` and `replay-final/` are intermediate harness runs: their
protein-only adapter omitted input proteins without annotations. `replay-v3/`
stopped on an incorrect assertion key for the browser's compact protein schema.
Do not use their protein totals as final evidence. The final adapter includes all
2,661 input proteins and reports 802 annotated proteins. Original production
prediction files were unchanged throughout.

Final gates pass 46 service tests and 73 visualization tests, including live
Chromium. Failing-before logs, the missing-FastAPI environment failure, package
versions, source snapshots, exact CLI commands, old/new portal summaries and
screenshots are retained. `verification/jar-verification.json` proves that the
new JAR changes only `visualize/make_viz.py`; all other ZIP entries are identical.
The immutable JAR is cached locally by its full SHA-256
`b410dadec00ead589f1f2fd0280ced8eea460755c202c41784f4de3df8a012d5`.

The source report gives reproduction commands and scope. Replay helpers retain
original workspace/cache paths in their receipts; restore or remap those inputs
when running elsewhere. `references/` holds the rendering ontology and EC map.
Large files use deterministic gzip. The storage index gives original sizes and
uncompressed hashes; run `sha256sum -c SHA256SUMS` to verify stored files.

This is reporting correctness, not biological validation. Batch provenance was
unavailable and remains empty. Protein-only inputs lack genomic coordinates.
No public service was deployed. No whole audit item is signed off.
