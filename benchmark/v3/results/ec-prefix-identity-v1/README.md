# EC identity and claim provenance diagnostics

These are real integration-defect observations, not an accuracy benchmark.

- `observed-prefix-split.json`: probe 10389 uses the real GCA_000007005.1 genome.
  All 1,244 DIAMOND proteins also occur in the native CLEAN input. The sources
  agree on 511 protein/EC identifiers after prefix normalization, but the engine
  exports separate candidates for their bare and `EC:` forms. The observed
  combined count is 30,023; the canonical count should be 29,512. Two original
  claims, their contracts and input hashes are retained as the regression case.
- `parent-supervision-identity.json`: all 4,262 homology and 2,444 family-excluded
  parent supervision rows have distinct canonical identities. No duplicate-label
  merging or relaxed supervision check is needed for the spelling migration.
- `observed-contract-loss.json` and `roundtrip.groovy`: a JVM round trip of that
  real claim pair through the JSONL writer converts both predictor contracts
  from maps into strings. Both exported identities become null; numerical
  feature digests and source names remain unchanged. Numerical parity alone
  therefore misses the provenance loss.

The EC ingestion correction passes all 29 focused cases and the full 1,472-case
JVM run in job 10401 (six established skips, no failures). This follows updating
two output assertions that expected bare identifiers. The additional structured
metadata correction passes 30 focused cases and all 1,473 full-suite cases in
10402 (six established skips, no failures). The 105 focused Python cases also
pass. Corrected real-genome probe 10403 is running; its result remains pending. Raw acquisition files are unchanged.
These checks do not demonstrate improved enzyme accuracy or complete any of the
twelve whole acceptance requirements.
