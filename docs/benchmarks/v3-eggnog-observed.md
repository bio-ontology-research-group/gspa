# Observed eggNOG contracts and annotation route

The old packaged CLI ignored `--threads 1` for eggNOG and always passed the
version-2 `--go_evidence` option. A replay with four visible JVM processors
launched version 2 with four CPUs; version 3 rejected the removed option and
produced a typed `FAILED` status with zero annotations. These failures are
retained as evidence, despite the CLI's overall zero exit code.

The repaired adapter executes installed eggNOG-mapper **2.1.13** and
**3.0.0-beta6**, with real upstream miniature databases, through annotation,
annotation with integration, and direct integration. This verifies these
versions and routes; it does not establish biological accuracy or support for
every future eggNOG release.

## Changes

- Typed `predictors.eggnog` configuration supplies `dataDir`, `executablePath`,
  `searchMode`, and `threads`. The existing `predictors.enable: [eggnog-mapper]`
  enables the adapter. CLI `--threads` overrides its configured CPU budget.
- The version probe reads past database warnings, requires a successful process,
  and recognizes the reported version. Version 2 keeps the non-electronic GO
  donor filter; version 3 omits that removed flag. Unrecognized major versions,
  invalid CPU counts, and unsupported search modes fail before the search.
  HMMER search mode is accepted only for version 2; the live checks here use
  DIAMOND for both versions.
- Named columns handle the different version-2 and version-3 tables. Duplicate
  column names fail. Headerless input is restricted to exactly 21 columns for
  the legacy version-2 layout; a version-3 table cannot silently use that layout.
- Every annotation retains the reported tool version, seed ortholog, E-value,
  bit score, ortholog group, and the legacy score transformation. The raw score
  is explicitly marked uncalibrated. Version-3 `annotation_confidence` is
  retained as its original per-source string, without treating it as a probability.
  Missing optional provenance remains missing; shifted headers cannot turn an
  E-value into a seed identifier.

The existing `min(1, -log10(E-value)/30)` transform, clipped below at zero and
with zero E-value mapped to one, is preserved. This change does not retrain the
integration calibration or certify its probabilities.

Example configuration:

```yaml
predictors:
  enable: [eggnog-mapper]
  eggnog:
    executablePath: /path/to/emapper.py
    dataDir: /path/to/eggnog-data
    searchMode: diamond
    threads: 2
```

## Inputs and independent reference execution

Version 2 uses the upstream [2.1.13 source and test fixtures](https://github.com/eggnogdb/eggnog-mapper/tree/2.1.13/tests/fixtures),
with two protein queries and its small eggNOG 5 database. Version 3 uses source
commit [b3757a6d226047527a729546e58ff530d76a5d7d](https://github.com/eggnogdb/eggnog-mapper/tree/b3757a6d226047527a729546e58ff530d76a5d7d)
and the author's [miniature eggNOG 7 self-test bundle](https://data.cgmlab.org/eggnog-mapper/emapper-3.0/selftest/eggnog-mapper-selftest.tar.gz),
containing six nif queries and their annotation closure. Each tool was installed
in a separate Python 3.10 environment. Both used DIAMOND 2.1.11.

Reference annotation tables were captured by direct tool execution before the
GSPA repair. They are stored byte-for-byte with input FASTA, execution log,
source commit, version and hashes under
`gspa-predictors/src/test/resources/observed-output/eggnog/v2` and `v3`.
These are genuine tool outputs from selected miniature databases, not independent
experimental function labels. The version-3 dataset intentionally self-matches
its reference seeds. No model-training or evaluation split is inferred from it.

An upstream version-reporting quirk appends the current working directory's Git
hash. Thus the printed suffix `498560f` in these outputs identifies the GSPA
working directory, not the eggNOG source. The separately recorded upstream
commit and installed-code hashes are authoritative; the original printed strings
are preserved.

## Validation

Ten of eleven initial contract regressions fail against the old source. The
final suite contains 15 passing cases, adding shifted optional provenance and
header-loss/ambiguity checks. Existing adapter failure/zero/missing-output tests
continue to pass with an explicit version-2 process stub. CLI configuration tests
cover both YAML defaults and the thread override. The complete JVM gate passes
1,293 tests with six optional skips, and the packaged JAR builds successfully.

`benchmark/v3/check_eggnog_route.py` executes three real CLI routes per version.
It verifies each claim identity and raw score against the separately captured
reference table, successful predictor status, actual CPU/database arguments,
version-specific flags, and preservation of score/seed provenance. Its final
receipts bind every database-directory file, including version 3's GO ontology,
the installed eggNOG source/extensions, DIAMOND, Java, JAR, inputs and outputs.
Installed dependency versions are retained separately.

The six final CLI runs pass:

| Installed version | Input proteins | Raw claims | Route result |
|---|---:|---:|---|
| 2.1.13 | 2 | 67 | Original claim identities/scores preserved; annotation/direct integration byte-identical |
| 3.0.0-beta6 | 6 | 26 | Previous failed invocation now succeeds; annotation/direct integration byte-identical |

Raw claims, including the new provenance, are identical between the two
annotation routes for each version. The new engine SHA-256 is
`705b865d9591746227c1e8d46475240dcb78426f49f5aa54f1153a695ee944f5`.
The old engine is
`256121146db83406f5e02bb20920f793d722da977404f37aaad2f1b932ef5cf9`.

Source snapshots, tests, all failure and successful logs, upstream source bundles,
the downloaded self-test bundle, generated database-cache files, tool manifests,
JVM XML and final route receipts are retained in
`benchmark/v3/results/eggnog-observed-v1/`. Earlier route receipts that omitted
the GO file from their input inventory remain preliminary; the final receipts
include the full database directory and were rerun after that correction.

Broader predictor accuracy, learned-model version/distribution eligibility,
other search backends, nucleotide dispatch, and the remaining original audit
requirements remain open. No whole item is accepted by this adapter repair.
