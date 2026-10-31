# HMMER: observed execution and annotation route

The 2026-10-26 audit reproduced three defects using HMMER 3.4:

* Availability used `hmmsearch --version`, which exits 1 because that option is
  unsupported. A valid executable and database therefore produced `UNAVAILABLE`.
* The parser accepted non-finite statistics and missing query accessions, and
  could produce a negative confidence from a finite negative bit score.
* `annotate --threads 1` did not reach HMMER; the observed command requested 16 CPUs.

The adapter now probes `hmmsearch -h`, requiring a successful exit, the program
header and the HMMER version line. Invalid E-values, non-finite scores and missing
accessions produce typed output-contract failures. Finite negative bit scores
remain in metadata, with confidence bounded to zero. The domain configuration
has a HMMER thread count, overridden by the CLI and passed through the pipeline.
Other adapters' thread handling has not been certified by this change.

## Evidence and limits

`benchmark/v3/generate_hmmer_fixture.py` invokes real `hmmbuild` and `hmmsearch`
on a deterministic **synthetic** alignment and two synthetic proteins. One is
an alignment member; the other is a polyalanine negative control. The model ID
`GSPA_TEST_001.1` is explicitly a test identity, not a genuine Pfam accession.
The fixture records commands, tool versions and binary hashes, inputs, model,
stdout/stderr and artifact hashes. This establishes execution and parsing
behavior, not sensitivity, specificity or biological validity on real Pfam data.

Seventeen adapter cases cover availability/version, observed columns and hashes,
malformed statistics/accession, bounded confidence, live successful/zero-hit runs
and failure on a malformed model. Three require `GSPA_HMMSEARCH_BIN`; all three
were executed against `/usr/bin/hmmsearch`. CLI configuration tests additionally
cover YAML thread defaults and command-line overrides; the pipeline regression
checks the actual `--cpu` argument. Failing before-fix logs are retained.

`benchmark/v3/check_hmmer_route.py` replays the old and repaired production JARs:

| Route | Predictor status | Raw claims |
|---|---|---:|
| Old annotation CLI | UNAVAILABLE | 0 |
| Repaired annotation CLI | SUCCESS | 1 |
| Repaired annotation CLI with integration | SUCCESS | 1 |

The old CLI exits zero despite the unavailable requested predictor, which is why
exit code alone cannot certify coverage. The repaired routes retain the expected
query accession and independent domain E-value/score, omit the negative control,
honor one CPU and preserve the pre-integration claim. Annotation-driven and
direct integration produce byte-identical one-candidate TSVs with priors disabled.
The replay hashes its inputs before and after execution and hashes all outputs.
It does not establish parity for other predictors or context-stage combinations.

The repaired JAR SHA-256 is
`8a9dedbd9358b1e1249dba3c3a2e1006e23c71e601e4d68beeb36fd00d27e65a`.
The old JAR is
`58dbebf1cf6ad97cd216d676b0cf1221600a97c7350f6f45abf27b54c9669654`.
Full JVM validation passes 991 tests with two skips. The benchmark suite passes
525 tests with three skips, including all eight live Nextflow cases. Build and replay logs,
observed artifacts, source snapshots and before-fix failures are archived under
`benchmark/v3/results/hmmer-observed-v1/`.

The older `mock-output/` fixtures lack verifiable generating inputs, executables
and run records. Their README and test descriptions now identify them as format
fixtures. They cannot substantiate the earlier blanket claim of genuine output
fixtures for every adapter. Items 2 and 3 remain open for the other execution
routes, genuine tool outputs and coverage-loss cases.

## AntiFam follow-up

The wider check reproduced the same availability defect in AntiFam, another
`hmmsearch` adapter. Both adapters now use the same executable probe. AntiFam also
requires a configured database, rejects non-finite/negative E-values, non-finite
bit scores and missing accessions, bounds its E-value confidence, and honors the
CLI thread count. Twelve new cases failed before the repair and now pass. They
include real matching/no-hit execution and failure when gathering thresholds
are requested from the synthetic model that has none. Successful observed runs
use E-value gating; genuine AntiFam models and detection accuracy remain open.

The follow-up build passes 1,003 JVM tests with two skips; all eight live Nextflow
cases also pass against that exact build. Actual old/new CLI runs
change AntiFam from UNAVAILABLE/zero annotations to SUCCESS/one expected test
PSEUDOGENE annotation. However, **both emit zero raw integration claims**:
AntiFam is absent from the current source-to-evidence map. This work does not
establish pseudogene-aware functional rescoring. The four-run direct HMMER
annotation/integration comparison also passes again after the shared extraction.

The follow-up JAR is
`0812564ad19764a65234f418bf2c539eb16578fb66ff62034421fd0ad39e4caa`;
`benchmark/v3/results/antifam-observed-v1/` records failed/passing tests, both
replays and frozen source. The preceding archive remains bound to its original
build and is not silently relabelled with this newer JAR.

## Reproduction

Use a new output directory for each invocation. Set `TMPDIR` and the Java temporary
directory to an actual disk-backed directory before running either command.

```bash
python benchmark/v3/generate_hmmer_fixture.py --out "$OBSERVED_FIXTURE" \
  --hmmbuild /usr/bin/hmmbuild --hmmsearch /usr/bin/hmmsearch
GSPA_HMMSEARCH_BIN=/usr/bin/hmmsearch ./gradlew test :gspa-cli:shadowJar
python benchmark/v3/check_hmmer_route.py \
  --before-jar "$OLD_JAR" --after-jar "$NEW_JAR" \
  --fixture "$OBSERVED_FIXTURE" --out "$REPLAY_OUTPUT" \
  --hmmsearch /usr/bin/hmmsearch
```
