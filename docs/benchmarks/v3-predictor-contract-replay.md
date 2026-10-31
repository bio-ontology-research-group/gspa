# Verified predictor identities for the learned-model migration

The current learned-model input gate checks actual source names, context mode,
and binary function/role/availability indicators. It does **not** check predictor
versions, reference versions, parameter identities or continuous score ranges.
This checkpoint prepares verified data for that outstanding requirement; it does
not claim that the model gate is repaired or that the models are recalibrated.

The homology training claims contained only `origin: swissprot_reference_hit` in
their metadata. Their separate historical report recorded the DIAMOND version,
command, reference hashes and hit-file hash. Attaching a guessed version to the
old claims would not establish the identity of a current executable. Instead,
`benchmark/v3/replay_predictor_contract.py` performs a fresh predictor execution
and requires it to reproduce the historical claims before adding metadata.

## Execution and preservation proof

The script verifies the recorded reference FASTA, reference annotations, old
alignment output, summary tables and every original per-genome claim file. It
accepts only the historical DIAMOND command structure, with its named output
columns, sensitivity, E-value and bounded thread count. Unknown or changed
options fail instead of being normalized away.

The fresh run rebuilds the reference database, runs DIAMOND, and independently
reconstructs both existing transfer modes from the new alignments using the
unchanged homology-transfer implementation. Every reconstructed file must match
its historical SHA-256 and record count. The complete genome set must match.

The observed replay uses DIAMOND 2.1.11, four threads, an E-value limit of 1e-5,
up to ten targets, and sensitive search. Its 164,787 alignment rows are
byte-identical to the historical hit file. All 18 genomes match in both modes:

| Transfer mode | Reconstructed claims | Preservation result |
|---|---:|---|
| Homology | 484,202 | Every historical file byte-identical |
| Family-excluded stress test | 272,366 | Every historical file byte-identical |

Each new claim receives `metadata.predictor_contract` with schema
`gspa.predictor-contract/1`. The contract identifies:

- Tool name, reported version and executable SHA-256.
- A combined identity for the reference sequences and reference annotations.
- Search and transfer parameters, including the family-exclusion thresholds
  when applicable.
- The homology adapter's source SHA-256 and the explicit score transformation:
  maximum supporting identity times capped query coverage, rounded to six decimals.

The two modes have different parameter identities. Query/genome identity is
bound in the run receipt rather than the reusable predictor contract. The
separate contract manifest retains the unhashed parameter/reference details.

Removing only the new metadata field from every enriched row restores the
original claim file byte-for-byte. All inputs are hashed again after execution.
Historical files remain untouched; rollback consists of continuing to use them.
Existing output paths are refused. Failed runs retain a failure receipt and log
and cannot receive `status: verified`.

## Checks and evidence

Seven tests cover complete reconstruction, reversible enrichment, separate mode
identities, missing genomes, changed historical inputs, changed fresh predictions,
failed search, and refusal of unknown command options. Contract hashes change
when the tool, reference annotations, parameters, transfer mode or adapter change.
The complete benchmark suite passes 659 tests with 13 optional skips.

The verified replay is in
`~/.cache/gspa-v3-data/predictor-contract-v1/`, with `inputs.json`,
`execution.json`, `contracts.json`, `run.json`, the new search output, exact
reconstructions, and enriched claim files. Reproduce into a fresh directory:

```bash
export TMPDIR="$HOME/tmp"
"$HOME/tmp/gspa-v3-testenv/bin/python" benchmark/v3/replay_predictor_contract.py \
  --baseline "$HOME/.cache/gspa-v3-data/training-cohort/baseline" \
  --prepared "$HOME/.cache/gspa-v3-data/training-cohort/prepared" \
  --out "$HOME/tmp/predictor-contract-replay"
```

Source, input/output hashes, original and enriched claims, reference/search
inputs, executable, test logs and verification receipts are archived in
`benchmark/v3/results/predictor-contract-v1/`. No JVM implementation or learned
artifact was changed in this checkpoint; the current engine remains
`705b865d9591746227c1e8d46475240dcb78426f49f5aa54f1153a695ee944f5`.

## Remaining work

The new metadata still needs to survive candidate extraction, feature export,
training-table preparation and reporting, and enter train/calibration-frozen
model eligibility in both Java and Python. Missing or changed identities need
explicit abstention, and continuous score support needs separate treatment.
Model refitting, current-route parity and unchanged valid numerical parameters
must then be verified against these freshly reproduced inputs. No model may be
declared version-compatible merely because this data-preparation step passed.

The reference labels are predictor inputs, not independent target truth. The
replay reads no target-label files and does not choose thresholds, splits or
claims from observed correctness. It preserves the earlier reference-overlap,
closed-world supervision, null-training and biological-generalization limitations.
No original audit item is accepted at this checkpoint.
