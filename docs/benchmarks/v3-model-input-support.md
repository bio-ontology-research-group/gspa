# Learned-model input coverage

The previous learned-stage implementation applied a model to every candidate
whose numerical feature contract matched. That included predictors and roles
absent from its training/calibration data. A contract match establishes column
meaning and order; it does not establish calibration on a new input distribution.

The 4,262-row development table illustrates the problem: every row uses the
`swissprot_diamond_baseline` source, represented numerically by `seq.src.other`.
An unrelated unnamed tool could populate exactly the same numerical columns.
The earlier poor CLEAN transfer result remains evidence of a generalization
failure; this change does not turn it into an accuracy success.

## Source identity and frozen coverage

The production extractor now preserves sorted, normalized **actual source names**
in `CandidateFeatures.sequenceSources`. Only usable baseline claims contribute;
invalid scores and non-sequence/context claims do not. Pooling and numerical
features remain unchanged. Feature exports carry `sequence_sources`, the table
builder preserves it through truth joining/masking, and the trainer includes it
in its train/calibration digest. Missing metadata remains unknown; an explicitly
empty source list denotes a context-only candidate.

Trainer version 3 writes `input_support`, schema `gspa.model-input-support/1`.
Each profile records the exact source-name set and function-type, GO-aspect,
role, source/group-reporting and context-availability indicators. Both binary
classes must occur in both training and calibration for that profile to support
inference. Test rows never determine a profile or its class coverage.

These are necessary coverage checks, **not a calibration or generalization
guarantee**. They do not verify predictor versions, parameter settings, reference
database versions, continuous score distributions or taxonomic transfer.
Aliases with different actual names are conservatively distinct. More predictor
data and independent transfer validation remain necessary under item 4.

## Production behavior

An unsupported candidate has a reason and null learned margin/probability. It
retains its preceding annotation/posterior and pinned floor. It is never removed
by the learned retract threshold or promoted from an unsupported role. Supported
candidates retain the existing scoring behavior. If any member of a function/role
competition is unsupported, the entire competition receives no locus/null
probabilities; filtering cannot create a smaller apparently complete competition.

The stage report and feature export preserve every candidate and its support
status. Partial coverage declares `scoring_mode=mixed`, including supported
scores retained only in the report with no final annotation changes. The outcome
and change ledger say which scores were actually projected. With no supported candidates,
the learned stage is `skipped`, with the reason and the unchanged legacy scores.
The browser shows unsupported inputs explicitly, without describing them as
low-confidence predictions below the display threshold. A skipped report is
trusted for this purpose only when the model hashes agree and it contains no
learned scores.

Historical production artifacts without `input_support` remain inspectable but
cannot produce production learned probabilities. Re-export features and refit
from the recorded training/calibration splits; do not add a hand-written coverage
claim to an old artifact. Mechanism fixtures remain usable only through the
existing test-only route.

Python probability prediction applies the same gate. Frozen test evaluation
reports total/supported rows and exclusions; score metrics apply only to supported
rows, and recall over all test positives is reported separately. With zero
supported rows it returns unavailable metrics, not fabricated performance.
Locus evaluation excludes entire competitions containing unsupported members.

## Real-data replay

`benchmark/v3/replay_model_input_support.py` runs the built CLI over the original
18-genome cohort, checks every original numerical feature and candidate identity,
rebuilds the separately reviewed ThrRS table, and refits both model kinds with
their original settings. It requires unchanged numerical scoring/calibration
parameters, freezes artifacts before evaluating test coverage, and checks input
hashes again afterwards. This reuses previously examined data; it is not a new
blind benchmark.

The first replay omitted cohort context and was rejected by numerical-feature
comparison. A subsequent execution reproduced all 484,202 numerical rows, but
its requested RLC stage had no gaps and was correctly skipped; the strict batch
verifier refused to certify RLC execution. The final command retains cohort
context and does not request that empty RLC workload. These rejected runs remain
in their cache directories and are not completion evidence.

Reproduce with a new output directory and an already-built JAR:

```bash
export TMPDIR="$HOME/tmp"
export JAVA_TOOL_OPTIONS="-Djava.io.tmpdir=$TMPDIR"
"$HOME/tmp/gspa-v3-testenv/bin/python" benchmark/v3/replay_model_input_support.py \
  --original-manifest "$HOME/.cache/gspa-v3-data/reviewed-thrs-locus-v1/table/training_table_manifest.json" \
  --original-models "$HOME/.cache/gspa-v3-data/locus-truth-v2-final/models" \
  --jar "$HOME/.cache/gspa-v3-data/engines/fd05605d3f935227da63efc5ff3ada63bfaddd467665278737684f57b01cb5d1.jar" \
  --out "$HOME/tmp/model-input-support-replay"
```

The source-identity/coverage gate advances item 4 and the reporting requirements
of items 2/10. It does not complete their independent biological or route gates.

## Observed replay and production behavior

The successful replay preserves all **484,202** exported numerical vectors and
all **4,262** training rows' previous supervision and splits. Both model kinds
have exactly unchanged scoring/calibration parameters. Eight input profiles were
observed; six meet the binary-class coverage requirement. On the previously
examined test split, **462/608** rows are supported; 146 are withheld because
their profiles lack both classes in train/calibration. Metrics below apply only
to the supported subset and measure reference-annotation support, not experimental
correctness. They are not directly comparable to the old full-test metrics.

| Frozen model | Supported-row AUROC | Brier | Threshold precision | Recall over all 75 test positives |
|---|---:|---:|---:|---:|
| Logistic | 0.7576 | 0.1245 | 15/32 | 15/75 |
| Boosted trees | 0.7154 | 0.1277 | 21/59 | 21/75 |

The reviewed ThrRS competition remains one positive/negative training group.
There are still zero all-negative training groups and no complete calibration/test
groups. Both null layers explicitly remain untrained.

`check_model_input_support.py` runs six actual CLI jobs: matched baseline,
logistic and tree runs for real homology claims and compact CLEAN claims. The
homology genome, GCA_000007005.1, supplies 26,328 candidates: 18,698 supported
and 7,630 unsupported. Both models agree with Python on every coverage decision
and score (maximum absolute errors below 4e-15). Every unsupported function's
baseline TSV row is unchanged. This single-genome check omits cohort context;
the full-cohort replay above checks that context separately.

The CLEAN input for GCA_000006965.1 supplies 94,338 candidates. Both models
withhold every score, record `learned_integration=skipped`, clear no pinned floors,
and produce prediction TSVs byte-identical to the matched baseline. These are the
previously examined **compact/truncated** Price inputs. No full-proteome coverage,
accuracy improvement or successful predictor transfer is claimed.

The production mixed report was also rendered in Chromium: all 26,328 candidates
join actual NCBI protein IDs, the 7,630 unsupported rows are counted, an unsupported
candidate and its protein trace show the reason, and there are no JavaScript
errors. NCBI CDS metadata was adapted to the browser's input schema; the absent
batch provenance sidecar was explicitly empty, so this does not validate complete
provenance display. Fourteen browser regressions additionally cover the skipped
unsupported state and the existing interactions.

The final review reproduced a reporting defect with a focused regression: partial
coverage with no projected scores was labelled legacy, causing the browser to
withhold valid scored-only records. The workflow now bases mixed status on actual
candidate coverage and keeps projection/change counts separate.

Reproduce the six-run comparison after building the JAR:

```bash
"$HOME/tmp/gspa-v3-testenv/bin/python" benchmark/v3/check_model_input_support.py \
  --replay "$HOME/.cache/gspa-v3-data/model-input-support-v2-final" \
  --clean-claims "$HOME/.cache/gspa-v3-data/price-full-inputs/claims/clean/GCA_000006965.1.claims.jsonl" \
  --out "$HOME/tmp/model-input-support-cli"
```

Compact evidence is archived in `benchmark/v3/results/model-input-support-v1/`.
Large features, prediction files and browser HTML stay in the cache paths recorded
in its manifests. The replay engine is preserved by SHA-256 under
`~/.cache/gspa-v3-data/engines/`; later rebuilds do not retroactively certify the
earlier engine. Current CLI checks record their own engine and input hashes.

Engineering checks: 968 JVM tests (zero failures/errors, two skips), 47 Python
model tests, 37 browser-payload tests and 14 Chromium tests pass. The benchmark
suite passed 504 tests with three optional skips; all eight live Nextflow checks
also pass against the rebuilt JAR. These checks establish mechanics
and compatibility within their stated inputs, not acceptance of biological claims.
