# Production anomaly comparison, version 2

Status: implementation and execution verification; **biological accuracy is not
established**. This replaces the invalid case generator used in the
[historical comparison](v3-stage3-agent-evaluation.md). It does not complete
recommendation 11 or narrow its acceptance requirements.

## Information boundary

`benchmark/v3/run_agent_evidence.py` invokes the production CLI on all three
calibration and four test genomes of the existing development cohort. The
baseline and both learned arms consume identical claims and reference/context
inputs. Each arm has a content-verified run receipt and schema-2 stage ledgers.
This preparation reads split metadata, never evaluation labels.

`evaluate_stage3_agent.py prepare` reads the baseline's actual
`AnomalyExtractor` exports. It takes every exported case, generates opaque IDs,
and joins the two learned scores from real production TSVs. It rejects legacy
unverified runs and truth-bearing fields. Missing model scores remain missing;
no decoys, scores or annotations are constructed from truth. Every comparator
receives the same candidates, evidence and learned scores. The context comparator
uses the predeclared single feature `in_pathway_operon`.

The `suggest` and `agent` commands have no truth argument. All methods pass
through the same final candidate, function, ontology and taxonomic constraints.
Constraint rejection produces abstention, not an invented correct demotion.
The agent gets two turns, temperature zero, seed 42 and four concurrent workers;
its inspection tools only expose evidence already available to every comparator.
HTTP attempts include retries; the report records tokens returned by the endpoint
and explicitly identifies unavailable usage or billing cost.

`prepare_agent_truth.py` is a separate private evaluation step. It joins
experimental Swiss-Prot positives and their GO ancestors to the frozen public
cases. It preserves source accession, evidence code and identity-match basis.
It does not turn missing annotations into negative labels.

`calibrate` rejects truth from test genomes. The rule was fixed before inspecting
test outcomes: require at least ten active calibration decisions and at least
90% demonstrated precision, counting unknown decisions against that precision;
maximize correct decisions, then minimize harm, then choose the highest tied
threshold. If no threshold qualifies, abstain. Policies are immutable and bound
to the public inputs and calibration predictions. This is an empirical threshold,
not a guarantee of population precision.

`evaluate` accepts a frozen policy and separate test truth. It reports correct,
harmful, unknown, unchanged, abstained and failed decisions. Promoting an already
called function is not a recovery. Relaxation needs evidence for both removing
the original assertion and adding the broader one. Missing/duplicate resolutions
fail validation rather than becoming convenient abstentions.

## Executed evidence and limitations

Archived public inputs, production receipts, private truth provenance, proposals,
policies and reports are in
[`benchmark/v3/results/stage3-agent-evaluation-v2/`](../../benchmark/v3/results/stage3-agent-evaluation-v2/).
`validation.json` distinguishes successful execution from unavailable biological
validation; `sha256.json` covers the archived artifacts.

Fresh baseline, logistic and GBDT production inference succeeded on all seven
genomes: **21 genome/arm runs**. Their wall times were 105.790, 127.806 and 125.584
seconds respectively, or **359.180 seconds shared production cost**. Since every
method sees both learned scores, that shared cost belongs to every method;
incremental decision time is reported separately. Predictor generation, model
training, data acquisition and artifact verification are outside these timings.

| Split | Genome | Exported anomalies |
| --- | --- | ---: |
| calibration | GCA_000009085.1 | 1 |
| calibration | GCA_000027305.1 | 37 |
| calibration | GCA_000145295.1 | 0 |
| test | GCA_000007065.1 | 0 |
| test | GCA_000008765.1 | 1 |
| test | GCA_000009985.1 | 7 |
| test | GCA_000147055.1 | 150 |

There are **196 cases**, comprising 182 pathway gaps, ten borderline conflicts
and four borderline essential functions. One genome reaches the extractor's
150-case limit. These cases are not a representative sample of genomes or an
exhaustive sample of unresolved loci. Repeated pathway contexts can also produce
dependent case decisions; case counts must not be called independent recoveries.

The private join inspected 361 experimental assertions in these evaluation
genomes, but **none supports a selected candidate–function assignment**. Only one
selected candidate protein overlaps the labelled proteins at all. Thus neither
calibration precision nor test decision accuracy is identifiable here. All six
simple methods correctly produce an unsupported-threshold policy and abstain on
all 158 test cases. This is a diagnostic of insufficient validation labels, not
evidence that any method is accurate or superior.

The live `qwen3.8-27b` run also completed all 196 cases without request failures.
It made **392 HTTP requests**, with no retries, in **589.932 seconds** of decision
wall time. The endpoint reported **730,626 prompt tokens and 71,169 completion
tokens**. Before calibration it returned 193 abstentions and three demotions;
those are actions, not verified correct decisions. Its unsupported-threshold
policy abstains on all 158 test cases, just as the simple methods do. No monetary
cost is inferred without a billing rate. Shared production plus incremental
agent decision wall time is 949.112 seconds, with the timing exclusions above.

The models still lack a trained null layer, and this input set has no reaction
gaps or operons. It is the previously examined development cohort, not a new
blind test. The original independent-outcome, unseen-family, correction/harm and
complete-cost acceptance requirements remain in
[the acceptance record](../../tasks/v3-acceptance.md). The next empirical step
requires independently measured outcomes for production-selected unresolved
cases; selecting cases or manufacturing negatives from the desired answers is
not a substitute.

## Reproduction

Use the Python environment documented in `tasks/v3-acceptance.md`, set `TMPDIR`
under `$HOME/tmp`, and set `JAVA_TOOL_OPTIONS=-Djava.io.tmpdir=$TMPDIR`.

```bash
agent_data="$HOME/.cache/gspa-v3-data/agent-evaluation-v2"
python benchmark/v3/run_agent_evidence.py --out "$agent_data"
python benchmark/v3/evaluate_stage3_agent.py prepare \
  --sources "$agent_data/sources.json" --out "$agent_data/public-v2.json"
python benchmark/v3/prepare_agent_truth.py \
  --public "$agent_data/public-v2.json" \
  --labels benchmark/v3/results/training-cohort/labels.tsv \
  --go-obo "$HOME/.cache/gspa-v3-data/training-cohort/go-basic.obo" \
  --out "$agent_data/private-v2"
```

For each mode in `abstain_only`, `rank_by_raw_seq`,
`rank_by_stage1_posterior`, `rank_by_context_support`, `learned_logistic`, and
`learned_gbdt`, run `suggest --public ... --mode ... --go-obo ... --out ...`.
Run `agent` with the same public and ontology files for the live comparator.
Then run `calibrate --public ... --proposals ... --truth private-v2/calibration.json
--out policy.json`, followed by `evaluate --public ... --proposals ... --policy
policy.json --truth private-v2/test.json --out report.json`.

Choose fresh output paths for changed protocols; the decision, policy and report
writers deliberately refuse to overwrite existing artifacts. Production runs
can reuse only matching, verified receipts. Tests are in
`test_evaluate_stage3_agent.py`, `test_prepare_agent_truth.py` and the production
runner suites. Fixture tests prove contracts, not biological accuracy.
