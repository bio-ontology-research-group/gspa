# Locus/null training truth integrity

Status (2026-10-26): the preparation/training protocol is repaired and tested.
The replay described here has no complete, independently labelled locus groups.
A subsequent [GOA source audit](v3-gaf-negative-audit.md) produced one reviewed
two-candidate training competition in a separate augmented table; it has no
complete calibration/test or all-negative groups. **Items 4 and 7 remain open**;
no null model or biological validation is claimed by this work.

## Reproduced failures

The former builder retained a locus ID when an unlabelled member disappeared
during label joining. It also trusted a declaration that omitted another exported
candidate altogether. Either could train on a partial competition. Family masking
could similarly remove a positive while leaving its negatives, manufacturing an
apparently all-negative group. The trainer grouped solely by a user-supplied name,
so equal names in different genomes merged into one competition.

The label reader had a separate polarity bug: a supplied `label=0` was ignored
and read as positive truth. Five failing preparation regressions and a failing
cross-genome trainer regression were captured before the fixes. The final tests
also cover family masking, conflicting assays, explicit negative provenance,
GO negative scope, malformed certificates and direct trainer entry points.

## Label contract

Existing labels without a polarity field remain positive. Explicit negative
assays use `label=0` or `qualifier=NOT`; if both appear they must agree. Labels must
be exactly binary, not fractions or booleans coerced to integers. A negative row
needs a nonempty `reference` and the dataset/function policy needs a written
`explicit_negative_truth` source. For example, using **fixture identifiers only**:

```text
genome_id  protein_id  function_id  evidence  label  reference
g1         p1          EC:1.1.1.1  IDA       1      fixture:assay-A
g1         p2          EC:1.1.1.1  IDA       0      fixture:assay-B
```

Actual input is tab-separated. A source identifier must resolve to the actual
assay in a real dataset; writing a string is not proof of its validity. The same
evidence-code and role policies apply to positive and negative assertions.
Negative GO assertions stay exact, including primary-ID normalization; they are
not propagated to ancestors or other roles. Conflicting positive/negative truth,
including a positive propagated to an explicitly negated ancestor, fails.
Unsupported qualifiers also fail rather than being treated as positive.

Closed-world candidate labels remain available for explicitly declared marginal
reference-support models. They are stamped `closed_world_unconfirmed`, not
`explicit_negative_truth`. A protein with only a negative assay does not become
positively characterized for the purpose of inventing other closed-world labels.

## Complete competition contract

The production locus layer compares all supplied candidates for one function and
role within one genome. The declaration TSV has `genome_id`, `locus_group_id`,
`candidate_key`, and `truth_complete`. The policy's `locus_truth_complete` must be
a written source/justification, not a boolean. Duplicate candidates, conflicting
completeness flags, or groups mixing functions/roles fail.

The builder records the entire exported candidate universe **before** label
filtering. After label joining and family masking, a group is eligible only when:

- its declaration equals that complete exported universe;
- every declared member remains in the final labelled table;
- every negative is explicit assay truth, not merely unconfirmed.

Ineligible groups keep their useful marginal rows, lose locus supervision, and
retain the declared/exported/surviving memberships and exclusion reasons in
`locus_audit`. Eligible rows carry identical `locus_members` lists of example IDs
and a `locus_truth_source`. The trainer verifies certificates after loading and
before production training or evaluation. A later filter cannot silently turn a
complete group into a smaller one. Local group names are scoped by genome.
Multiple-positive groups still remain marginal examples and are excluded from
the one-candidate-or-null layer; isoenzymes are not forced into exclusivity.

Trainer version 2 hashes truth metadata along with train/calibration inputs.
Its artifacts record train/calibration label-basis counts and the probability
target. They never add test label summaries to the fitting artifact. Synthetic
mechanism tests remain explicitly synthetic; they cannot authorize a production
model with fabricated biological truth.

## Real-data regression replay

`benchmark/v3/replay_locus_truth_training.py` rebuilds the existing 18-genome
development table from its validated production exports, labels, policies and
family/genus split. It checks original input hashes, pins all feature files and
sources, uses the original fitting settings, freezes each refitted artifact
before test evaluation, and verifies inputs again afterwards. This is a replay
of previously used data, not a new blind benchmark.

All **4,262** original candidate identities, features, labels, and split memberships
remain unchanged: **541** positive rows and **3,721** closed-world unconfirmed rows.
These are candidate rows, including GO propagation, not independent assay counts.
There are **zero** explicit negative rows and **zero** declared complete locus
groups. Both logistic and tree models correctly retain `locus_model.trained=false`.
All scoring parameters, calibration values and thresholds equal their original
artifacts exactly. The new artifact hashes reflect the new truth metadata and
trainer protocol; their probability target is explicitly
`reference_annotation_support_under_declared_closed_world_policy`.

Both old/new artifact pairs were also loaded by the production CLI against the
same real claims and layout for `GCA_000007005.1`, using lite mode without priors.
For each model kind, all **26,328** prediction rows are byte-identical and all
learned candidate records agree. These four runs check artifact compatibility
and unchanged scoring, not biological performance under that reduced context.
Their commands and input/output hashes are retained in `cli-validation.json`.

Reproduce in a new output directory:

```bash
export TMPDIR="$HOME/tmp"
"$HOME/tmp/gspa-v3-testenv/bin/python" benchmark/v3/replay_locus_truth_training.py \
  --original-manifest "$HOME/.cache/gspa-v3-data/training-cohort/training_table/training_table_manifest.json" \
  --original-models benchmark/v3/results/training-cohort/learned_models \
  --out "$HOME/tmp/locus-truth-new-replay"
```

Ten new preparation regressions and six added trainer regressions pass. The full
benchmark suite ran **492 tests**, with three optional skips and no failures,
including eight live Nextflow tests. All **40 trainer tests** pass. The repaired
protocol changes no JVM implementation.

Inputs, source hashes, rebuilt table manifest, refitted models, test reports and
verification logs are archived in
[`benchmark/v3/results/locus-truth-v2/`](../../benchmark/v3/results/locus-truth-v2/).
Large feature exports and the 26 MB rebuilt table remain at the recorded cache
paths. Actual experimental negative/role acquisition, complete locus labels,
calibrated null fitting and supported-predictor transfer validation remain open.
