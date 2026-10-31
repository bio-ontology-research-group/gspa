# Native CLEAN, DIAMOND and context comparison

Status: the EC identity and claim-serialization corrections pass cluster
regression tests and real-genome probe 10403 completes with exit zero. Runtime
preflight 10400 is complete. Full-cohort source/context exports, refits and accuracy
results await the complete native CLEAN acquisition audit. Job 10406 is
queued behind both prerequisites; no full-comparison result is available. This extends the completed
[DIAMOND context comparison](v3-context-comparison-results.md), whose default
control still contained context. No model is promoted and no whole acceptance
item is complete.
An interim conditional evaluation of the 10406 fits is reported in the penultimate
section; it is not the primary endpoint and is superseded where it overlaps by the
last section, the verified full production replay.

Probe 10389 supplied the same 1,244 protein identifiers to both evidence sources,
but DIAMOND used bare EC numbers and CLEAN used `EC:` prefixes. GSPA retained
511 shared protein/EC pairs as separate candidates. The combined sequence-only
export therefore contained 30,023 candidates; resolving prefix equivalence gives
29,512. This is an identity defect, not experimental evidence that either
prediction is correct. The old export audit checked literal strings and missed
it. Its passing results are insufficient for the combined-source comparison.

The correction now normalizes EC spelling at claim ingestion, with unchanged raw
source files and explicit original-ID provenance. Python checks expect canonical
`EC:` identifiers and reject split duplicate exports. Supervision is compared by
genome, protein, function and role, with unchanged labels, references, locus groups
and splits; spelling changes are counted separately. All 6,706 original supervision
rows have distinct canonical identities. No duplicate labels are discarded.
A second correction preserves nested predictor contracts when claims are written
and reloaded; numerical feature equality alone had missed that provenance loss.

Job 10402 passes 30 focused JVM cases and the full 1,473-case suite (six skips,
no failures). All 105 focused Python checks pass. The tested JAR is
`/data/hohndor/gspa-v3/ec-prefix-identity-v3/jobs/10402/gspa-3.0.0.jar`, SHA-256
`37aeee2b565951dd3a4e049a23aacd2e8155f29bd8b6216fed7006e4273965e6`.
These are integration checks, not evidence of improved enzyme accuracy.

The new design has five sources: DIAMOND homology, DIAMOND with family exclusion,
native CLEAN, and CLEAN combined with each DIAMOND variant. Each receives four
controls: sequence only, UNIVERSAL, PATHWAY and ORGANISM context. Logistic and
tree fits give 40 models, followed by 720 genome/model production runs. The
sequence-only control removes actual context inputs and disables priors before
the JVM exports features. Every `ctx.*`, `cohort.*`, `asm.*` and `fba.*` value must
then be zero; columns are not masked after extraction. The executed-stage ledger
must also show no context stages. Function identity, GO complex classification,
and declared taxon/profile metadata remain available in every control.

The runner retains all 18 original genomes and 49,463 proteins, including CLEAN's
explicit unsupported-input abstention. Source files are concatenated unchanged,
with candidate-specific source and score contracts. The entire original
supervision, genome splits and train/calibration settings remain fixed, while
the candidate tables may legitimately change with the source. A genome without
label-matched candidates is still part of the full reference denominator. All
40 models must be frozen before test evaluation. An unexplained extra exported
candidate is an error, not a reason to drop it or assign it sequence evidence.

Production replay checks the actual CLI predictions against the independent
Python model reader. It compares candidates by identity, preserves distinct
roles, verifies feature digests and checks projected probabilities. Unsupported
scores must preserve the preceding output. Additions, removal of annotations
without a projectable role, and structural taxon restrictions remain explicit;
this comparison does not request threshold-based retraction. These are execution
and numerical checks, not evidence of correct biological assignments.
Empty and wholly unsupported genomes stay in the actual joint JVM batch; their
completed learned report and unchanged baseline output must verify. Large learned
reports are streamed by the receipt checker instead of failing an 8 MB size gate.

The primary endpoint uses every experimental-reference EC assertion in the
registered test genomes, including functions never proposed. The existing test
population has 39 assertions on 37 proteins. Sixty final-output arms (five
sources, four contexts, three scorers) are compared with native CLEAN's unchanged
maxsep confidence. The primary cutoff maximizes reference F1 across all positive
EC assertions in the calibration genomes, including unproposed functions. A
mixed candidate-table cutoff is secondary, alongside the common 0.5 cutoff.
All 61 arms' calibration thresholds are saved before test metrics. CLEAN's native
operating point additionally includes every maxsep-emitted call, including zero
confidence, and receives direct paired comparisons against the other arms.
Reports include reference coverage, ranking, matched and unconfirmed calls,
missing references, and paired gains and losses by genome. Missing predictions
are explicit no-calls even at threshold zero. Conditional candidate-table
metrics use each arm's own candidate population and cannot establish matched
cross-source accuracy. Exact EC identities are used; obsolete identifiers are
not silently merged.

This is a reused development cohort, not an untouched sequence test. CLEAN has
3,730 exact reference-sequence overlaps in these inputs; combining it with
family-excluded DIAMOND does not establish unseen-family performance. Positive
references are incomplete, and unmatched calls are unconfirmed rather than
experimentally disproven. The small, clustered EC test population also limits
claims of general improvement. The full native acquisition and independent
collection audit remain prerequisites; see the [CLEAN acquisition record](v3-clean-cohort-acquisition.md).

Implementation:

- `benchmark/v3/run_source_arm_context_comparison.py`: prepare, export, train and conditional evaluation.
- `benchmark/v3/check_source_arm_context_models.py`: complete production replay and score/role audit.
- `benchmark/v3/analyze_source_arm_context_models.py`: fixed-reference EC comparison and native CLEAN control.

The current 37 source/context checks use controlled fixtures and mocked training
calls. Alongside 22 source-input and 46 shared-runner/analysis checks, they test
source combinations, context removal, changed inputs, supervision and split
preservation, all 40 model plans, unsupported fallback, role projection, no-call
handling, full reference denominators and threshold freezing. No model was fitted
on the laptop. Real production probes and biological comparisons remain separate.

Probe history is retained rather than replacing failed attempts:

| Job | Outcome and scope |
|---|---|
| 10379 | Passed 30 checks; failed before inference because the wrapper named the older parent-cohort directory. Its never-run dependent 10380 was cancelled. |
| 10385 | Passed 35 checks and verified all 26,328 DIAMOND sequence-only candidates with zero context features; then rejected metadata for 18 genomes supplied with a one-genome manifest. |
| 10389 | All nine exports and both historical numerical model audits completed, but the literal-key audit missed the 511 duplicate EC identities. This is insufficient combined-source validation. |
| 10403 | Completes all nine exports and both historical numerical audits with exit zero. Confirms 29,512 combined candidates and both source contracts on all 511 shared EC assignments. |

The one-genome probe binds the original metadata and an unchanged copy of the
selected genome's entry, with the full pathway definitions. It compares DIAMOND,
native CLEAN and their combination under sequence-only, context and RLC-disabled
controls. It does not test cross-genome recurrence or biological accuracy. The
full benchmark retains the original 18-genome population.

The remote execution wrapper is `benchmark/v3/source_context_comparison.sbatch`,
with stage orchestration in `benchmark/v3/run_source_context_remote.py`. It
restores the original Python 3.10.19/numpy 2.2.6 runtime and the seven supplemental
supervision files, keeping the original registered source paths intact. Current
comparison code and the pinned streaming parser occupy a separate directory.
The wrapper preserves every attempt's results and stage logs, including failures.
Full mode requires the complete native CLEAN audit and a successful three-source
production probe before preparing any comparison arm. It then executes all 20
exports, 40 fits, 720 learned genome runs and the 61-arm EC analysis in order.

Runtime preflight 10392 passed 73 checks and stopped because the bundle omitted
an unchanged legacy-analysis helper. Preflight 10397 passed 81 checks and found
ten omitted parent-template files. Their original registered hashes match the
recovered files (45.9 MB). Fresh immutable preflight 10400 includes them, passes
its 104 packaged checks and completes parent-chain and runtime verification with
exit zero. The current wrapper suite has 105 cases after an additional raw-ID
provenance-conflict regression. No preflight fits models or evaluates accuracy.

Full execution requires a fresh immutable bundle binding the tested JAR, the
completed CLEAN audit and the independently checked schema-2 probe receipt.
It rechecks probe outputs and shared pipeline sources, then requires every source,
context, model and genome in final population checks. The native emitted-call
endpoint retains threshold zero; each EC panel retains all 39 assertions on 37
proteins in four test genomes. Completed stages receive durable partial
checkpoints. Only final verification and successful publication certify execution.

The [frozen launch checkpoint](../../benchmark/v3/results/source-context-launch-v1/)
records job 10406 at `/data/hohndor/gspa-v3/source-context-full-v1`. It requests
12 CPUs and 48 GiB on node004, with a two-day limit, no requeue, and successful
job dependencies on corrected probe 10403 and independent CLEAN audit 10368.
A failed dependency prevents inference. After both jobs finish, the launcher
requires their recorded exit codes to be zero, the 511 merged EC identities and
29,512 combined candidates, and CLEAN's complete 18-genome/49,463-protein
accounting. It exclusively records the prerequisite hashes before calling the
full runner, which independently checks their contents and all bound outputs.

All 60 bundled files are frozen now; the 37 files shared with the probe are
byte-identical. Seven focused checks cover failed or missing prerequisites,
incorrect merged counts, a reduced CLEAN population, changed or duplicate source
inputs and attempted overwriting of an existing sealed run. Sonnet reviewed the
launcher, sealing and downstream execution without finding an implementation
issue. The queued job is scheduling progress, not completed fitting, evaluation,
model promotion or acceptance of any whole requirement.

Corrected probe 10403 is now complete. All three source arms pass sequence-only,
UNIVERSAL and RLC-disabled exports; the sequence-only features have all context
values zero and only the intended stages execute. The combined arm retains all
511 shared EC pairs with both sources and contracts, producing 29,512 candidates.
Both historical learned-output reader audits agree to below 7e-16 probability
error, retaining 7,700 unsupported candidates each. These historical checks do
not refit a model. [The probe checkpoint](../../benchmark/v3/results/multisource-probe-corrected-v1/)
retains its successful exit, complete receipt, logs and exact source snapshot.
Independent job 10411 verifies all 145 published output hashes and the first two
feature exports, then fails because its wrapper assumes the wrong path for the
RLC-disabled export. Fresh audit 10413 uses the nine paths confirmed against the
published receipt and recomputes source identities, contracts and score features
with a separate stdlib reader. Job 10413 completes with exit zero: all 145
published hashes and 178,605 feature rows pass, including all 511 combined-source
EC identities under each control. No full-cohort accuracy or cross-genome
recurrence conclusion follows from this one-genome probe.

## Conditional candidate-table evaluation of the 10406 fits (interim, 2026-10-29)

> Superseded where it overlaps by the section "Full production replay of the
> 10406 fits (verified, Shaheen 2026-10-29)" below. The text is kept as written;
> statements that the replay and the reference analysis are pending no longer
> hold.

Job 10406 completed preparation, the 20 exports, the 40 model fits and the
conditional evaluation, then stopped in production replay on a display-threshold
tie; the replay is being resumed (see the
[supervision log](../../tasks/v3-supervision-plan-2026-10-28.md)). This section reports only the
conditional evaluation, which ranks the candidates each arm exported and scores
them against the reference labels after the models were frozen. As recorded in
the evaluation itself, it is conditional on the exported candidates, uses a
reused cohort that overlaps the CLEAN reference (exact known sequences, not
unseen families), and supports no claim of biological accuracy, calibration or
promotion. Production replay and the fixed-reference EC analysis are pending.

Tables are in
[`benchmark/v3/results/source-context-v1/evaluation-10406/`](../../benchmark/v3/results/source-context-v1/evaluation-10406/)
(`source_context_evaluation.md` and `.tsv`, tabulated from the unchanged
`evaluation.json` by `benchmark/v3/summarize_source_context_evaluation.py`). Each cell
below gives AUROC / log-loss (supported test rows / all test rows); abstained rows
are excluded from both metrics.

| source | kind | SEQONLY | UNIVERSAL | PATHWAY | ORGANISM |
|---|---|---|---|---|---|
| `diamond_homology` | logistic | 0.723 / 0.432 (463/608) | 0.761 / 0.402 (462/608) | 0.759 / 0.413 (461/608) | 0.766 / 0.400 (460/608) |
| `diamond_homology` | gbdt | 0.688 / 0.462 (463/608) | 0.719 / 0.424 (462/608) | 0.710 / 0.432 (461/608) | 0.715 / 0.426 (460/608) |
| `diamond_family_excluded` | logistic | 0.682 / 0.344 (300/416) | 0.751 / 0.309 (299/416) | 0.737 / 0.300 (297/416) | 0.752 / 0.307 (299/416) |
| `diamond_family_excluded` | gbdt | 0.713 / 0.329 (300/416) | 0.745 / 0.311 (299/416) | 0.737 / 0.302 (297/416) | 0.756 / 0.305 (299/416) |
| `clean` | logistic | 0.768 / 0.755 (38/38) | 0.042 / 0.707 (36/38) | 0.022 / 0.722 (35/38) | 0.032 / 0.656 (30/38) |
| `clean` | gbdt | 0.525 / 0.744 (38/38) | 0.628 / 0.727 (36/38) | 0.536 / 0.866 (35/38) | 0.582 / 0.725 (30/38) |
| `clean_diamond_homology` | logistic | 0.778 / 0.385 (483/631) | 0.776 / 0.366 (459/631) | 0.768 / 0.366 (461/631) | 0.759 / 0.369 (456/631) |
| `clean_diamond_homology` | gbdt | 0.776 / 0.392 (483/631) | 0.741 / 0.416 (459/631) | 0.723 / 0.413 (461/631) | 0.713 / 0.413 (456/631) |
| `clean_diamond_family_excluded` | logistic | 0.808 / 0.329 (333/449) | 0.809 / 0.305 (317/449) | 0.801 / 0.310 (312/449) | 0.797 / 0.309 (311/449) |
| `clean_diamond_family_excluded` | gbdt | 0.830 / 0.330 (333/449) | 0.814 / 0.350 (317/449) | 0.814 / 0.339 (312/449) | 0.818 / 0.331 (311/449) |

The same results as context minus SEQONLY, each side on its own supported rows
(U = UNIVERSAL, P = PATHWAY, O = ORGANISM; negative log-loss differences are
improvements):

| source | kind | ΔAUROC U / P / O | Δlog-loss U / P / O |
|---|---|---|---|
| `diamond_homology` | logistic | +0.038 / +0.036 / +0.043 | −0.030 / −0.019 / −0.032 |
| `diamond_homology` | gbdt | +0.031 / +0.022 / +0.027 | −0.039 / −0.031 / −0.036 |
| `diamond_family_excluded` | logistic | +0.069 / +0.055 / +0.070 | −0.035 / −0.043 / −0.037 |
| `diamond_family_excluded` | gbdt | +0.031 / +0.024 / +0.043 | −0.018 / −0.027 / −0.024 |
| `clean_diamond_homology` | logistic | −0.001 / −0.009 / −0.018 | −0.019 / −0.018 / −0.016 |
| `clean_diamond_homology` | gbdt | −0.035 / −0.053 / −0.063 | +0.024 / +0.021 / +0.021 |
| `clean_diamond_family_excluded` | logistic | +0.001 / −0.006 / −0.011 | −0.024 / −0.019 / −0.020 |
| `clean_diamond_family_excluded` | gbdt | −0.016 / −0.016 / −0.012 | +0.020 / +0.009 / +0.001 |
| `clean` (not comparable) | logistic | −0.726 / −0.746 / −0.736 | −0.048 / −0.033 / −0.099 |
| `clean` (not comparable) | gbdt | +0.103 / +0.011 / +0.057 | −0.017 / +0.122 / −0.018 |

What the tables support is limited. The counts are tiny: the positive labels per
cell number a few dozen, the evaluation supplies no intervals and no paired test,
and differences of a few hundredths of AUROC lie within plausible sampling noise.

- **DIAMOND-only arms.** Every learned context (UNIVERSAL, PATHWAY, ORGANISM)
  gives a higher AUROC and a lower log-loss than SEQONLY under both model kinds,
  for both DIAMOND variants. The direction is consistent and the gains are small.
- **CLEAN+DIAMOND arms.** Context gives no AUROC gain. With gbdt, AUROC is lower
  than SEQONLY in every context, and log-loss is also higher. With logistic
  models AUROC is flat to slightly lower while log-loss is slightly lower. For
  these arms this is a negative result for context. The three learned contexts
  differ little from each other in any arm.
- **Rows are not paired across contexts.** Abstention changes with context, so
  each context is evaluated on a different set of rows, and in the CLEAN+DIAMOND
  arms the abstentions remove positives disproportionately. A difference
  therefore mixes a model effect with a change of the evaluated population.
  Display thresholds are also fitted per model, so precision and recall (in the
  TSV and Markdown tables) come from different operating points.
- **No locus model was trained in any of the 40 fits.** The CLEAN arm has no rows
  carrying a locus group; in the other four arms the training split has no
  no-correct-candidate groups and only one single-positive group, short of the
  minimum for an honest null outcome. The `locus` fields carry no result.
- **The `clean` arm is not interpretable.** It has only 38 test rows and a
  positive rate several times that of the other arms, so it is a different
  population as well as a smaller one. Its near-zero logistic AUROC under the
  three learned contexts is an artefact, not a context effect: the supervisor
  inspected the frozen model artefacts and found that the Platt calibration slope of
  these three logistic models is a = 0.001 (SEQONLY: a = 0.93), fitted on a
  65-row calibration set. The raw model ranks the test rows inversely, and the
  calibration collapses to a near-constant output, a degenerate fit. This arm
  should not be cited as evidence for or against context.

The final comparison, after the resumed replay and analysis, will add:

- coverage over the full reference denominator, including functions never proposed;
- paired corrections and harms, on the rows shared across contexts;
- the experimental EC reference panel with native CLEAN's unchanged confidence as the control.

## Full production replay of the 10406 fits (verified, Shaheen 2026-10-29)

This section supersedes the interim section above wherever they overlap, and
completes the three items that section left open: coverage over the full
reference denominator, paired changes against named controls, and the
experimental EC reference with native CLEAN's unchanged confidence as control.
The interim section is kept unchanged as the record of what was known before the
replay. Its conditional-table results (AUROC, log-loss) are not repeated here.

**What was run.** Job 10406 had completed the exports, the 40 model fits and the
conditional evaluation, then stopped in production replay on a display-threshold
tie. The replay was resumed on Shaheen in jobs 15691804 (24 arms) and 15692096
(the remainder and the analysis), reusing the 40 frozen fits without refitting.
In total the run covers 61 analysis arms (5 sources x 4 contexts x 3 kinds, plus
native CLEAN), 20 exports, 40 fitted models and 720 replayed genome cells (5 sources x 4
contexts x 2 learned kinds x 18 genomes). The final receipt reports status
verified and `model_promotion: none`. One source drift is declared: the model
checker `check_source_arm_context_models.py` has a different hash at
analysis time than the one registered with the probe, and the receipt records both
hashes instead of hiding the change.

**What the replay verifies.** It verifies execution, input preservation and
agreement between the Python and JVM scores. It does not verify biological
accuracy, and no table below is evidence of accuracy beyond the stated
reference. Additions and removals are explicit: no replay added or removed a
function. Display-threshold ties, the failure that stopped 10406, are now
recorded per cell instead of aborting the run.

**What the analysis scores.** Every arm is scored against the same fixed
population: all positive exact EC assertions (ECO:0000269) on the reference
proteins of the registered test genomes, including proteins for which the arm
emits nothing. The EC identifiers are compared exactly, with no merging of
obsolete numbers. A call that matches no reference assertion is *unconfirmed*;
it is not a negative, and a missing reference assertion is not biological
absence. The reference is small: all arms share one denominator, shown in the
tables (37 proteins, 39 assertions, one of the four test genomes contributing
none). The `clean` source rests on these same 37 proteins and, in the
conditional candidate table of the interim section, on 38 rows, so any
difference of a few calls in it is within what one protein can change. The
cutoff per arm comes from three operating points: `calibration_selected` (the
frozen calibration cutoff of each model), `ec_calibration_selected` (a cutoff
chosen on the calibration genomes' EC reference) and `fixed_0.5`. Native CLEAN
is shown at each of those points and, separately, with every emitted call kept
(`native_emitted_calls`).

### Replay integrity

Per source, summed over 4 contexts x 2 learned kinds x 18 genomes (candidates
are replayed candidate functions; error is the largest absolute Python-JVM
difference in probability):

| source | genomes | cells | candidates | rescored | scored only | unchanged unsupported | added | removed | display ties | max abs Python-JVM error |
|---|---|---|---|---|---|---|---|---|---|---|
| `diamond_homology` | 18 | 8 | 3873616 | 2716922 | 38488 | 1118206 | 0 | 0 | 1053 | 8.88e-16 |
| `diamond_family_excluded` | 18 | 8 | 2178928 | 1548300 | 21838 | 608790 | 0 | 0 | 2913 | 5.00e-16 |
| `clean` | 18 | 8 | 482400 | 358816 | 0 | 123584 | 0 | 0 | 0 | 2.22e-16 |
| `clean_diamond_homology` | 18 | 8 | 4274200 | 2935700 | 38488 | 1300012 | 0 | 0 | 246 | 1.44e-15 |
| `clean_diamond_family_excluded` | 18 | 8 | 2623304 | 1805110 | 21838 | 796356 | 0 | 0 | 0 | 1.55e-15 |

Candidates left unchanged because the learned stage could not score them, by reason:

| source | `profile_lacks_both_classes_in_train_and_calibration` | `sequence_score_outside_training_calibration_range` | `unseen_source_function_role_or_availability_profile` |
|---|---|---|---|
| `diamond_homology` | 1051018 | 16388 | 50800 |
| `diamond_family_excluded` | 568978 | 11566 | 28246 |
| `clean` | 43406 | 79634 | 544 |
| `clean_diamond_homology` | 1160854 | 78100 | 61058 |
| `clean_diamond_family_excluded` | 687498 | 71056 | 37802 |

The agreement between Python and JVM scores is at floating-point rounding for
every source. In every source a quarter to a third of the candidates are left
unchanged because the model has no support for them, so the learned arms
rescore only part of the candidate universe; the unchanged candidates keep the
preceding annotation and are counted in the analysis below like any other call.
For the `clean` source the unsupported share is zero under SEQONLY and
substantial under every learned context, and the dominant reason differs by
context (calibration range for PATHWAY, missing profile classes for ORGANISM), so contexts do not
act on the same subset of candidates.

### Native CLEAN baselines on the full EC reference

| native CLEAN at operating point | cutoff | proteins | reference pairs | matched | precision | recall | micro-F1 | top-1 | unconfirmed calls |
|---|---|---|---|---|---|---|---|---|---|
| `calibration_selected` | 0.8519 | 37 | 39 | 18 | 0.692 | 0.462 | 0.554 | 24/37 | 8 |
| `ec_calibration_selected` | 0.8519 | 37 | 39 | 18 | 0.692 | 0.462 | 0.554 | 24/37 | 8 |
| `fixed_0.5` | 0.5 | 37 | 39 | 19 | 0.679 | 0.487 | 0.567 | 24/37 | 9 |
| `native_emitted_calls` | 0.0 | 37 | 39 | 24 | 0.632 | 0.615 | 0.623 | 24/37 | 14 |

Thresholding native CLEAN trades recall for precision: with every
emitted call kept it recovers the most reference assertions, and any calibration cutoff
removes some of them. All arms below are compared against these rows. The
top-1 count is identical across the four rows because the thresholds only
remove lower-ranked calls.

### Full EC reference per arm

Micro-F1 over the 39 reference assertions; the number of recovered
assertions is in parentheses. `rule` is the unlearned arm for that source (the
source's own confidence) and does not depend on context. `precision` is
undefined when an arm emits no call; those cells count as no recovery.

**`calibration_selected`**, micro-F1 (matched reference pairs of 39):

| source | kind | SEQONLY | UNIVERSAL | PATHWAY | ORGANISM |
|---|---|---|---|---|---|
| `diamond_homology` | rule | 0.276 (8) | 0.276 (8) | 0.276 (8) | 0.276 (8) |
| `diamond_homology` | logistic | 0.328 (11) | 0.379 (11) | 0.379 (11) | 0.250 (6) |
| `diamond_homology` | gbdt | 0.182 (5) | 0.345 (10) | 0.269 (7) | 0.373 (11) |
| `diamond_family_excluded` | rule | 0.043 (1) | 0.043 (1) | 0.043 (1) | 0.043 (1) |
| `diamond_family_excluded` | logistic | 0.111 (3) | 0.085 (2) | 0.045 (1) | 0.093 (2) |
| `diamond_family_excluded` | gbdt | 0.048 (1) | 0.204 (5) | 0.128 (3) | 0.174 (4) |
| `clean` | rule | 0.554 (18) | 0.554 (18) | 0.554 (18) | 0.554 (18) |
| `clean` | logistic | 0.554 (18) | 0.640 (24) | 0.583 (21) | 0.523 (17) |
| `clean` | gbdt | 0.540 (17) | 0.585 (19) | 0.562 (18) | 0.571 (18) |
| `clean_diamond_homology` | rule | 0.562 (18) | 0.562 (18) | 0.562 (18) | 0.562 (18) |
| `clean_diamond_homology` | logistic | 0.667 (24) | 0.644 (19) | 0.545 (15) | 0.610 (18) |
| `clean_diamond_homology` | gbdt | 0.523 (17) | 0.545 (18) | 0.459 (14) | 0.537 (18) |
| `clean_diamond_family_excluded` | rule | 0.562 (18) | 0.562 (18) | 0.562 (18) | 0.562 (18) |
| `clean_diamond_family_excluded` | logistic | 0.609 (21) | 0.621 (18) | 0.571 (16) | 0.618 (21) |
| `clean_diamond_family_excluded` | gbdt | 0.600 (21) | 0.522 (18) | 0.476 (15) | 0.537 (18) |

**`ec_calibration_selected`**, micro-F1 (matched reference pairs of 39):

| source | kind | SEQONLY | UNIVERSAL | PATHWAY | ORGANISM |
|---|---|---|---|---|---|
| `diamond_homology` | rule | 0.281 (9) | 0.281 (9) | 0.281 (9) | 0.281 (9) |
| `diamond_homology` | logistic | 0.281 (9) | 0.379 (11) | 0.386 (11) | 0.174 (4) |
| `diamond_homology` | gbdt | 0.258 (8) | 0.369 (12) | 0.381 (12) | 0.418 (14) |
| `diamond_family_excluded` | rule | 0.042 (1) | 0.042 (1) | 0.042 (1) | 0.042 (1) |
| `diamond_family_excluded` | logistic | 0.043 (1) | 0.164 (5) | 0.164 (5) | 0.164 (5) |
| `diamond_family_excluded` | gbdt | 0.172 (5) | 0.185 (5) | 0.189 (5) | 0.189 (5) |
| `clean` | rule | 0.554 (18) | 0.554 (18) | 0.554 (18) | 0.554 (18) |
| `clean` | logistic | 0.554 (18) | 0.640 (24) | 0.583 (21) | 0.640 (24) |
| `clean` | gbdt | 0.576 (19) | 0.585 (19) | 0.562 (18) | 0.639 (23) |
| `clean_diamond_homology` | rule | 0.630 (23) | 0.630 (23) | 0.630 (23) | 0.630 (23) |
| `clean_diamond_homology` | logistic | 0.625 (25) | 0.600 (24) | 0.545 (15) | 0.622 (23) |
| `clean_diamond_homology` | gbdt | 0.608 (24) | 0.545 (18) | 0.459 (14) | 0.605 (23) |
| `clean_diamond_family_excluded` | rule | 0.622 (23) | 0.622 (23) | 0.622 (23) | 0.622 (23) |
| `clean_diamond_family_excluded` | logistic | 0.609 (21) | 0.621 (18) | 0.571 (16) | 0.629 (22) |
| `clean_diamond_family_excluded` | gbdt | 0.600 (21) | 0.522 (18) | 0.476 (15) | 0.537 (18) |

**`fixed_0.5`**, micro-F1 (matched reference pairs of 39):

| source | kind | SEQONLY | UNIVERSAL | PATHWAY | ORGANISM |
|---|---|---|---|---|---|
| `diamond_homology` | rule | 0.000 (0) | 0.000 (0) | 0.000 (0) | 0.000 (0) |
| `diamond_homology` | logistic | 0.250 (6) | 0.174 (4) | 0.136 (3) | 0.091 (2) |
| `diamond_homology` | gbdt | 0.192 (5) | 0.000 (0) | 0.049 (1) | 0.000 (0) |
| `diamond_family_excluded` | rule | 0.000 (0) | 0.000 (0) | 0.000 (0) | 0.000 (0) |
| `diamond_family_excluded` | logistic | 0.000 (0) | 0.000 (0) | 0.000 (0) | 0.049 (1) |
| `diamond_family_excluded` | gbdt | 0.000 (0) | 0.000 (0) | 0.000 (0) | 0.000 (0) |
| `clean` | rule | 0.554 (18) | 0.554 (18) | 0.554 (18) | 0.554 (18) |
| `clean` | logistic | 0.567 (19) | 0.640 (24) | 0.640 (24) | 0.639 (23) |
| `clean` | gbdt | 0.576 (19) | 0.620 (22) | 0.600 (21) | 0.588 (20) |
| `clean_diamond_homology` | rule | 0.554 (18) | 0.554 (18) | 0.554 (18) | 0.554 (18) |
| `clean_diamond_homology` | logistic | 0.667 (24) | 0.644 (19) | 0.633 (19) | 0.590 (18) |
| `clean_diamond_homology` | gbdt | 0.629 (22) | 0.567 (19) | 0.559 (19) | 0.537 (18) |
| `clean_diamond_family_excluded` | rule | 0.554 (18) | 0.554 (18) | 0.554 (18) | 0.554 (18) |
| `clean_diamond_family_excluded` | logistic | 0.648 (23) | 0.644 (19) | 0.623 (19) | 0.590 (18) |
| `clean_diamond_family_excluded` | gbdt | 0.592 (21) | 0.529 (18) | 0.559 (19) | 0.537 (18) |

Reading across sources: the DIAMOND-only arms are far below native CLEAN
on this reference, and the family-excluded variant recovers almost nothing.
That is expected if the reference proteins' families are removed from the
DIAMOND search, but this section does not test that explanation. The
`rule` rows of the DIAMOND-only sources emit no call at 0.5, which implies the
rule confidence of DIAMOND hits lies below 0.5; at that operating point a
learned arm that recovers anything wins by construction, and the comparison
says more about the scale of the rule score than about the learned model. The
CLEAN-containing arms stay close to native CLEAN and none of them clearly
exceeds native CLEAN with every call kept (see below).

### Paired comparisons

Each cell is the number of exact reference EC calls the first arm gains and loses
relative to the second, pooled over the same 37 proteins (per-genome counts are
in the TSV). Gains and losses are not netted in the TSV; a cell with gains and
losses can hide a larger turnover than the net suggests.

**`calibration_selected`**. Learned minus rule, same source and context (exact reference EC calls gained / lost, over 37 proteins):

| source | kind | SEQONLY | UNIVERSAL | PATHWAY | ORGANISM |
|---|---|---|---|---|---|
| `diamond_homology` | logistic | +3 / −0 | +3 / −0 | +4 / −1 | +3 / −5 |
| `diamond_homology` | gbdt | +0 / −3 | +3 / −1 | +2 / −3 | +4 / −1 |
| `diamond_family_excluded` | logistic | +2 / −0 | +1 / −0 | +0 / −0 | +1 / −0 |
| `diamond_family_excluded` | gbdt | +0 / −0 | +4 / −0 | +2 / −0 | +3 / −0 |
| `clean` | logistic | +0 / −0 | +6 / −0 | +4 / −1 | +1 / −2 |
| `clean` | gbdt | +0 / −1 | +1 / −0 | +1 / −1 | +2 / −2 |
| `clean_diamond_homology` | logistic | +6 / −0 | +1 / −0 | +1 / −4 | +0 / −0 |
| `clean_diamond_homology` | gbdt | +2 / −3 | +1 / −1 | +0 / −4 | +0 / −0 |
| `clean_diamond_family_excluded` | logistic | +3 / −0 | +1 / −1 | +1 / −3 | +3 / −0 |
| `clean_diamond_family_excluded` | gbdt | +3 / −0 | +1 / −1 | +1 / −4 | +0 / −0 |

**`calibration_selected`**. Context minus SEQONLY, same source and kind (gained / lost, 37 proteins):

| source | kind | UNIVERSAL | PATHWAY | ORGANISM |
|---|---|---|---|---|
| `diamond_homology` | rule | +0 / −0 | +0 / −0 | +0 / −0 |
| `diamond_homology` | logistic | +1 / −1 | +2 / −2 | +1 / −6 |
| `diamond_homology` | gbdt | +5 / −0 | +2 / −0 | +6 / −0 |
| `diamond_family_excluded` | rule | +0 / −0 | +0 / −0 | +0 / −0 |
| `diamond_family_excluded` | logistic | +1 / −2 | +0 / −2 | +1 / −2 |
| `diamond_family_excluded` | gbdt | +4 / −0 | +2 / −0 | +3 / −0 |
| `clean` | rule | +0 / −0 | +0 / −0 | +0 / −0 |
| `clean` | logistic | +6 / −0 | +4 / −1 | +1 / −2 |
| `clean` | gbdt | +2 / −0 | +2 / −1 | +3 / −2 |
| `clean_diamond_homology` | rule | +0 / −0 | +0 / −0 | +0 / −0 |
| `clean_diamond_homology` | logistic | +0 / −5 | +0 / −9 | +0 / −6 |
| `clean_diamond_homology` | gbdt | +3 / −2 | +0 / −3 | +3 / −2 |
| `clean_diamond_family_excluded` | rule | +0 / −0 | +0 / −0 | +0 / −0 |
| `clean_diamond_family_excluded` | logistic | +0 / −3 | +0 / −5 | +0 / −0 |
| `clean_diamond_family_excluded` | gbdt | +0 / −3 | +0 / −6 | +0 / −3 |

**`calibration_selected`**. Arm minus native CLEAN at the same cutoff (gained / lost, 37 proteins):

| source | kind | SEQONLY | UNIVERSAL | PATHWAY | ORGANISM |
|---|---|---|---|---|---|
| `diamond_homology` | rule | +2 / −12 | +2 / −12 | +2 / −12 | +2 / −12 |
| `diamond_homology` | logistic | +3 / −10 | +4 / −11 | +5 / −12 | +1 / −13 |
| `diamond_homology` | gbdt | +1 / −14 | +3 / −11 | +2 / −13 | +4 / −11 |
| `diamond_family_excluded` | rule | +0 / −17 | +0 / −17 | +0 / −17 | +0 / −17 |
| `diamond_family_excluded` | logistic | +2 / −17 | +1 / −17 | +0 / −17 | +1 / −17 |
| `diamond_family_excluded` | gbdt | +0 / −17 | +4 / −17 | +2 / −17 | +3 / −17 |
| `clean` | rule | +0 / −0 | +0 / −0 | +0 / −0 | +0 / −0 |
| `clean` | logistic | +0 / −0 | +6 / −0 | +4 / −1 | +1 / −2 |
| `clean` | gbdt | +0 / −1 | +1 / −0 | +1 / −1 | +2 / −2 |
| `clean_diamond_homology` | rule | +0 / −0 | +0 / −0 | +0 / −0 | +0 / −0 |
| `clean_diamond_homology` | logistic | +6 / −0 | +1 / −0 | +1 / −4 | +0 / −0 |
| `clean_diamond_homology` | gbdt | +2 / −3 | +1 / −1 | +0 / −4 | +0 / −0 |
| `clean_diamond_family_excluded` | rule | +0 / −0 | +0 / −0 | +0 / −0 | +0 / −0 |
| `clean_diamond_family_excluded` | logistic | +3 / −0 | +1 / −1 | +1 / −3 | +3 / −0 |
| `clean_diamond_family_excluded` | gbdt | +3 / −0 | +1 / −1 | +1 / −4 | +0 / −0 |

**`fixed_0.5`**. Learned minus rule, same source and context (exact reference EC calls gained / lost, over 37 proteins):

| source | kind | SEQONLY | UNIVERSAL | PATHWAY | ORGANISM |
|---|---|---|---|---|---|
| `diamond_homology` | logistic | +6 / −0 | +4 / −0 | +3 / −0 | +2 / −0 |
| `diamond_homology` | gbdt | +5 / −0 | +0 / −0 | +1 / −0 | +0 / −0 |
| `diamond_family_excluded` | logistic | +0 / −0 | +0 / −0 | +0 / −0 | +1 / −0 |
| `diamond_family_excluded` | gbdt | +0 / −0 | +0 / −0 | +0 / −0 | +0 / −0 |
| `clean` | logistic | +1 / −0 | +6 / −0 | +6 / −0 | +5 / −0 |
| `clean` | gbdt | +1 / −0 | +4 / −0 | +3 / −0 | +2 / −0 |
| `clean_diamond_homology` | logistic | +6 / −0 | +1 / −0 | +1 / −0 | +0 / −0 |
| `clean_diamond_homology` | gbdt | +4 / −0 | +1 / −0 | +1 / −0 | +0 / −0 |
| `clean_diamond_family_excluded` | logistic | +5 / −0 | +1 / −0 | +1 / −0 | +0 / −0 |
| `clean_diamond_family_excluded` | gbdt | +3 / −0 | +1 / −1 | +1 / −0 | +0 / −0 |

**`fixed_0.5`**. Context minus SEQONLY, same source and kind (gained / lost, 37 proteins):

| source | kind | UNIVERSAL | PATHWAY | ORGANISM |
|---|---|---|---|---|
| `diamond_homology` | rule | +0 / −0 | +0 / −0 | +0 / −0 |
| `diamond_homology` | logistic | +1 / −3 | +1 / −4 | +2 / −6 |
| `diamond_homology` | gbdt | +0 / −5 | +0 / −4 | +0 / −5 |
| `diamond_family_excluded` | rule | +0 / −0 | +0 / −0 | +0 / −0 |
| `diamond_family_excluded` | logistic | +0 / −0 | +0 / −0 | +1 / −0 |
| `diamond_family_excluded` | gbdt | +0 / −0 | +0 / −0 | +0 / −0 |
| `clean` | rule | +0 / −0 | +0 / −0 | +0 / −0 |
| `clean` | logistic | +5 / −0 | +5 / −0 | +5 / −1 |
| `clean` | gbdt | +3 / −0 | +2 / −0 | +2 / −1 |
| `clean_diamond_homology` | rule | +0 / −0 | +0 / −0 | +0 / −0 |
| `clean_diamond_homology` | logistic | +0 / −5 | +0 / −5 | +0 / −6 |
| `clean_diamond_homology` | gbdt | +0 / −3 | +0 / −3 | +0 / −4 |
| `clean_diamond_family_excluded` | rule | +0 / −0 | +0 / −0 | +0 / −0 |
| `clean_diamond_family_excluded` | logistic | +0 / −4 | +0 / −4 | +0 / −5 |
| `clean_diamond_family_excluded` | gbdt | +0 / −3 | +0 / −2 | +0 / −3 |

**`fixed_0.5`**. Arm minus native CLEAN at the same cutoff (gained / lost, 37 proteins):

| source | kind | SEQONLY | UNIVERSAL | PATHWAY | ORGANISM |
|---|---|---|---|---|---|
| `diamond_homology` | rule | +0 / −19 | +0 / −19 | +0 / −19 | +0 / −19 |
| `diamond_homology` | logistic | +2 / −15 | +1 / −16 | +1 / −17 | +1 / −18 |
| `diamond_homology` | gbdt | +1 / −15 | +0 / −19 | +0 / −18 | +0 / −19 |
| `diamond_family_excluded` | rule | +0 / −19 | +0 / −19 | +0 / −19 | +0 / −19 |
| `diamond_family_excluded` | logistic | +0 / −19 | +0 / −19 | +0 / −19 | +0 / −18 |
| `diamond_family_excluded` | gbdt | +0 / −19 | +0 / −19 | +0 / −19 | +0 / −19 |
| `clean` | rule | +0 / −1 | +0 / −1 | +0 / −1 | +0 / −1 |
| `clean` | logistic | +0 / −0 | +5 / −0 | +5 / −0 | +5 / −1 |
| `clean` | gbdt | +0 / −0 | +3 / −0 | +2 / −0 | +2 / −1 |
| `clean_diamond_homology` | rule | +0 / −1 | +0 / −1 | +0 / −1 | +0 / −1 |
| `clean_diamond_homology` | logistic | +5 / −0 | +0 / −0 | +0 / −0 | +0 / −1 |
| `clean_diamond_homology` | gbdt | +3 / −0 | +0 / −0 | +0 / −0 | +0 / −1 |
| `clean_diamond_family_excluded` | rule | +0 / −1 | +0 / −1 | +0 / −1 | +0 / −1 |
| `clean_diamond_family_excluded` | logistic | +4 / −0 | +0 / −0 | +0 / −0 | +0 / −1 |
| `clean_diamond_family_excluded` | gbdt | +2 / −0 | +0 / −1 | +0 / −0 | +0 / −1 |

The comparisons give a mixed picture. Learned arms gain more reference
calls than they lose against the rule arm of the same source and context in
most cells, with the clearest gains in the DIAMOND-only sources, but the
CLEAN+DIAMOND sources lose calls under PATHWAY and, at the calibration cutoff,
also elsewhere. Context relative to SEQONLY is not consistent: DIAMOND-only gbdt
gains at the calibration cutoff and loses at 0.5, and the CLEAN+DIAMOND
sources lose reference calls under every learned context at 0.5. The direction of a context
effect therefore depends on the operating point, which is a negative result for
a robust benefit of context on this reference. This agrees with the interim
section, where context gave no AUROC gain in the CLEAN+DIAMOND arms.

Against native CLEAN at the same cutoff, the DIAMOND-only arms lose many
calls; the CLEAN-containing arms differ from native CLEAN by a few calls in
either direction. Against native CLEAN with every call kept, no arm gains
more than a single call and many lose several:

**`calibration_selected`** versus every native CLEAN emitted call (gained / lost):

| source | kind | SEQONLY | UNIVERSAL | PATHWAY | ORGANISM |
|---|---|---|---|---|---|
| `diamond_homology` | rule | +1 / −17 | +1 / −17 | +1 / −17 | +1 / −17 |
| `diamond_homology` | logistic | +1 / −14 | +1 / −14 | +1 / −14 | +0 / −18 |
| `diamond_homology` | gbdt | +0 / −19 | +1 / −15 | +0 / −17 | +1 / −14 |
| `diamond_family_excluded` | rule | +0 / −23 | +0 / −23 | +0 / −23 | +0 / −23 |
| `diamond_family_excluded` | logistic | +0 / −21 | +0 / −22 | +0 / −23 | +0 / −22 |
| `diamond_family_excluded` | gbdt | +0 / −23 | +0 / −19 | +0 / −21 | +0 / −20 |
| `clean` | rule | +0 / −6 | +0 / −6 | +0 / −6 | +0 / −6 |
| `clean` | logistic | +0 / −6 | +0 / −0 | +0 / −3 | +0 / −7 |
| `clean` | gbdt | +0 / −7 | +0 / −5 | +0 / −6 | +0 / −6 |
| `clean_diamond_homology` | rule | +0 / −6 | +0 / −6 | +0 / −6 | +0 / −6 |
| `clean_diamond_homology` | logistic | +0 / −0 | +0 / −5 | +0 / −9 | +0 / −6 |
| `clean_diamond_homology` | gbdt | +0 / −7 | +0 / −6 | +0 / −10 | +0 / −6 |
| `clean_diamond_family_excluded` | rule | +0 / −6 | +0 / −6 | +0 / −6 | +0 / −6 |
| `clean_diamond_family_excluded` | logistic | +0 / −3 | +0 / −6 | +0 / −8 | +0 / −3 |
| `clean_diamond_family_excluded` | gbdt | +0 / −3 | +0 / −6 | +0 / −9 | +0 / −6 |

No arm in the table recovers more reference assertions than native CLEAN keeps
when it emits everything. Where a learned arm has a higher micro-F1 than that
unthresholded baseline, the difference comes from fewer unconfirmed calls,
not from more recovered assertions; with 37 proteins, one or two calls move
the micro-F1 by several hundredths.

**The `clean` logistic arms.** Under the calibration cutoff the UNIVERSAL
logistic arm recovers as many assertions as native CLEAN with every call kept.
This should not be read as a context benefit: the interim section reports that
the Platt slope of these logistic models is near zero, so their output is almost
constant and a cutoff either keeps or drops nearly all calls together. The
recovery is consistent with that degeneracy; this section did not re-inspect
the model files.

**Limits.** The reference is a few dozen proteins from one cohort that overlaps
the training data (exact known sequences, not unseen families); the analysis supplies no
intervals or paired test; unconfirmed calls are not errors; the operating points
are not independent of each other; and nothing here supports promoting a model. The
`conditional_candidate_metrics`, per-genome comparisons and
`reference-protein-scores.jsonl` are in the results directory.

Tables are written by
`benchmark/v3/summarize_source_context_replay.py` from the unchanged
`analysis/results.json` and `replay/verification.json` into
[`benchmark/v3/results/source-context-resume-v2/summary/`](../../benchmark/v3/results/source-context-resume-v2/summary/)
(`source_context_replay_summary.md` with every arm and comparison for all four
operating points, `source_context_replay_results.tsv` with every numeric field,
and `source_context_replay_integrity.tsv`). Two runs produced byte-identical
outputs, with these SHA-256 values:

| file | SHA-256 |
|---|---|
| `source_context_replay_results.tsv` | `7ef08872cd89921e071e80bed425023ceb5966879ccf2d7a4cbefa7c318a7f3b` |
| `source_context_replay_integrity.tsv` | `78ac42648cd2edb6cd72f77d0100346b49c8e4014efc226ad446aa2ad3a6e65c` |
| `source_context_replay_summary.md` | `123229bde3503c4a9c8cf88e9a69b2e5aa349210a93a9448364148f04569f88e` |
