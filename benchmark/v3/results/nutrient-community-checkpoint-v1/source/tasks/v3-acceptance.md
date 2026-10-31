# GSPA: implementation and acceptance of all twelve audit requirements

Objective: fully implement, validate, test, and document **every** requirement in
`tasks/todo.md`, including the detailed item-12 proposal. This is the current
acceptance record. The release claims at `498560f` are not completion evidence.
AllTheBacteria remains excluded as specified in the existing implementation plan.

Status: **in progress**. No whole item is signed off yet. Passing mechanism tests,
running a stage, or producing a report alone does not establish biological validity.
Negative empirical results must remain visible; success cannot be redefined as
an easier simulation or a smaller benchmark after inspecting test labels.

## Acceptance matrix

| Item | Implementation and regression proof required | Empirical acceptance required | Current evidence / remaining work |
|---|---|---|---|
| 1. Actual application | Production inference for every arm; content-addressed execution evidence; all failures retained; no label-dependent input generation. | Whole-proteome claims; matched stage/full-system comparisons; frozen family/genome splits; calibration-only thresholds; calibration, correction/harm counts and paired genome uncertainty. | Verified-batch receipts, invalidation tests and process-group timeout cleanup pass. All 132 full-matrix parts have verified publication receipts (1,526,919,081 claims); collection-wide job 10109 independently rehashed all 851 bound files successfully. TopEC/EnzBERT and complete DeepProZyme baseline/default-prior capacity pass. Indexed dense job 10080 and independent comparison 10108 completed with exit 0: all 9,790,188 identities and full prediction bytes match, with 37 files rehashed twice; larger/cohort capacity, complete comparisons and biological validation remain pending. |
| 2. Routes and state | Single CLI, annotation, web, batch and Nextflow share scores, provenance and candidate identity, including joint cohort execution. Requested/missing/no-work/failed states tested. | Matched-input comparisons through each actual route, including multi-genome context and mixed suggesters. | Shared workflow and no-work RLC state repaired. Nextflow now runs one joint batch; 18 real genomes match direct CLI predictions, full cohort reports, learned candidates and stage states. Eight live Nextflow regressions cover parity, cache/input bindings and failures. The typed ClaimKey parser now preserves pipe-containing protein IDs; old/new CLI replay proves the targeted-rescue crash fixed with ordinary-ID scores unchanged. The AntiFam quality path now preserves typed observations and scan coverage through annotation/direct integration, with unchanged heuristic functional scores and explicit learned abstention. Other upstream annotation, route combinations and reporting remain pending. |
| 3. Predictors | eggNOG columns/scores, HMMER database, genome dispatch, InterPro output types, BP-only operon transfer; typed failures and partial output. | Versioned genuine tool outputs and complete parser/dispatch paths, including injected failures and coverage loss. | Real HMMER 3.4 exposed repaired availability, malformed-statistic/accession and CLI-thread defects. Versioned synthetic-input observed outputs, three live adapter tests and old/new annotation CLI replay now pass; annotation/direct integration predictions agree. Older mock files are explicitly format fixtures. The separate AntiFam adapter now passes equivalent availability/contract repairs and real synthetic-model execution, and its quality observations now survive raw export, integration and annotation projection; current learned models explicitly abstain on flagged proteins. See the AntiFam quality repair record for scope and verification. Observed InterProScan 5.72 output now reproduces and verifies repaired statistic types, malformed-value rejection and thread forwarding; the actual annotation/direct-integration replay passes with nine identical raw claims and byte-identical predictions. A genuine gapsmith training-genome search now exposes and verifies repaired input, configuration, parser, candidate identity and reporting contracts. All 13,147 observed rows survive seven production CLI replays, producing 3,469 candidates; annotation/direct predictions agree and both learned models abstain on the new source. Other adapters, real-family accuracy, nucleotide attribution and full dispatch/coverage-loss paths remain to be certified. |
| 4. Learned inference | Revisable calibrated scoring, correlated channels, trained no-candidate outcome; no irreversible floor or candidate-count artifact. | Honest complete locus labels; OOF training; family/genus holdouts; logistic/tree/rule comparisons; calibration and harm; supported predictor distributions. | The legacy 100-candidate scoring discontinuity, negative polarity, per-genome locus grouping and post-filter completeness checks are repaired. Both models now use contract v3, excluding raw cohort counts, and were refitted with all 4,262 supervision rows preserved; their target is reference support. A separately augmented table now has one reviewed experimental negative and one complete two-member training competition; no calibration/test or all-negative groups. Null remains untrained. Actual source identity and frozen train/calibration profile checks now prevent unsupported rescoring, with matching JVM/Python and real CLI verification. Predictor-version/score-distribution checks, broader null/role data and independent evaluation remain pending; prior EC-specialist transfer failures remain unresolved. |
| 5. Reaction/interaction context | Explicit shared candidates, alternative EC mappings, direction, organism feasibility, independent interaction roles. | Pathway-only versus universal-network versus organism-specific ablations on held-out enzyme loci with actual gaps/context. | Missing GO evidence incorrectly contributed probability 0.5 to RLC diversity, commitment and self penalties. Three failing-before/passing-after tests and ten production CLI runs now verify the repair. Functional TSVs remain unchanged because RLC emits separate suggestions. Prior RLC runs had no gaps. A verified training-genome gapsmith observation now supplies 50 candidate-gap records across 46 reactions, without held-out labels. The single-genome local-context bridge now preserves local/intrinsic features with recurrence unavailable: eight real CLI runs on 29,079 unchanged claims and the 50 observed gap hypotheses verify 1,182 local contexts, unchanged sequence evidence and exact peer-run feature compatibility. Both models abstain on these new local-context inputs. Observed proteome-search gaps now reach shared integration automatically: ten real CLI replays verify all 50 audited records, file/in-memory/explicit-gap feature parity, complete/zero/missing/failed distinctions, and unchanged raw claims. Standalone RLC anchor mapping and alternative preservation are now repaired and verified. Full ablations, calibrated contribution and broader input coverage remain pending. |
| 6. Joint recurrence | Two-pass cohort inference, singleton conservation, density-matched background, ancestry/completeness weighting, complementary gaps. | Clade holdouts; duplicate strains; shuffled neighborhoods; genuinely fragmented/missing sequence; candidate-level and genome-level effects. | Real nucleotide-loss controls remove sequences, claims and memberships together (54 verified genome/arm runs). Contract-v3 replay previously passed exact-copy probability invariance on 118 outputs but exposed a heterogeneous-cluster defect. Complete observed-input identities now repair that defect, with 27 added JVM cases and a successful replay of the original failure. Both models are refitted from 484,202 re-exported candidates with unchanged supervision; updated controls independently verify all 118 outputs with exactly invariant numeric probabilities for both models across all 18 original genomes. Independent validation, near-duplicate/fragmentation scope and complementary-gap acceptance remain open. |
| 7. Retention and roles | Weak-candidate retention independent of display; catalytic/accessory/transporter/regulator/paralog/subunit distinctions. | Real predictor scores and difficult same-operon alternatives; experimentally supported role truth; recall and locus ranking, including no-candidate cases. | Gapsmith now retains weak hits and preserves complex-subunit candidates separately from whole EC activity, including all source rows. The raw-claim role handoff is now repaired: 829 observed complex-subunit roles survive both learned models, with 14 mechanism tests and before/after production proof; both models still abstain on this source. Independent role training remains incomplete. The reviewed ThrRS catalytic/editing pair provides one complete GO-function training competition, not a role-labelled benchmark. Raw roles now also survive masked experiment construction; complex candidates require supported partners, and mixed catalytic/subunit absence candidates cannot bypass that requirement. Full context and role truth/training remain incomplete. |
| 8. Absence/environment | Competing missing-annotation, missing-sequence, pathway-loss and external-supply explanations; explicit medium/biomass/transport assumptions. | Gene attribution scored separately from growth/KO accuracy on complete genomes, auxotrophs and controlled sequence loss. | Explicit four-scenario absence comparison now runs through the production workflow with retained-candidate binding, bounded supply/transport, unavailable states and whole-condition budgets. Nine CLI runs and 147 independent GLPK comparisons pass; no gene attribution is promoted. A newly audited two-species nutrient dataset supplies 2,464 real assay inputs (2,240 primary), with raw-derived uncertainty and outcome-independent case construction. All 1,144 single-organism assays now run under six frozen arms (6,864 conditions), with one solver-unknown prediction retained. Paired rescue accuracy remains below the 80% constant-rescue comparator in both species for every arm; numerical verification and its solver-history diagnosis are documented in the nutrient-phenotype report. Biological absence discrimination, actual sequence-loss evaluation, empirical cohort-wide absence comparisons and separate gene-attribution accuracy remain pending. Automatic construction from run gaps now has shared budgets and preserves unavailable evidence, verified by 100 synthetic CLI/GLPK comparisons. |
| 9. Cross-feeding | Separate organism models, transport and extracellular exchange under common environmental constraints; conservative attribution. | Measured mono/coculture, exchanged metabolites and transporter KOs; separate false enzyme assignment effects. | Growth-rescue benchmark exists. New primary raw data include 1,200 primary-condition coculture cases alongside matched monocultures and supplementation; organism identities stay separate. Model evaluation on these data, exchanged-metabolite/transporter validation and annotation effects remain pending. |
| 10. Reporting | Stage-wise baseline/final changes, evidence, alternatives and uncertainty; ontology/EC enzyme classes; browser/export reconciliation. | Real production output rendered with correct no-work, failed, missing, learned/null/role and partial-data states. | Portal mixed-scoring and complete-abstention reports now agree with the results browser; nested stage assumptions and unknown values are readable. Seven actual stored production-output replays pass through the packaged CLI, FastAPI endpoints and Chromium, preserving all candidate identities/scores and complete input-protein denominators. Gates pass 46 service and 73 visualization tests. The rebuilt JAR differs only in the visualization resource. Broader route, provenance and production-state coverage remain pending. |
| 11. Agent comparison | Public evidence separated from truth before selection; opaque IDs; no invented scores; identical inputs/guardrails/budgets for simple and learned methods; unresolved-only calls. | Unseen genomes/families; real unresolved cases; calibration-frozen decisions; correctness/harm/abstention and complete inference/tool cost. | Invalid v1 generator replaced; 21 fresh production runs yield 196 real cases, but none has matching experimental truth. All seven calibrated methods therefore abstain. Live agent completed 392 calls with usage recorded; independent outcome data and full validation remain required. |
| 12. Experiments | Evidence-informed competing assignments; bounded carbon/nitrogen/oxygen/nutrient/KO/product experiments; alternative feasible fluxes; negative/inconclusive cases and all assumptions retained. | Mask target assignments from ALL evidence/model inputs; held-out organisms and measured media/KO/secondary-product labels; single-medium, fixed, random, simple evidence and agent under equal solver budgets; fixed-precision recovery/harm/abstention/cost. | Masked GPR scenarios run through Python/Java with null, paired phenotypes, FVA, products and common budgets. Guided construction now uses raw/context candidates and independent reaction mappings; five production CLI mechanism runs agree with GLPK on 54 comparisons and preserve hidden-assignment invariance. Earlier real-model runs still disagree with the measured phenotype. The new nutrient/supplementation benchmark also fails to beat a constant-rescue comparator, despite improved individual-assay agreement in the unprinted-inorganic sensitivity arms. Supported complex/multifunction clauses now pass 200 independent synthetic CLI comparisons; combined knockouts and oracle isolation are repaired. Real partner/architecture evidence, calibrated recovery and held-out biological/product validation remain required. |

## Work order and verification

1. Repair execution/evidence integrity and truthful reporting; preserve old artifacts
   as historical, with their limitations attached.
2. Verify/acquire complete prediction inputs and traceable independent labels for
   loci, roles, phenotypes, transport and products; freeze splits and exclusions.
3. Complete route parity, context/absence/role integration, train calibrated
   marginal and null models on supported channels, and run all production ablations.
4. Replace the agent evaluation and implement masked-function experiment recovery;
   run frozen equal-budget comparisons on held-out experimental truth.
5. Run regression, route, browser and empirical gates; inspect source/output hashes;
   reconcile documentation and perform a requirement-by-requirement completion audit.

Environment: use `TMPDIR=${TMPDIR:-$HOME/tmp}` for scratch work, with Java
`-Djava.io.tmpdir` set to that directory. Python test environment:
`$HOME/tmp/gspa-v3-testenv/bin/python`; set `GSPA_TEST_PYTHON` for JVM sidecars.
Do not infer a live process from old task notes; inspect the actual job/session.

## Verified progress

- 2026-10-27: [Predictor role handoff](../docs/benchmarks/v3-claim-role-bridge.md) now preserves supplied roles through the shared learned extractor and raw JSONL route. The preceding JAR fails the real-output role oracle; the repaired JAR passes all seven CLI replays. Exactly 829 complex-subunit identities change, while the other 2,640 feature digests and all functional TSVs remain unchanged. Fourteen added regressions, 1,161 JVM tests (two skips), and nine live Nextflow cases pass. This is evidence preservation, not trained role accuracy or whole-item acceptance.

- 2026-10-26: `IntegrationRunAudit` distinguishes a requested stage with no eligible
  work from executed inference. RLC with no gaps records `skipped`, reason and
  unchanged candidate counts; missing graph/layout still fails. Regression tests
  cover state preservation and the actual shared workflow.
- 2026-10-26: `verified_batch.py` replaces file-count reuse in cohort ablations and
  perturbations. Receipts cover command, executable/JAR, inputs/references/models,
  runner source, stage ledgers and output hashes. Failed/skipped required stages
  cannot certify an ablation. Old outputs are archived; changed inputs invalidate
  reuse. Thirteen process-level regression tests pass (stub engine, not biological evidence),
  including `--option=value`, multiple module files and directory dependency changes.
- 2026-10-26: the main `v3_runner.py` now checks a pre/post inference identity and
  verifies it on evaluation-only reuse. Changed launch scripts, references/models,
  auxiliary files, arguments and genome context invalidate reuse; legacy records
  without an identity are refused. Gradle build sources/resources are hashed.
  Schema-2 ledgers must identify the genome and execute requested ablation stages;
  changed sidecars also invalidate reuse. Six new hardening tests reproduced the
  failures before the fix. All 33 hardening tests pass after it. Indirect wrapper/
  manifest dependencies must be declared in `gspa.input_files`; package versions
  and remote service state remain outside the content receipt.
- 2026-10-26: JVM verification: 948 tests, zero failures/errors, two skipped;
  `:gspa-cli:shadowJar` succeeds. Python benchmark verification after the main
  runner changes: 446 tests, three optional checks skipped. The real-engine
  runner suite was also run with `GSPA_V3_ENGINE_SMOKE=1`: 22 tests pass,
  including actual Gradle-launched integration. The subsequent two additional
  verified-batch cases pass in the 13-test focused suite. These are engineering
  checks, not acceptance of the outstanding biological claims.
- Full Price input scale from the committed matrix audit: 1,517,016,382 nonzero
  complete-EC scores across six predictors. At an assumed 300 bytes per JSONL
  claim this would be about 455 GB before ensembles. This estimate is not an
  exported dataset or a production run. Whole-proteome inference still needs an
  explicit, validated memory/storage strategy; silently reusing truncated inputs
  is not an acceptable substitute.

The old Price context report used the same 154 claim-file hashes as the earlier
truncated-input evaluation. Both shipped models have `locus_model.trained=false`.
The old Stage-3 evaluation embeds answers and uses hard-coded model probabilities.
These artifacts must not be recertified without rebuilding the relevant inputs,
running the corrected pipeline and obtaining the required independent evidence.

## Replacement agent comparison (2026-10-26)

The v1 truth-selected generator has been removed from the executable evaluator.
The v2 pipeline separates verified production export, public prediction,
calibration-only policy selection and test scoring. Missing labels stay unknown,
constraints are identical across methods, unsupported thresholds abstain, failures
remain failures, and model scores come from actual matched production runs.
Sixteen comparator/evaluator regressions and the private positive-join regression
pass. The full benchmark suite after the evaluator replacement ran 462 tests
with three optional skips and no failures; the private-join test was added and
passed afterwards.

Fresh inference completed 21 genome/arm runs, yielding 38 calibration and 158 test
anomalies. The experimental truth join supports **zero** selected assignments.
Thus the run validates execution and exposes a data gap; it supplies no evidence
of accuracy. Every simple method's frozen policy abstains. The live agent also
completed all 196 cases: 392 calls, no failures, 589.932 seconds, 730,626 prompt
tokens and 71,169 completion tokens. Its calibrated policy likewise abstains.
Source receipts, public cases, private truth provenance, proposals, policies,
reports and checksums are archived in
`benchmark/v3/results/stage3-agent-evaluation-v2/`. See
[the version-2 report](../docs/benchmarks/v3-stage3-production-comparison.md).
Item 11 remains open for independent biological validation and the wider
requirements in the matrix above; the objective is unchanged.

## Masked metabolic assignment mechanism (2026-10-26)

The production experiment entry point now supports explicit competing assignments
after deleting the target's original GPR clauses and aliases. OR alternatives and
AND partners are preserved. Planning sees only opaque hypotheses, bounded numeric
evidence, structured conditions and previous conditional outcomes. A null
assignment is required. Complete candidate comparisons share one solver ledger;
growth, infeasibility, errors, alternative fluxes and conditional products remain
distinct. Java verifies requested scenarios/budgets and reports them without
promotions or calibrated attribution.

Eighteen focused Python tests, the existing 53 metabolic tests, the real Python
CLI and shared Java workflow pass. The full benchmark suite passes (463 tests,
three optional skips); the full Java suite passes (949 tests, two skips). The
Python script suite's Gradle-calling visualization test conflicts with concurrent
Gradle `test`; the Java suite was rerun separately after identifying that test
infrastructure interaction. The production JAR was rebuilt successfully.

On the stored iML1515 model, a documented hand-specified mechanism example ran
three assignments under two conditions: 20 solver calls for both fixed and live
Qwen selectors. COBRApy/GLPK agreed with all 12 control/knockout states and their
flux intervals for each arm. Qwen made two successful API requests (7,601 input
tokens, 16 output tokens). All assignments predict alanine growth despite the
registered experimental no-growth label: that mismatch is retained, and no
function recovery is claimed. See
[the mechanism report](../docs/benchmarks/v3-masked-assignment-experiments.md) and
`benchmark/v3/results/masked-assignment-mechanism/`.

Automatic generation from real candidate evidence, independent hidden function
truth, calibration and the full held-out biological/product experiment benchmark
remain outstanding. The explicitly supplied-manifest path is implemented; the
`evidence_guided` builder rejects this new format until those mappings are
implemented. This work advances item 12 without accepting it as complete.

## Actual sequence-loss control (2026-10-26)

`prepare_sequence_loss.py` deletes deterministic 10 kb blocks from real assemblies,
keeps only complete CDSs in retained fragments, and filters protein FASTAs, claims
and orthogroups consistently. Unknown-position proteins stay explicitly outside
the intervention. It reads no labels. The matched intact/coordinate-loss/sequence-
loss runner completed 54 production genome/arm runs with content-verified receipts.
The independent input audit verified exact nucleotide slices, retained protein
translations, byte-identical claims and filtered memberships. Eleven focused tests
and the full benchmark suite pass (474 tests, three optional skips).

The evaluator preserves genome/protein identity instead of merging identical
sequences. Removed proteins remain in the original reference denominator and are
also reported separately from retained proteins. Calibration is frozen before test
scoring, GO ancestor dependencies and incomplete positive-only truth are explicit,
and unmatched functions remain unverified rather than experimental negatives.

EC top-1 reference matches are 14/37 intact, 12/37 under coordinate loss and 6/37
under actual sequence loss. Among 13 retained labelled proteins, the corresponding
counts are 7, 6 and 6. No improvement is demonstrated. The EC policy meets its
90% reference-support rule on calibration (25/27 calls) but not intact test (6/12);
GO has no qualifying threshold and abstains. The policy is not retuned to this
negative result. These data advance items 6 and 8 and expose a calibration issue
under item 4; they do not complete any of those items. See
[the protocol and report](../docs/benchmarks/v3-sequence-loss.md) and the archived
`benchmark/v3/results/sequence-loss-v1/` artifacts.

## Nextflow joint cohort route (2026-10-26)

`gspa-nf/modules/cohort.nf` stages all genomes' claims/layouts and shared references,
then invokes the production two-pass batch workflow once. Both `main.nf` after
claim merging and the integration-only `cohort.nf` entry point use it. Per-genome
provenance, cohort/learned reports, stage ledgers and optional feature exports are
preserved. Required/malformed inputs fail explicitly; runtime cohort failures leave
failed ledgers for every genome in the Nextflow work directory.

Live tests reproduced two orchestration bugs before fixes: swapping claim files
between genomes could reuse an old cache entry, and an empty manifest could exit
successfully without inference. Explicit file bindings now enter the task hash,
deep hashing catches same-size/same-timestamp content edits, and empty collections
reach the required-input validator. The official Nextflow launcher is necessary;
direct framework-JAR invocation omitted JVM options needed for cache persistence.

The final Nextflow 26.04.6 execution matches the verified direct batch CLI over all
18 development genomes: 484,202 prediction rows byte-for-byte, all cohort reports,
all learned candidate records/model identity, and stage states. Inputs/engine hashes
and per-genome staged bindings verify. Eight live Nextflow regressions pass; the
full benchmark suite passes (482 tests, three optional skips) and 15 service config
tests pass. The main annotation entry point is preview-checked with the existing
Groovy parser, not a full upstream predictor run. See
[the route protocol](../docs/benchmarks/v3-nextflow-cohort.md) and
`benchmark/v3/results/nextflow-cohort-v1/`. Item 2 remains open for its other route
and reporting requirements; no biological accuracy claim follows from parity.

## Evidence-guided masked assignments (2026-10-26)

The production experiment workflow now constructs masked assignment hypotheses
from the requested protein's raw EC claims and measured neighborhood/cohort
context. It retains weak candidates, records the raw sources/scores, validates
independent reaction-level mappings, separates alternative reactions, and includes
null. Original target GPRs, final probabilities, metabolic/learned inferred claims
and refiner reaction tables cannot select the hypotheses. Conditions come from
exact chemical identity or explicit supported exchange hints; all inputs still
pass the sidecar's structured-condition and boundary checks. Explicit non-catalytic
and complex-only roles are not silently made catalytic. Automatic complex and
multifunction clauses remain open; supported explicit clauses remain available.

Fourteen added JVM cases cover construction and the real Java/Python path. The
full JVM suite passes: 963 tests, zero failures/errors, two skips; the JAR builds.
A cleanup regression also reproduces and fixes an old construction report
surviving a later run with experiments disabled. A concurrent Gradle recompile
invalidated one intermediate test run; that evidence was discarded, and the final
full suite was executed serially against the final sources.
The final benchmark suite passes (482 tests, three optional skips, including eight
live Nextflow checks against the rebuilt engine); all 18 Python masking tests pass.

Five built-CLI executions on a declared synthetic network cover single-medium,
fixed, seeded-random and evidence-guided selection, plus a changed hidden original
assignment. All use the same generated public problem and the same 20-call cap;
the single-medium arm uses 10 calls and the others 20. Changing the original target
assignment leaves the public problem and comparison cells unchanged. Independent
COBRApy/GLPK checks agree on all 54 control/knockout states and their checked flux
intervals. These are mechanism checks, not an accuracy benchmark or an agent
comparison. The earlier real-model negative phenotype result remains unchanged.

See [the updated protocol](../docs/benchmarks/v3-masked-assignment-experiments.md)
and `benchmark/v3/results/guided-masked-assignment-v1/`. This advances item 12;
source-level masking/splits for real predictor data, complex/multifunction
construction, frozen confidence policies and held-out biological/product
validation remain required. No whole item is accepted by this work.

## Locus/null truth integrity and real-data replay (2026-10-26)

Five failing preparation regressions and a failing trainer regression reproduced
three problems: an explicit negative label was read as positive; incomplete groups
could survive label filtering; and equal locus names in different genomes merged.
The repaired builder preserves explicit polarity and referenced negative truth,
keeps GO negatives exact, rejects conflicting assertions, and verifies each group
against the complete exported function/role/genome candidate universe after
label joining and family masking. Incomplete groups retain marginal rows but
cannot teach null. Production training/evaluation verify the membership/source
certificate again. Unconfirmed closed-world labels cannot supervise null.

Ten preparation and six trainer regressions were added. The full benchmark suite
passes (492 tests, three optional skips, including eight live Nextflow tests);
all 40 trainer tests pass. No JVM implementation changed. Rebuilding the real
18-genome cohort preserves every original feature, label and split for all 4,262
rows. Both models refit with exactly unchanged scoring/calibration parameters;
the new artifacts explicitly identify their closed-world reference-support target.
Four production CLI runs verify old/new logistic and tree artifacts on one real
genome: each pair has 26,328 byte-identical prediction rows and identical learned
candidate records.

The cohort has 541 positive candidate rows, 3,721 unconfirmed rows, no explicit
negative rows and no complete locus declarations. Both locus layers therefore
remain correctly untrained. This fixes prerequisites for items 4/7, not the
missing empirical validation. See [the protocol](../docs/benchmarks/v3-locus-truth.md)
and `benchmark/v3/results/locus-truth-v2/`. An available next data source to audit
is the local 11 GB `benchmark/neural/cafa6_recon/goa_uniprot_all.gaf.gz` snapshot;
its negative assertions and suitability for complete groups have not been checked.

## GOA negative audit and first reviewed competition (2026-10-26)

The complete GAF scan now verifies all 667,198,446 lines, gzip integrity and
pre/post source/code/mapping hashes. It finds 10,526 negative annotations,
4,117 with experimental codes. Thirteen rows join cohort accessions; one is
nonexperimental. Exact intersection with all 484,202 production candidates in
18 genomes finds four candidate matches. Relations, isoforms and assay context
are preserved; neither missing annotations nor accession matches create truth.

Full-text review of PMID:15240874 and exact sequence/versioned-accession checks
support one complete two-protein GO:0004829 competition in the existing training
genome GCA_000007005.1. The standard builder and trainer validator accept its
complete membership after family masking. A separate augmented table preserves
all 4,262 original identities/features/binary labels/splits; one unconfirmed zero
now has explicit negative support. There are no complete calibration/test groups
and no all-negative outcomes, so no null model was fitted or promoted.

Other assertions remain excluded or pending, including a source/reaction mismatch
for Q6LZY8 and broader process negatives that cannot be copied into generic GO
labels. The audit does not claim an exhaustive review of every reference.
See [the source audit and reproduction protocol](../docs/benchmarks/v3-gaf-negative-audit.md),
`benchmark/v3/results/gaf-negative-audit-v1/` and
`benchmark/v3/results/reviewed-thrs-locus-v1/`. Items 4 and 7 remain open.

The final benchmark suite passes: 504 tests, no failures, three optional skips,
including eight live Nextflow tests. Twelve new tests cover inventory/coverage
integrity and exact sequence/version checks. No JVM implementation changed.

## Learned-model input coverage (2026-10-26)

Actual source names now survive extraction, labelled-table preparation and the
training/calibration digest. The trainer freezes exact source/function/role/
availability profiles using only train/calibration data. Both classes must occur
in both splits. Unsupported production candidates receive null learned scores,
retain their preceding annotations/floors and cannot form a filtered locus
competition. Reports, exports and the browser preserve the reasons. Partial
coverage records mixed scoring even when all supported scores remain report-only;
an entirely unsupported workload records a skipped learned stage.

The real 18-genome replay preserves 484,202 numerical vectors and 4,262 training
rows, with unchanged fitted parameters for both models. Only 462/608 previously
examined test rows meet the coverage rule; 146 remain explicit exclusions. Null
is still untrained: one complete positive/negative training group does not supply
an all-negative outcome or calibration/test groups.

Six actual CLI executions compare both models to matched baselines. For the real
homology genome, 18,698 scores agree with Python to below 4e-15, and 7,630
unsupported function rows remain unchanged. All 94,338 compact CLEAN candidates
remain unsupported; both model runs produce TSVs byte-identical to baseline.
These truncated, previously examined Price inputs are not a whole-proteome or
independent accuracy benchmark. Real mixed output also renders in Chromium with
26,328 joined candidates, correct unsupported reasons and no JavaScript errors.

The final JVM suite passes 968 tests (two skips, zero failures/errors), including
a regression that first reproduced the report-only mixed-state defect. The
Python model tests pass 47 cases; browser payload tests pass 37 and real Chromium
tests pass 14. The benchmark suite passes 504 tests with three optional skips;
all eight live Nextflow regressions also pass after the final build.
See [the protocol and limitations](../docs/benchmarks/v3-model-input-support.md)
and `benchmark/v3/results/model-input-support-v1/` for receipts and checksums.

This supplies a necessary categorical coverage guard. Predictor versions,
reference databases, score distributions and taxonomic transfer remain unchecked;
profile coverage is not proof of calibration/generalization. Items 2, 4, 7 and 10
remain open, as do all other whole-item acceptance gates.

## Lossless Price matrix export and capacity measurement (2026-10-26)

The historical Price runner still points at the earlier compact input directory;
its full-context completion claim has now also been corrected in the data/export
protocols. The 132 verified source matrices contain 1,526,919,081 nonzero complete
and preliminary EC scores. A frozen all-matrix plan exports on node-local scratch
and publishes only independently verified, hash-matching copies to persistent
storage. Raw scores/identifiers and weak/unplaced/unresolved candidates survive.
No threshold or top-k is introduced to fit the benchmark into available resources.

The first real run exposed a false target-statistics check: the compact extractor
had converted to float32 and remapped obsolete ECs using maximum reduction, while
the full exporter compared untransformed values. Array 9890 was stopped and its
failures retained. The corrected comparison reproduces only that historical
projection and now checks complete target-vector hashes as well as counts/maxima.
Its mapping reference matches the archived extraction snapshot exactly. Raw export
values and identifiers remain unchanged. Previously failing TopEC and EnzBERT
parts now pass every-claim verification (7,163 and 349,400 claims respectively).

The corrected full plan is `/data/hohndor/gspa-v3/price-full-matrix-v2/`: array 9965
completed the two regression cases; array 9967 processes the other 130 matrices,
two at a time. This remains an ongoing acquisition, not a completed benchmark.
Inspect live Slurm state and receipts before continuing; do not infer liveness
from this note. The persistent store had 7.2 TB available at launch.

Production capacity job 9913 completed stage 0 and stage 1 for the complete TopEC
GCF_000236665.1 matrix. Both preserve all 3,481 nonzero candidates from 3,494 matrix
rows. Peak RSS was 4.68/4.77 GiB; elapsed time was 44.77/1,010.67 seconds. Stage 1
rebuilt the ontology dependency cache, dominating its first-run time. This small,
sparse matrix does not prove dense-input or cohort capacity. Job 10011 measures
the already verified 349,400-claim EnzBERT matrix separately.

Eleven new regressions cover independent matrix verification, corrupted or missing
claims, wrong accession versions, changed inputs, immutable publication, damaged
copies and compact-statistics projection. The focused suite passes 37 tests; the
full benchmark suite passes 515 with three optional skips, including live Nextflow.
No JVM source changed. See [the protocol](../docs/benchmarks/v3-price-streaming-export.md)
and `benchmark/v3/results/price-full-matrix-v1/` for the partial evidence snapshot.

Remaining under item 1: finish/reconcile all exports, measure and support dense
inference/cohort memory requirements, consistently handle obsolete/current EC
identity in biological evaluation, and run the complete frozen comparison.
Original predictor FASTA identity, independent calibration/truth and every other
whole-item acceptance gate remain open. These engineering results do not change
that status.

## Larger capacity check and bounded candidate comparison (2026-10-26)

Job 10011 completed both production capacity arms for the full 349,400-candidate
EnzBERT GCF_000236665.1 matrix. The no-priors combiner takes 94.47 seconds and
7,782,496 KiB peak Java RSS; default priors take 403.07 seconds and 6,492,324 KiB.
Both outputs retain every candidate and are byte-identical. This is retention and
resource evidence, not an accuracy improvement. These capacity arm names must not
be confused with the raw-predictor maximum-score accuracy baseline.

The earlier checker materialized the expected key set and entire output in Python.
It now externally sorts identity keys on job scratch with a 64 MiB sort buffer and
compares every key, including duplicate detection. Eight regressions cover exact
identity, missing/extra/changed/duplicate candidates, malformed inputs, bounded
Python allocation, and retained failure evidence. Independent job 10040 compares
both real EnzBERT outputs with the new checker and verifies 33 recorded file hashes.
The broader benchmark run reports 523 tests with no failures and 11 skips: the
eight live Nextflow opt-ins plus three other optional checks. Existing live results
apply to the unchanged JAR/route. No JVM source changed.

Job 10037 runs complete DeepProZyme GCF_000236665.1 baseline integration on
9,790,188 candidates with a 64 GiB heap/80 GiB allocation. At the latest check it
has loaded all claims and entered inference, but has not completed. The frozen
code/launch metadata and EnzBERT receipts are archived under
`benchmark/v3/results/price-capacity-v2/`. Its export progress snapshot records
70/132 parts and 55,334,600 claims at 17:44 UTC, with no failures. This snapshot
counts worker receipts; it does not rehash all published bulk bytes.

Dense/cohort memory, complete export reconciliation, full ablations and independent
biological validation remain open. The disk-based candidate comparison does not
make the JVM or the separate full-JSON stage-ledger validator memory-bounded.

## Remove the legacy 100-candidate scoring switch (2026-10-26)

The legacy refiner still gated its weak-evidence positive-prior ceiling on the
total number of protein/function candidates. A regression reproduced a score
change from 0.92765 to 0.52410 when an unrelated candidate crossed that boundary.
The candidate's own likelihood and prior contribution now determine the ceiling;
the rule already used for runs of at least 100 candidates remains unchanged.
Strong positive and negative context contributions retain their prior behavior.

Eight actual production CLI genome/build runs independently reproduce the defect
using a declared synthetic metabolic gap: the old target score changes from
0.96095 at 99 candidates to 0.52434 at 100. The repaired score is 0.52434 for
1/99/100/101 candidates, and the complete 100/101-candidate outputs are
byte-identical between builds. Content receipts bind both JARs and every input
and output. These are mechanism fixtures, not biological labels or evidence
of improved accuracy.

The complete JVM suite passes 972 tests (two skips, no failures/errors); the JAR
build succeeds. The benchmark suite passes 523 tests with three optional skips,
including eight live Nextflow checks against the rebuilt JAR. The learned-context
overview no longer incorrectly says that real training has never occurred, and
documents the version-2 contract, supported-input guard and untrained null layer.
See [the repair and reproduction protocol](../docs/benchmarks/v3-candidate-count.md)
and `benchmark/v3/results/candidate-count-v1/`. The preserved repaired JAR is
`58dbebf1cf6ad97cd216d676b0cf1221600a97c7350f6f45abf27b54c9669654`.

The cap itself remains a legacy heuristic. Null/role training, full biological
calibration and generalization remain outstanding; item 4 is not complete.
Running remote capacity jobs retain their frozen older engine and are not
silently recertified against this build.

## First dense baseline capacity gate completed (2026-10-26)

Job 10037 has completed with its original frozen `c469273...` JAR. All 9,790,188
DeepProZyme GCF_000236665.1 candidates survive production no-priors integration
and the final exact disk-based key comparison. JVM elapsed time is 930.31 seconds,
with 56,327,604 KiB peak RSS (53.72 GiB). Output consists of an approximately
859 MiB prediction TSV and 2.5 GiB full stage ledger. Job 10053 independently
rehashes 25 recorded input/output files after measurement and confirms every hash.
Receipts, resource records, code identity and the recheck are archived under
`benchmark/v3/results/price-capacity-v3/`; bulk outputs stay at the recorded paths.

Array 9967's throttle is now four, following a live resource check (40 logical
CPUs, 128,722 MiB configured node memory; each export requests two CPUs/8 GiB).
This scheduling change does not alter the frozen export code, sources or inputs.
The all-matrix export is still incomplete. Baseline capacity for this single dense
matrix does not establish larger-matrix, default-prior or joint-cohort capacity.
IBEX currently offers batch nodes with at least 336,896 MiB and large-memory nodes
with at least 3,018,752 MiB; these are observed resources, not inference results.

## Observed HMMER and timeout regressions (2026-10-26)

[HMMER verification](../docs/benchmarks/v3-hmmer-observed.md) records real
HMMER execution on synthetic inputs, failed before-fix cases, and old/new CLI
outputs. A valid installed tool previously produced UNAVAILABLE because
`hmmsearch --version` is unsupported. Corrected probing, strict numeric/accession
contracts and CLI thread forwarding now pass. Seventeen adapter cases include
three live-tool checks; two added CLI cases and the pipeline assertion cover the
thread setting. Annotation and direct integration agree byte-for-byte under the
same no-prior workflow; this is not a biological accuracy result. Full JVM tests:
991, zero failures/errors, two skips. Full benchmark tests with the final cached
JAR: 525, zero failures, three optional skips; all eight live Nextflow cases run.

The timeout regression reproduced a surviving SIGTERM-ignoring descendant after
the old runner timed out. The repaired runner terminates the separate process
group and writes no success receipt; it cannot guarantee cleanup of a process
that deliberately creates a new session. All 23 focused runner/capacity tests
pass. `benchmark/v3/results/batch-timeout-v1/` retains the old code, failed log,
passing gates and frozen dense-prior launch. Capacity results now record both
per-arm inference deadlines and total elapsed time including verification.
Slurm 10058 uses a four-hour inference deadline for the full 9.79-million-claim
DeepProZyme stage1 run, on the original capacity JAR. Results remain pending;
its running state does not establish dense-prior capacity.

AntiFam follow-up: the same unsupported version probe and malformed-statistic
contracts were reproduced in the second hmmsearch adapter. Twelve new cases
failed before repair and pass afterwards, including real matching/no-hit and
missing-gathering-threshold failure paths. Both adapters share the verified
probe and honor CLI thread settings. The final JVM suite passes 1,003 tests
(two skips), and all eight live Nextflow tests pass against that exact JAR.
Old/new AntiFam CLI runs change from UNAVAILABLE/zero annotations to
SUCCESS/one expected synthetic PSEUDOGENE annotation. Both emit zero raw
integration claims; the source-to-evidence map lacks AntiFam. Functional
rescoring from that quality flag is still unverified. The four-run direct-HMMER
route comparison passes again after extracting the shared probe. Evidence:
`benchmark/v3/results/antifam-observed-v1/`; no biological accuracy claim.

## Essentiality candidate indexing and protein identity (2026-10-26)

[Verified investigation](../docs/benchmarks/v3-essentiality-index.md) found a
second key parser splitting UniProt-shaped protein IDs at the wrong delimiters.
Three new round-trip tests and actual old/new CLI inference reproduce and repair
the targeted-rescue crash; ordinary-ID scores remain unchanged.

A live stack from dense-prior job 10058 samples repeated essentiality candidate
scans and expensive Groovy array conversion. A fresh per-iteration GO index and
direct key-array construction preserve all 240 recorded prior states (90 with
rescues) across both modes and repeated/changed inputs. The 2,000-EC/32-category
regression reduces key visits from 64,000 to 2,000. Synthetic 100,000-candidate
prior initialization falls from median 4.493 s to 0.072 s; this is not an
end-to-end or biological improvement claim. Existing heuristic floors remain.

Final local validation: 1,010 JVM tests, zero failures/errors, two skips; 525
benchmark tests, zero failures, three optional skips, including eight live
Nextflow checks on the indexed JAR. Frozen full-input jobs: EnzBERT stage1 10079
completed (349,400 candidates, 8 GiB heap/12 GiB allocation), DeepProZyme stage1 10080
(9,790,188 candidates, 64 GiB heap/80 GiB allocation, waits for original 10058).
Original frozen jobs and outputs are unchanged. EnzBERT inference takes 224.56 s
with 3.37 GiB peak RSS; total measurement including verification is 243.40 s.
Independent job 10082 compares every candidate and rehashes 36 recorded files
twice: predictions are byte-identical to the prior complete stage1 run. Both
jobs report COMPLETED/exit 0. Heap and concurrent-load differences prevent
treating this as a controlled end-to-end speedup. Dense/cohort results remain
required before broader capacity acceptance. Evidence archives:
`benchmark/v3/results/claim-key-v1/` and `benchmark/v3/results/essentiality-index-v1/`.


## Supported complex/multifunction construction and combined knockouts (2026-10-26)

Independent measured evidence now supplies explicit partner AND clauses and named
joint activity groups in the guided route. Alternative EC mappings remain separate,
complex-subunit roles cannot silently become monomers, and non-catalytic roles stay
excluded. Deterministic per-group enumeration and final hypothesis caps expose
omissions; no hidden GPR, function confirmation or calibrated score is introduced.
Sixteen added JVM cases cover semantics, roles, aliases, grouping, deduplication,
ordering, truncation and hidden-assignment invariance. The final JVM suite passes
1,026 tests, two skips, zero failures/errors.

Two failing real LP regressions showed that separately applied gene knockouts
could leave OR isozyme alternatives active. All condition deletions and the paired
target deletion now evaluate together. Twenty masked-assignment tests and all 53
metabolic tests pass. The independent COBRA checker also leaked an unused gene's
knockout state across cells. Its two failing isolation/order regressions now pass
with a separate model copy per simulation; three oracle regressions pass in total.
The first failed random-arm check and independent isolated comparisons are retained.

Five new actual CLI runs construct null, complex, second-activity and joint
hypotheses against six conditions. All 200 phenotype/FVA comparisons agree with
independent GLPK, within the shared 256-call cap; changing hidden target assignments
preserves the public problem and every cell. Five earlier simple-fixture CLI runs
also pass all 54 comparisons on the same JAR. These synthetic tests establish
mechanisms, not biological recovery or planner superiority. The full benchmark
suite passes 528 tests with 11 opt-in/optional skips; all eight live Nextflow tests
pass separately on this build. No further JVM changes occurred after the build.

Evidence: `benchmark/v3/results/supported-assignments-v1/` and the
[updated protocol](../docs/benchmarks/v3-masked-assignment-experiments.md).
JAR: `006ecc5339125b5845c1a27cceb175077a4f24975dc02c9f2285f414455664cc`.
Real partner/architecture evidence acquisition, source-level masking/splits,
frozen confidence rules and held-out phenotype/product validation remain open.
No whole item is accepted by this work. Dense-prior job 10058 remains live and
optimized job 10080 remains dependent on it; the all-matrix export is still running.


## InterProScan score contracts and genuine output (2026-10-26)

The parser treated every method statistic as an E-value. Its fallback also turned
malformed statistics into 0.9, and valid profile scores could become negative.
Seventeen focused regressions reproduced score/provenance/batch-failure/thread
problems. The corrected parser distinguishes known E-values from profile and
unreported statistics, preserves original values and semantics, and marks its
legacy score transforms as uncalibrated. The annotation route now honors the
thread budget and exposes the adapter's local-database-only option.

Completed Slurm job 10095 ran installed InterProScan 5.72-103.0 with four analyses,
two public proteins and a synthetic control, rechecking 3,951 bound files. The
actual output reproduces three negative scores in nine old-engine annotations;
all nine identities survive repair, while the six valid E-value-derived scores
remain unchanged. Seven historical TSV parser replays preserve all 379,740
annotation identities and repair 33,758 negative scores. Their historical execution
provenance and biological correctness are not newly established by this replay.

The final JVM state passes 1,046 tests with two skips, zero failures/errors; 19
score/observed-output cases and the pipeline regression pass. The benchmark suite
passes 528 tests with three optional skips, including all eight live Nextflow
checks on JAR c0def1656fc17170869b1d8f25b80e1bf433ecbc583794f63c27aff23b1bc821.
The two CLI test helpers now default to ~/tmp rather than the forbidden RAM tmpfs.
Frozen job 10096 completed with exit 0: actual annotation runs produce nine
identical raw claims, and annotation-plus-integration matches direct integration
byte for byte. All 24 copied output hashes verify independently after transfer. See
[the protocol](../docs/benchmarks/v3-interpro-observed.md) and
`benchmark/v3/results/interpro-observed-v1/` for the observed fixture and evidence.

Per-method confidence remains a legacy heuristic; explicit metadata does not
calibrate the provisional integration model. Other adapters, full route/coverage
failure combinations and independent biological acceptance remain open. No whole
item is accepted. Original dense-prior 10058 and optimized dependency 10080 retain
their frozen engines; no result from them is silently attributed to this JAR.


## Verified shuffle and duplicate-strain controls (2026-10-26)

The historical no-argument perturbation entry point now dispatches to explicit
prepare/run/calibrate/evaluate commands. Frozen preparation preserves all 49,463
proteins and claims, shuffles 49,418 within their original contigs, and creates
five exact copies of a training strain with its original clade/cluster metadata.
Independent checks rehash all inputs/outputs and verify metadata, memberships,
FASTA coverage and coordinate-slot conservation. All 59 production genome/arm
outputs have executed cohort and learned stages and valid receipts.

Calibration uses intact calibration genomes only. EC selects 0.44 (25/27 supported);
GO abstains. Test EC top-1 agreement is 14/37 in each arm. Intact and duplicate EC
calls are 6 supported/6 unverified; shuffle gives 7/5, with net supported change
+1 (paired genome-bootstrap interval −3 to +3). This does not validate a benefit
from intact neighborhoods or transfer of 90% calibration support.

All 156,101 test candidate identities survive duplication. Of these, 114,697 have
numeric probabilities in both arms and 41,404 are unsupported in both; no support
transitions occur. Maximum learned probability change is 0.0017367852562789055.
Raw peer/background counts and duplicate-sensitive permutation features still
reach the fitted model, so exact score invariance remains a repair requirement.
These observations were not used to retune thresholds.

Three pre-fix evaluator cases fail on null probabilities, missing coverage and
malformed/duplicate records. Ten new focused tests and eleven sequence-loss
regressions pass; the full benchmark gate passes 538 tests, three skips, including
eight live Nextflow checks. No JVM changes were made in this slice. Evidence and
source snapshots: `benchmark/v3/results/cohort-controls-v2/`;
[protocol](../docs/benchmarks/v3-cohort-controls.md). No whole item is accepted.


## Strain-cluster permutation correction (2026-10-26)

The label-free duplicate replay exposed a false significance gain: independent
background draws for exact copies reduced the strain cluster's null variance.
One shared random stream per cluster and permutation now keeps corresponding
member draws dependent and prevents copies consuming other units' random draws.
Five before-fix checks fail and pass after repair. The production synthetic replay
now preserves permutation p-value/excess across zero/one/three/five copies. Raw
peer/background count columns still change model margins and require a separate
feature-contract/training repair.

The final JVM gate passes 1,050 tests, two optional skips, zero failures/errors.
The frozen rebuilt JAR is a792849b501e7ad7a9989488ea3b2ce2bd85242de44b0f279b213272bd5a8547.
The complete benchmark gate passes 538 tests, three skips, including eight live
Nextflow checks. The fresh 59-output production replay completed with independently verified
hashes/stages and unchanged original candidate identities. Each arm retains all
484,202 original candidates. Across all 18 genomes, duplicate-probability change
falls from 0.00196057 to 0.00164723; 346,269 numeric pairs and 137,933 unsupported
pairs remain, with no support transitions. No function labels were read and no
policy was refitted. The nonzero residual remains a failed invariance requirement.
Model calibration for corrected features and full duplication invariance remain
unaccepted. [Protocol](../docs/benchmarks/v3-cohort-permutation.md), evidence under
`benchmark/v3/results/cohort-permutation-v1/`.


## Cohort feature contract migration and broader duplicate failure (2026-10-27)

Version 3 removes six legacy unweighted cohort columns and two raw peer/background
count columns. The exporter now reports weighted-context availability, and two
metadata coverage fractions use clade/strain weights. Three failing regressions
now pass; a genuine version-2 fixture is rejected. The rebuilt engine passes
1,055 JVM tests with two optional skips. The complete benchmark gate passes 543
tests, three skips, including nine live Nextflow checks and actual rejection of
an incompatible model. The Python trainer passes 40 tests.

The production re-export preserves all 484,202 candidate identities/source profiles
and non-cohort values, drops exactly eight columns, and changes only repaired
permutation statistics among retained columns. Both models are refitted under the
original training/calibration protocol with all 4,262 supervision rows unchanged.
The 462/608 supported test rows remain a reference-support evaluation; 146 are
withheld. At frozen display thresholds, reference support is 15/32 for logistic and
21/59 for the tree. Both null layers remain untrained; no biological accuracy or
generalization claim follows.

Both models complete all three real controls: 118 independently verified
genome/arm outputs. Across all 18 original genomes, duplication preserves 346,269
numeric probabilities and 137,933 unsupported pairs, with no missing identities
or support transitions. Maximum difference is 2.22e-16 for logistic and zero for
the tree. Separate whole-genome policies select logistic EC 0.44; GO and both tree
policies abstain. Test EC top-1 agreement remains 14/37 for logistic and 13/37 for
the tree in each arm. Logistic intact/duplicate calls are 6 supported/6 unverified;
shuffle is 7/5. The calibration-support failure remains visible.

A broader synthetic replay still fails: when a strain cluster contains differing
observations, exact copies of one member shift the cluster's average. Effective
N stays 1 but smoothed support rises from 0.5 to 0.678571 with five copies. This
requires deduplicating identical evidence within heterogeneous clusters while
preserving distinct members. The real-control success cannot certify general
duplication invariance.

Old engines/models/results remain preserved; current local cohort, sequence-loss
and agent runners default to the version-3 models, with an explicit model-directory
override. Engine ec9be608a606c05fa5d543fddec2434e9d7ce51f185465f7d151a30def33e136.
[Protocol](../docs/benchmarks/v3-cohort-contract.md); evidence/models/source snapshots
under `benchmark/v3/results/cohort-contract-v3/`. No whole item is accepted.

## Heterogeneous-cluster repair and complete export coverage (2026-10-27)

The previous review turn produced fresh test/reproduction evidence and a persistent
review record; it was progress, not a blocked or no-progress turn. The active
objective remains all twelve original requirements.

The complete observed-input identity repair fixes duplicate reweighting inside
heterogeneous strain clusters, consistently across all weighted cohort summaries.
Four cases failed before repair. Twenty-seven added JVM cases cover the failure,
metadata/role/coordinate distinctions, multiple contigs/paralogs, identifier/order
invariance, subset representatives and preservation of every output candidate.
The archived production replay now gives support 0.5 for 0/1/3/5 copies with zero
model-column or attributable margin change. Full JVM gates pass 1,082 tests
(six optional skips); the benchmark gate passes 551 tests (three optional skips,
including nine executed Nextflow cases).

Fresh 18-genome export preserves 484,202 identities/source profiles and all
non-cohort values. Canonical permutation sampling changes 5,043 excess values and
893 p-values. Both models are refitted with all 4,262 supervision rows/splits
unchanged. Test coverage remains 462/608, and null layers remain untrained.
Updated logistic/tree control replays independently verify all 118 outputs. Across
all 18 original genomes, 346,269 numeric probabilities per model are exactly
invariant under duplication; 137,933 unsupported candidates retain their status,
with no missing candidates. Frozen call policies and held-out reference-support
counts are unchanged, including failed calibration transfer. No biological
acceptance claim follows from the correctness repair. See
[the repair protocol](../docs/benchmarks/v3-cohort-evidence-identity.md).

All 132 full Price matrix parts now have verified publication receipts, totaling
1,526,919,081 claims with no count mismatches. Independent collection audit 10109
is rehashing all 851 bound files. Its seven integrity tests pass. Dense
DeepProZyme default-prior job 10058 completed with 9,790,188 candidates preserved,
2:44:34 Java time and 49.91 GiB peak RSS; optimized job 10080 is running and
comparison 10108 is queued. These are capacity/integrity results, not annotation
accuracy or completion of recommendation 1.

## AntiFam quality repair (2026-10-27)

[Typed quality preservation](../docs/benchmarks/v3-antifam-quality.md) addresses review finding 5.
Raw observations have no fabricated functional probability, regional evidence is
preserved, and current learned models abstain on flagged proteins. Training and
label-join readers reject unsupported quality inputs. This does not establish
biological detection accuracy or a learned penalty, and no whole item is signed off.

## Reaction-local missing evidence (2026-10-27)

The [repair record](../docs/benchmarks/v3-rlc-missing-evidence.md) archives
three reproduced regressions, ten before/after production CLI runs, the 1,113-test
JVM gate (two optional skips), and nine passing live Nextflow checks. Missing
predictions contribute no diversity, commitment or self penalty; explicit 0.5
claims are preserved. This closes a scoring defect, not item 5. Genuine gaps and
suggestion-aware biological comparisons remain open.

## Reaction-table gap eligibility (2026-10-27)

The [converter repair](../docs/benchmarks/v3-gap-table-contract.md) groups all
hits before gap selection, excludes strong/unknown/spontaneous/complex evidence,
normalizes EC mappings and rejects malformed tables before opening output. Eleven
of twelve initial CLI regressions fail before repair; all seventeen final checks
pass. The 569-test benchmark gate passes with 12 optional skips. Genuine installed
tool validation, pathway namespaces and subunit-aware gap inference remain open.

## Completed dense-prior and collection checks (2026-10-27)

The [dense-prior archive](../benchmark/v3/results/price-capacity-v4/README.md)
now includes the optimized run and successful independent full-candidate/TSV
comparison. The [collection archive](../benchmark/v3/results/price-full-matrix-v2/README.md)
now includes the successful 851-file audit of all 132 full export parts. These
close the recorded engineering checks, not biological evaluation or larger/cohort
capacity.

## Genuine pathway output and adapter repair (2026-10-27)

The [observed gapsmith record](../docs/benchmarks/v3-gapsmith-observed.md) binds
a complete training-proteome search to the repaired adapter. Job 10111 completed
with exit 0 and unchanged references/inputs. The old parser yields zero claims
from 13,147 rows; the repair retains 3,469 candidates for 802 proteins, separating
EC candidates from reaction/subunit candidates and retaining every raw row.
Nine input and fourteen initial output regressions fail before repair; the final
change adds 34 JVM cases. The full 1,147-test gate passes (two optional skips), as
do nine live Nextflow cases and seven CLI replays. An added failure-before check
also prevents stale observation files when reusing an output directory. Both
learned models explicitly abstain on the unseen source. The complete real-table
gap audit verifies 50 candidate records covering all 46 eligible reactions.

Intermediate build/parser and external staging failures are preserved. These
checks establish acquisition, parsing, retention and route behavior, not pathway
accuracy, calibrated role assignment or the required held-out context comparison.
No whole item is accepted.

## Independent RLC anchor mappings (2026-10-27)

The [mapping repair](../docs/benchmarks/v3-rlc-mapping.md) removes the hidden
pathway-database dependency, retains all supplied EC-to-GO associations, and
prevents correlated labels from multiplying diversity. Ten added JVM cases
pass; the full gates are 1,221 JVM tests (two optional skips), 573 benchmark
tests (twelve optional discovery skips), and nine live Nextflow checks.
Twenty-nine CLI replays cover real mapping recovery, preserved automatic gaps,
and missing-evidence controls. Real standalone mapping now yields the same
820 feature rows and 41 suggestions as the previous pathway-backed route,
with identical functional predictions. The real map is single-valued, so
multiple-association behavior is supported by synthetic regressions only.
RLC remains heuristic; item 5 and every whole audit item remain open.

## Full-cohort pathway acquisition checkpoint (2026-10-27)

The [acquisition protocol](../docs/benchmarks/v3-cohort-pathway-acquisition.md)
freezes all 18 acquired proteomes before searching, without function-label
selection. Array 10112 finished with eight verified bacterial searches, retaining
147,818 reaction rows and 15,384 pathway rows, and ten explicit archaeal preflight
failures. The installed archaeal database was incomplete; real execution also
exposed a directory-marker requirement. Both are covered by added regression
checks. The full benchmark gate passes 591 tests, with twelve optional skips.

Reference acquisition 10120 is publishing 46,477 files from checksum-verified
downloads; its time limit was extended to three hours because shared-storage metadata
I/O is slow. Retry array 10131 waits on `afterok:10120`. The independent audit
still reports incomplete (8/18), and preserves all failed attempts. Complete
search acquisition, production gap replay, reference-overlap audit and matched
context/biological evaluation remain required. No whole item is accepted.

## Partial cohort replay and confirmed reference overlap (2026-10-27)

The [replay and audit](../docs/benchmarks/v3-cohort-pathway-replay.md) completes
32 production CLI calls for all eight verified bacterial searches. It preserves
42,795 raw candidates and checks all 372 eligible reaction groups. Automatic and
explicit gap routes have identical features, gap records and suggestions for
264,687 homology candidates, 389 gap records and 282 suggestions. Functional TSVs
remain byte-identical. These are local-context execution checks, not joint
cohort, trained-model or biological comparisons. Ten archaeal searches still
depend on reference publication; the replay receipt remains `partial_verified`.

Full-sequence audit 10133 confirms 395 target proteins identical to named
reference sequences, including 93 test proteins. This lower-bound finding rules
out treating these searches as proof of unseen-sequence performance. The first
scan job's scheduler failure is retained. A repeated training-genome search also
adds two weak-hit EC candidates with identical query/tool/reference hashes; its
pathway table and 50 production gaps remain unchanged. The cause is unresolved.
No accuracy thresholds or labels are adjusted in response to these diagnostics.

Nine added benchmark regressions pass; the complete gate is 600 tests with
twelve optional skips. No JVM code changed. All replay outputs, exact commands,
hashes, source snapshots and negative findings are preserved under
`benchmark/v3/results/cohort-pathway-replay-v1/`. Full acquisition, reproducibility,
reference independence and matched context/biological acceptance remain open.
No whole audit item is accepted.

## Untestable reaction signatures (2026-10-27)

The [signature-availability repair](../docs/benchmarks/v3-signature-availability.md)
prevents an empty signature, or one smaller than the required anchor count, from
creating negative recurrence evidence. Present families become unknown-context;
family presence/copy observations and layout availability survive. Local context
is withheld from learned inputs when structurally untestable. Testable negative
windows and certified/uncertified family absence remain distinct.

Three initial cohort failures and an intermediate local-bridge failure are
preserved. Eight added JVM cases, ten before/after CLI runs, the 1,229-test JVM
gate, all four live HMMER/AntiFam checks and the 600-test benchmark gate pass.
The latter includes all nine live Nextflow checks with three other optional skips.
Matched eighteen-genome exports preserve all 484,202 identities, non-cohort values
and functional TSV bytes. Context features change for 4,780 candidates; 4,779 lose
unsupported local context and 4,064 lose informative peer evidence. Both models
are refitted with all 4,262 supervision records and original splits preserved.
Null layers remain untrained and test coverage remains 462/608; these models
predict reference support, not independently established enzyme correctness.

Both models' 118 intact/shuffled/duplicated outputs now verify. For each model,
346,269 numeric probabilities are exactly invariant to duplication and 137,933
unsupported candidates remain unsupported. Logistic EC selects 0.45 on calibration
(25/26 supported) but still yields only 6 supported/6 unverified test calls;
shuffle yields 7/5, with net change +1 and interval −3 to +3. Other policies abstain.
Reference independence, actual context ablations, independent locus/role/null
truth and biological validation remain required. No whole item is accepted.
# Three context modes and complete acquisition (2026-10-27)

The [context-mode implementation and comparison](../docs/benchmarks/v3-context-modes.md)
adds universal-network, named-pathway and organism-restricted context to the
shared workflow, CLI and Nextflow. Correlated pathway alternatives remain one
anchor; indirect focal alternatives are excluded. Each organism peer uses its
own restricted graph, including intermediate-node restrictions. Missing input
fails; known empty context remains untestable. Models bind coverage to mode and
abstain on unseen modes. Numerical contract v3 and default universal behavior
are preserved.

The gates pass 1,245 JVM cases (two optional skips), 610 benchmark cases (three
optional skips, ten executed live Nextflow cases), 41 trainer cases and 13
model-input-support cases. Initial implementation failures and a separate
workbook-fixture timestamp defect remain preserved. The first real input
preparation included an empty-stoichiometry reaction that the loader excludes;
the production batch rejected it. Corrected inputs match all 34,703 production
graph IDs.

All 18 pathway searches now verify, including ten completed archaeal retries:
49,463 proteins, 175,272 reaction rows and 19,424 pathway rows. The complete
acquisition does not remove the previously confirmed reference overlap.
Prepared context mappings audit 145 unknown graph IDs and 4,732 groups without
usable mappings. Organism sets describe observed hypotheses, not flux feasibility.

All 72 genome/arm exports verify. Each arm preserves 484,202 candidate/source
identities, non-cohort features and functional TSV bytes. The default control
also reproduces every prior numerical feature. Six separately fitted models
preserve all 4,262 supervision records and splits. Their JVM/Python support and
probabilities agree on 25,572 combinations, with maximum difference 6.11e−16.
Mode-specific supported test coverage is 462/608, 463/608 and 461/608; the common
intersection has 460 rows. All paired genome-level Brier difference intervals
include zero, so this establishes no reliable context benefit. Null models
remain untrained. These are reference-support comparisons on an already examined
cohort, not independent enzyme correctness or biological validation.

Evidence is retained under `benchmark/v3/results/context-modes-v1/`.
Full reference/family exclusions, independent locus/role/null truth, biological
context and organism-feasibility evaluation, complete gap-route replay for the
new archaeal inputs, and the other original requirements remain open. No whole
audit item is accepted.

## Explicit absence comparisons and complete gap replay (2026-10-27)

The [absence comparison](../docs/benchmarks/v3-absence-comparison.md) now runs through
the production metabolic-experiment stage. It compares missing annotation,
unobserved sequence activity, focal-reaction loss and bounded external supply.
Candidate interventions are assignment-specific; they are not whole-gene knockout
predictions. Retained measured EC evidence and validated reaction mapping are
required before testing a candidate. Missing candidates/supply stay unavailable.
Explicit transport, alternative ports, withdrawal, maintenance, incomplete
assemblies, solver errors and whole-condition budgets are tested. Pairwise growth
patterns preserve ambiguity and cannot select an explanation or promote a gene.

Nine actual CLI calls preserve functional TSV and candidate feature bytes. All
147 independent COBRApy/GLPK phenotype comparisons agree. The gates pass 1,251
JVM cases (two optional skips), 613 benchmark cases (three optional skips, ten
executed live Nextflow cases) and 89 metabolic Python cases. Initial test/API
compatibility failures and final passing evidence are preserved in
`benchmark/v3/results/absence-comparison-v1/`.

Separately, the complete 18-genome pathway replay passes 72 production CLI calls:
175,272 reaction rows, 19,424 pathway rows, 56,613 retained raw candidates, 497
eligible reaction groups and 521 gap records (464 with GO). Automatic and explicit
gap routes match across all 484,202 homology candidates, 7,448 RLC rows and 373
suggestions, preserving functional predictions. Evidence is archived under
`benchmark/v3/results/cohort-pathway-replay-v2/`. This uses the prior context-mode
engine and does not evaluate the new absence scenarios on those genomes.

The new absence mechanism uses explicit single-reaction scenarios and synthetic
models; its fragmented variant changes metadata/candidate availability, not real
nucleotide sequences. Real auxotroph and fragmented-assembly evaluation,
independent gene-attribution truth, automatic cohort-gap scenario construction,
community-derived supply and the other original gates remain open. Previously
confirmed pathway-reference overlap remains unresolved. No whole item is accepted.

## Measured nutrient-environment inputs (2026-10-27)

The [new raw-data audit](../docs/benchmarks/v3-nutrient-environments.md) acquires
version 3 of the 2026 nutrient-environment preprint and thirteen files from its
pinned author repository. It reconstructs 16,608 wells into 2,464 assay summaries
across *E. coli* and *B. subtilis*. The 2,240 primary-condition cases include 480
unsupplemented monocultures, 480 supplemented monocultures, 1,200 cocultures and
80 wild-type controls. No-nitrogen controls remain explicit. The missing F+S
replicate and missing measurements are never converted to negative outcomes.

All author-export means reconstruct. The 616 supplemented-export SEM values
incorrectly repeat the mean; the notebook confirms their construction from the
mean columns. Every discrepancy is retained, source files remain unchanged, and
derived uncertainty uses raw replicates with separate pooled-well and
experiment-mean summaries. Simulation inputs contain only identities and assay
conditions; tests prove outcome and coverage changes cannot alter that input set.

Both supplied models map all twelve study genes. Model GPRs preserve the three
affected iML1515 serA reactions and the leuB isozyme bypass instead of forcing
single-reaction losses. Coculture inputs keep organisms separate. These mappings
are suitable for known-genotype phenotype checks, not independent function truth.

Thirteen focused tests and the complete 626-case benchmark gate pass, with three
optional skips and ten executed live Nextflow cases. Sources, models, all derived
versions, exact commands, initial acquisition failure and final tests are retained
in `benchmark/v3/results/nutrient-environments-v1/`. Production JVM code is unchanged.

Flux/media/transport assumptions and matched model runs on these observations
remain to be completed. No phenotype accuracy, independent gene attribution,
transporter knockout validation or whole-item completion is claimed.

## Matched nutrient and supplementation phenotypes (2026-10-27)

The [production benchmark](../docs/benchmarks/v3-nutrient-phenotypes.md) now runs
all 1,144 single-organism assays under six frozen medium/oxygen arms: 6,864
conditions and 13,200 production LPs. Planning and simulation cannot read the
observations. Canonical gene deletion preserves complete GPRs, including the
leuB isozyme bypass and all three affected serA reactions. Exact chemical
exchanges, unrepresented ingredients, maintenance, no-added-source controls,
WT controls, budgets and unknown solver outcomes remain explicit.

The run gives 1,889 growth predictions, 4,974 no-growth predictions and one
unknown. Under the printed medium, E. coli fails to grow in every primary
condition, including 37 measured growing WT assays. Unprinted inorganic supply
improves individual-assay agreement, but every arm in both species remains below
the 80% constant-positive supplementation-rescue comparator. No arm, method,
annotation or learned model is promoted on these results.

Twelve focused regressions and the full 638-test benchmark gate pass, with
three optional skips. The independent GLPK replay exposed solver-history
dependence in two finite objective values; fresh-model diagnosis agrees with
production, and the revised checker resets the solver basis and checks numerical
residuals. The separate HiGHS unknown remains unknown in the frozen predictions.
The final replay verifies all 13,200 GPR comparisons and 13,198 flux comparisons;
one production unknown and one finite discrepancy keep the strict oracle exit status at 1.
Original failures and diagnoses are retained alongside the final verification
in `benchmark/v3/results/nutrient-phenotypes-v1/`.

This validates known-genotype phenotype execution and records its empirical
limitations. It does not supply independent gene-function truth, discriminate
biological absence causes, validate coculture transport or recover a masked
enzyme assignment. Those requirements and whole-item acceptance remain open.

## Production reporting parity (2026-10-27)

The [reporting audit](../docs/benchmarks/v3-reporting-parity.md) reproduces and
repairs three defects: mixed scoring was unknown in the portal; valid complete
model-abstention reports were hidden as untrusted; nested stage records displayed
as `[object Object]`. Failed stages and hash-mismatched/scored abstention reports
remain rejected. Structured values are escaped and preserve unknowns and nulls.

Seven stored production-output sets pass packaged CLI rendering, actual FastAPI
summary/artifact endpoints and Chromium checks. Both models retain all 26,328
rows in each mixed run (7,630 unsupported) and all 3,469 rows in each gapsmith
abstention run. Three synthetic absence-mechanism outputs cover unknown metadata,
insufficient budget and a disabled baseline. The replay adapter was corrected
to retain the complete 2,661-protein input behind the 802 annotated proteins;
initial attempts and the compact-key assertion error remain documented.

The new gates pass 46 service and 73 visualization tests. JAR
`b410dadec00ead589f1f2fd0280ced8eea460755c202c41784f4de3df8a012d5`
differs from `a625288a…` only in `visualize/make_viz.py`; all other entries are
byte-identical. Evidence, source snapshots, screenshots and checksums are in
`benchmark/v3/results/reporting-parity-v1/`. This is reporting correctness on
stored inference, with batch provenance still unavailable. Broader reporting
coverage, biological validation and whole-item acceptance remain open.

## Automatic absence comparisons from run gaps (2026-10-27)

The [guided absence route](../docs/benchmarks/v3-guided-absence.md) now accepts
completed pathway-observation gaps and explicit gap records, retaining each raw
candidate and independently mapped reaction alternative. It rejects same-run
refiner topology and table-only automatic mappings, respects adapter roles,
preserves unavailable candidate/supply states, and shares solver, condition and
time budgets across child comparisons. Every child passes the existing absence
verifier. Failure preserves the verified prefix and stops subsequent work;
reruns clear owned artifacts without deleting unrelated child files.

Fifty focused tests pass. The full JVM suite reports 1,265 tests, zero failures
and six optional skips (including four live-HMMER checks without their optional
environment variable); all 89 Python mechanism tests pass. Six packaged CLI runs
include a baseline and five guided cases. All 100 executed phenotype comparisons
agree with independent GLPK; annotation and feature bytes remain identical to
baseline. Initial implementation/test failures and final receipts are archived
in `benchmark/v3/results/guided-absence-v1/`.

These are synthetic mechanism checks. Queue ordering is an unevaluated heuristic,
peer fractions do not establish target absence metadata, and supplied models are
not automatically independent. Biological absence discrimination, empirical
cohort-wide validation, community-derived supply and calibrated gene recovery
remain open. No whole audit item is accepted.

## Experiment candidate-role repair (2026-10-27)

The [experiment role audit](../docs/benchmarks/v3-experiment-role-bridge.md)
reproduces ignored raw roles in masked-assignment construction. The shared
resolver now preserves raw candidate roles, rejects malformed values, excludes
derived-role overrides and records role provenance for proposed hypotheses.
Subunits require explicit supported partner clauses, including catalytic/subunit
combinations. Explicit and guided absence paths now apply that same standalone
activity restriction; the absence format has no partner clause.

Thirteen regression cases include twelve failing-before cases; all 105 focused
cases pass after repair. The full JVM gate has 1,278 tests, zero failures and six
optional skips. Eighteen matched old/new CLI runs complete as expected and all
90 executed phenotype comparisons agree with independent GLPK. Five masked
cases now abstain, one absence request fails validation, and the valid default
and supported cases preserve their public problems and phenotype cells.
Successful pairs preserve annotation bytes. Initial harness metadata failure,
final results and source snapshots are retained under
`benchmark/v3/results/experiment-role-bridge-v1/`.

These fixtures establish implementation behavior, not biological role or partner
accuracy. Independent role-labelled training/evaluation and calibrated recovery
remain required; no whole item is accepted.

## Matched coculture execution checkpoint (2026-10-27)

The [new community evaluation](../docs/benchmarks/v3-nutrient-community.md)
freezes all 1,320 measured coculture identities, including 120 no-nitrogen
controls, across the same six medium/oxygen arms. Separate genotype/model blocks
preserve complete GPR effects. Every case/arm has matched permissive monocultures
and three community exchange policies under a five-call budget. Pool supplies,
removals, native maintenance and unknown transport support remain explicit.
The bulk endpoint is treated as a proxy, not proof of both strains growing or
measured metabolite exchange.

Seven focused tests and the 645-test benchmark gate pass, with 13 optional skips.
All 60 preliminary real-model LP comparisons agree with an independent GLPK
formulation. The immutable package passes checksums. Slurm array 10150 runs the
7,920 production plans and full independent replay in six bounded shards.
The final plan hash is eaa665f93ea9a2c9b4bbd943b28b441711222c8cc77d5c770ebea59c4a438b46.

Complete numerical verification and endpoint evaluation remain pending. No
whole-item acceptance, transporter validation or gene-attribution claim follows
from this checkpoint.
