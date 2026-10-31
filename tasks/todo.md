# Todo: GSPA improvements

The [supplied GFF3/protein route repair](../docs/benchmarks/v3-supplied-gene-models.md)
now preserves native Prodigal coordinates. Actual baseline 10523 placed none of
4,359 single-mode or 4,286 meta-mode proteins; repaired CLI 10527 places all and
restores 819/813 neighborhood groups. Pyrodigal exports remain byte-identical.
Final 10526 passes 67 focused cases and 1,569 full cases (six optional skips);
all ten actual CLI controls pass. This restores context evidence; it is not an EC
accuracy gain or whole-item completion. All twelve requirements remain open.

Current resumption checkpoint: [exact goal state and next experiments, 2026-10-28](v3-resume-2026-10-28.md).
All twelve whole requirements remain open; the checkpoint distinguishes verified
results, live jobs, and work not yet implemented. The Antigravity audit that followed
it (paralog/role truth cases, cluster state at 13:21) is saved in
[v3-agy-audit-2026-10-28.md; live supervision log with every pi task, cluster job and decision: `tasks/v3-supervision-plan-2026-10-28.md`](v3-agy-audit-2026-10-28.md).

## Audit follow-up — 2026-10-26 (completion reopened)

The [gene-caller repair](../docs/benchmarks/v3-gene-caller-repair.md) is implemented:
bounded processes, typed outcomes, complete-output validation and coordinate-aware
joins now replace the permissive order join. Intermediate build 10517 passes 53
focused and 1,548 full JVM cases (six optional skips). Native Prodigal 2.6.3
observation 10515 supplies actual version/zero/failure/output evidence. Final
build 10520 passes 57 focused and 1,552 full JVM cases (six optional skips).
Actual CLI comparison 10521 passes all 23 invocations: both native callers preserve
valid exports exactly, correct zero/failure/unavailable states and output rejection
are verified, and antiSMASH imports retain 65 identical unscored memberships on
this gene set. The independent 830-file repair archive verifies. This scoped
contract repair does not establish biological accuracy or whole-item acceptance.

The [native gene-caller pre-repair audit](../docs/benchmarks/v3-gene-caller-status.md)
now distinguishes observed successful zero calls from native failure. Jobs 10507,
10509, 10510 and 10511 all exit zero: the unchanged production pipeline misreports
the successful zero run, accepts incomplete output and mismatches reordered
proteins; the proposed strict reader also rejects valid native numeric-contig
output. Independent retained checks reproduce every finding. This is the historical pre-repair finding; the completed repair and validation
are recorded immediately above. The next unverified lead is the supplied-GFF3/
protein join described in the resumption checkpoint. All twelve whole requirements remain open.

The [antiSMASH repair](../docs/benchmarks/v3-antismash-observed.md) now registers
the requested adapter and preserves three native regions and 70 protein
memberships across a 145-record, 4,915-protein captured-output replay. All 124
focused cases and 1,502 JVM cases pass (six optional skips); four actual CLI
comparisons and independent native-feature membership checks pass. Membership
remains unscored and creates no enzyme calls. The [Nextflow follow-up](../docs/benchmarks/v3-antismash-nextflow.md)
now passes 46 focused cases, 1,521 JVM cases (six optional skips), and actual
annotation/import/Nextflow parity on all 70 memberships plus a zero-region sample;
missing native output fails the workflow. Fresh antiSMASH execution, learned BGC
contributions and biological function recovery remain open; no whole item is accepted.

The suspected Nextflow operon ID defect was **not reproduced with pinned
Pyrodigal 3.7.1**: actual single/meta runs retain matching IDs for all 4,359/4,291
proteins and 3,867/3,789 operon members (job 10484). The previous synthetic GFF
used numeric IDs, unlike the actual pinned caller. No ID-normalization patch is
warranted. The comments claiming a three-predictor operon ensemble were inaccurate
and now describe the actual same-strand distance grouping. [Actual Nextflow context
verification](../docs/benchmarks/v3-pyrodigal-operon-route.md) passes in both modes
on this assembly: all 819/812 groups match exhaustive adjacency checks, direct/NF
predictions agree, the controlled in-group candidate changes while other candidates
stay fixed, and numeric-ID controls lose the context effect. All nine command
checks and the portable archive verifier pass. This is not proof of operon or
function accuracy; no whole item is accepted.

The [full six-predictor Price reference comparison](../docs/benchmarks/v3-price-reference-evaluation.md#full-six-predictor-comparison-queued)
is queued as job 10432 after complete collection audit 10410. It retains all
132 cells, 148 eligible proteins and 151 reference functions per predictor.
The unchanged evaluator passes 15 tests; 73 completed-cell metadata checks pass.
Full scoring is pending. All twelve requirements remain open.

The [extracellular-product experiment repair](../docs/benchmarks/v3-extracellular-products.md)
now passes 168 experiment tests, 70 related benchmark tests and all independent
Sco-GEM comparisons (8 standard, 36 masked endpoints and 18 classifications).
Nine new Python checks reproduce the original omissions. The same repair now
extends to evidence-guided Java condition construction. Job 10422 passes 55
focused checks, the full 1,474-case JVM gate (six optional skips), and four actual
CLI runs with independent solver comparisons. These are mechanism checks, not
measured production accuracy or recovered unknown functions. Item 12 remains open.

Robert selected all eleven audit recommendations for implementation and requested
the agent-guided metabolic experiment extension below. The follow-up audit contradicted the blanket completion claim. All twelve original
requirements remain the acceptance scope; see [the current acceptance record](v3-acceptance.md).
The parenthetical claims below are historical implementation reports, not accepted
validation. Audit baseline: commit `7714377`; reviewed release: `498560f`.

The [secondary-product identity audit](../docs/benchmarks/v3-germicidin-identity.md)
now checks the actual primary tables and all four model germicidin products.
Five source-based tests expose a name/formula conflict and ambiguous formula
matches; no model or assay is relabelled. Full planner/function-recovery validation
remains open. Largest Price job 10295 completes both 45,290,880-candidate arms;
collection finalizer 10410 is queued behind the remaining arrays.

A real mixed-source probe found an EC identity defect: 511 shared protein/EC
pairs remained separate because DIAMOND used bare identifiers and CLEAN used
`EC:` prefixes. The ingestion and claim-serialization corrections now pass
cluster regression tests. Corrected real-genome probe 10403 completes with exit
zero and verifies all 511 shared EC assignments. Independent audit 10413 passes
all 145 published output hashes and 178,605 feature rows across nine exports.
The earlier literal-string audit passes do not establish correct evidence combination.
All twelve whole items remain open.

The [five-source context comparison](../docs/benchmarks/v3-source-context-comparison.md)
now has preparation, export, 40-model fitting, production replay and fixed-reference
EC analysis code. Its real probe exposed the EC-prefix defect above and a second
claim-roundtrip defect that loses structured predictor contracts. Ingestion and
serialization corrections pass 30 focused JVM checks and the full 1,473-case
JVM suite (six skips, no failures), with 105 focused Python checks passing.
Runtime preflight 10400 verifies the original parent inputs and numerical runtime.
The original supervision has 6,706 unique canonical identities, with no duplicate
labels to discard.
Full comparison job 10406 is queued behind successful corrected probe 10403
and complete CLEAN audit 10368, with exact prerequisites bound before fitting.
These are integration and scheduling findings, not accuracy results or whole-item
sign-offs.

The [complete EnzBERT reference panel](../docs/benchmarks/v3-price-reference-evaluation.md)
now adds all 22 genomes, 11,553,500 acquired claims and 148 eligible proteins.
Raw/no-prior/default-prior top-1 agreement is unchanged at 53/148; candidate
coverage is 96/151 and both integration comparisons have zero corrections or
harms. Job 10313 passes 15 tests and rehashes all 506 bound files. The independent
41-file archive checks 888 records and the unchanged top-1 identities. This is a
third completed scoped panel, not native predictor accuracy or full-system
validation; all twelve whole requirements remain open.

The [native CLEAN cohort acquisition](../docs/benchmarks/v3-clean-cohort-acquisition.md)
retains all 18 genomes and 49,463 proteins. Setup 10330 installed its environment
but failed before inference because one protein contains J outside ESM's alphabet.
Its dependent jobs were cancelled after that terminal failure. A fresh schema-2
snapshot records an explicit abstention for this protein and runs the other
49,462 unchanged sequences, retaining the original denominator. Thirty focused
checks pass after Sonnet review. New setup 10365 passes all 44 remote checks
and exits zero. Native pilot 10366 completes all 2,994 proteins with 3,695
calls, preserving 801 zero-confidence calls. Array 10367 is running and audit
10368 is pending. Preparation finds 3,730 exact reference-sequence matches, which must remain
visible in the comparison. Completed native results, joint-source feature exports,
training and accuracy remain pending. The multi-source export extension separately
passes 18 new and six compatibility tests. All whole requirements remain open.

The [complete TopEC/MAPred reference evaluation](../docs/benchmarks/v3-price-reference-evaluation.md)
now scores all 44 original-collection cells and 148 eligible proteins. Default
priors produce zero top-1 gains or losses: MAPred remains 50/148 and TopEC 5/148.
Candidate coverage is 114/151 and 5/151 reference functions respectively. Exact
and unambiguous-current EC analyses agree. Job 10306 passes 15 tests and rehashes
all 792 bound files; independent reduced-table checks verify all 1,776 records.
This is reused Price evidence for two predictor panels, not full v3 validation.
All twelve whole requirements remain open.

The [local execution and publication checkpoint](../docs/benchmarks/v3-price-local-execution.md)
passes all 64 remote tests and both complete TopEC pilot comparisons without
changing prediction bytes or non-timing ledger content. Final pilot 10292 uses
the tested compact-audit JAR; all 132 cells and 1,526,919,081 claims are registered
and covered by submitted jobs 10295–10298, gated on dense parity and largest-cell
success. This is scheduling and storage progress, not completed benchmarking.
Original diagnostic 10230 now verifies all 45,290,880 CLEAN candidates in its
no-prior arm. The independent 96-file archive verifies. All twelve items remain open.

The [compact audit storage repair](../docs/benchmarks/v3-compact-audit.md) passes
11 focused tests and the 1,443-case JVM gate (six optional skips). Both complete
TopEC pilot arms retain identical predictions and ledger content apart from
durations. Actual-helper measurements reduce one-million-record audit retention
from 583 MB to 22 MB with identical serialized records; this is not total process
memory or biological accuracy. The 261-file archive verifies independently.
Complete DPZ parity job 10279 passes against both baseline 10235 arms. Predictions
and non-timing ledger content match for all 9,790,188 candidates. Peak RSS drops
15.4%/4.2% in these single observations; default-prior wall time increases. Jobs
10295, 10296 and 10297 have started their dependent production waves.
The largest CLEAN no-prior diagnostic now also verifies every candidate identity.
All twelve whole items remain open.

The [EC identity audit](../docs/benchmarks/v3-ec-identity.md) checks all 8,472
entries in a pinned ENZYME release and all 1,309 legacy replacements. Of those,
101 choose a split branch and six lead to deleted entries; 1,202 resolve uniquely.
Eight tests and complete independent source recomputation pass. No labels, raw
exports, predictions or scores changed. This prevents an unsafe normalization
shortcut in the forthcoming full comparison; it is not an accuracy result.

The [shared ledger-verifier repair](../docs/benchmarks/v3-generic-ledger.md) closes
the generic runner’s remaining full-document read and nonfinite trailing-JSON
acceptance. Local and remote gates each pass 105 cases with one optional JVM skip.
Job 10244 reproduces both dense artifact hashes through both entry points at
24,400 KiB peak RSS. Frozen live benchmark code is unchanged; this is verification
progress, not biological accuracy or whole-item acceptance.

Current Price execution: [full production collection and streaming verification](../docs/benchmarks/v3-price-production-collection.md)
registers all 132 parts and both prior settings. Final replay 10226 preserves the
original dense output/ledger hashes at 23,424 KiB peak RSS; 44 local focused tests
and 24 remote focused tests pass. The first 66-cell wave and largest-CLEAN capacity
diagnostic are running after pilot 10228 verified both arms. These are ongoing runs, not a completed
accuracy comparison; all twelve whole items remain open.

Latest reporting repair: [complete EC reference populations and candidate trace](../docs/benchmarks/v3-ec-candidate-coverage.md)
traces all 335 experimental EC assertions across both regimes and 756,568 raw
claims. Of 39 test assertions, 19 lack a qualifying annotation in the fenced
reference and five lack a matching retained hit; family exclusion removes ten
more. There is no observed integration loss. Schema-2 routine comparison reports
now include every reference protein/function and zero-reference genomes at the
unchanged cutoffs. Job 10224 verifies all 36 panels against the independent audit,
with all 18 cutoffs, 36 original summaries and 48 paired statistics unchanged.
Nine focused tests and the portable 53-file archive pass. This controlled
DIAMOND-source comparison is not a full-ensemble accuracy claim; no predictor or
model was changed, and all twelve whole items remain open.

Latest point-12 implementation: [intracellular product experiments](../docs/benchmarks/v3-intracellular-products.md)
connect explicit synthesis targets to standard and masked evidence-guided routes,
linked condition/deletion experiments, budgets and typed agent feedback. Final
validation passes 280 Python/benchmark tests, 54 focused Java checks, the 1,436-case
Java gate (six optional skips), both CLI routes and six exact legacy-output pairs.
All 12 real standard capacity comparisons, 72 masked phenotype comparisons and
36 demand classifications agree with the independent solver. Six final live
Qwen selections are valid, but the standard agent does not reach a product
deletion and still overinterprets feasible flux. These are controlled known-model
diagnostics, not independent function recovery or demonstrated agent benefit.
All twelve whole items remain unchecked.

Latest deletion repair: [gene identity across FBA and community experiments](../docs/benchmarks/v3-knockout-identity.md)
resolves labels, rejects unknown/ambiguous identifiers and combines all deletions
per organism. All 262 affected tests pass. The [real Sco-GEM diagnostic](../docs/benchmarks/v3-germicidin-diagnostic.md)
verifies 18 native-ID before/after pairs, 90 independent LP comparisons and two
label regressions. Growth misses germicidin loss; maintenance permits synthesis
capacity without growth under phosphate closure. All 24 primary production
observations remain retained but unscored because medium, strain, transport and
congener mappings differ. The complete 81-file archive verifies. Intracellular
agent targets were still open at that checkpoint; independent function-recovery
validation remains open. No whole
item is accepted.

Latest agent-loop repair: [metabolic planner feedback](../docs/benchmarks/v3-planner-feedback.md)
restores prior evidence, product ranges, FVA and missing-result status to the
standard Qwen prompt. All 162 affected tests, six exact non-agent comparisons,
eight live Qwen requests and independent 91-file archive checks pass. The richer
prompt changes test order but increases token cost; Qwen still overinterprets
optional flux and growth-only knockouts. No biological gain is established.
The [secondary-metabolism source intake](../docs/benchmarks/v3-secondary-source-intake.md)
verifies 26 pinned Sco-GEM files and 104 raw table rows, but retains conflicting
strain headers, missing/censored measurements, intracellular-demand semantics
and measured-product fitting leakage. Its 39-file archive verifies; zero cases
are accepted for benchmarking. All twelve whole items remain unchecked.

Current CLEAN audit: [score semantics and integration](../docs/benchmarks/v3-clean-contract.md)
repairs ambiguous native distance/confidence interpretation, shared native caches,
and metadata loss across import/Nextflow. All 21 controlled route invocations pass;
the complete 541-protein native execution preserves all 646 selected calls. Fixed
reference micro-F1 is 0.498 on new and 0.426 on Price, without a GSPA comparison.
Native GMM values remain raw evidence; the existing downstream probability
transform is still heuristic. The complete 768-file archive and every native/import comparison independently
verify. No EC improvement or whole item is accepted.

Latest biological diagnostic: [sulfonate transport and medium dependence](../docs/benchmarks/v3-sulfonate-diagnostic.md)
retains all 247 published measurements and 13 unreported outcomes, with the disputed
genotype explicitly unresolved. The final replay verifies 780 plans, 414 LP pairs
against independent COBRApy/GLPK, two no-sulfur controls, and 14 tests. Six low-yield
butanesulfonate phenotypes remain inconsistent with predicted no growth. Sulfate
alone distinguishes none of the twelve executable genotypes; varying sulfur
sources separates four response patterns but leaves genotype ambiguity. Known
model GPRs, strain mismatch and assumed trace-ion supply prevent a function-recovery
claim. The complete 116-file archive verifies; all twelve items remain open.

Latest dbCAN integration repair: [annotation and Nextflow routes](../docs/benchmarks/v3-dbcan-route.md)
registers the adapter, supports explicit version-4/5 configuration, preserves
machine-readable unscored provenance and rejects unknown protein IDs. The full
Nextflow replay exposed and repaired an empty-builtin-stream abort. All 27 focused
checks and 502 predictor/CLI cases pass (four optional skips); 909 unchanged core
cases passed separately (two skips). Ten controlled route invocations and
independent retained-output checks verify raw-claim parity and missing-output
failure. No real dbCAN search, biological improvement or whole-item acceptance
is claimed; all twelve items remain open.

The [completed context comparison](../docs/benchmarks/v3-context-comparison-results.md)
now retains all twelve fits, 216 production runs and 40,236 table-reader checks.
Independent job 10223 verifies all 6,981 bound files and reproduces the thresholds,
scores and statistics. The registry was already bound by the parent registration;
previous wording claiming otherwise was incorrect. The expanded EC panel includes
all 37 experimental-reference proteins and 39 assertions: only 15 assertions are
candidates with homology and five after family exclusion. At frozen selected
cutoffs, the best arms match 11/39 and 5/39; at 0.5 recovery is much lower. There
is no consistent context advantage or matched CLEAN comparison. The portable
203-file archive verifies all 36 EC panels; all twelve whole items remain open.

Latest genomic contract repair: [strict TSV and geNomad native-output handling](../docs/benchmarks/v3-genomic-contract.md)
rejects malformed/nonfinite records, clears stale full-region records and prevents
missing native summaries from masquerading as zero hits. Thirteen JVM regressions
and fifteen Python subcases fail before repair; six Python methods, 33 focused
checks, 1,384 JVM tests (six optional skips), sixteen CLI/native-process fixture
runs and independent checks pass. Real tool accuracy and remaining native adapters
are not certified.

Latest dispatch repair: [empty-input dispatch and state](../docs/benchmarks/v3-empty-input-dispatch.md)
lets nucleotide predictors run without proteins, records skipped protein scans,
clears previous-genome state and retains requested integration on empty candidates.
All seven regressions fail before the repair and pass afterward. Remote job 10166
passes 88 focused checks and all 1,364 JVM tests (six optional skips). The subsequent
[genomic GFF repair](../docs/benchmarks/v3-genomic-region-output.md) exports original
regions even without CDS and preserves distinct colocated-CDS annotations.
Four intended before/after failures, nine focused checks, 1,369 JVM tests
(six optional skips), twelve production CLI runs and independent retained-output
verification pass. Broader real-tool/route and biological acceptance remain open.

Latest repair: [community result validation](../docs/benchmarks/v3-community-result-contract.md)
reproduces acceptance of missing ledgers, over-budget runs and missing organism
summaries. The repaired stage passes 76 focused tests, 1,356 JVM tests (six
optional skips) and fourteen old/new CLI runs. All malformed reports now fail;
valid predictions, features and raw community results remain byte-identical.
No biological attribution or whole-item acceptance is claimed.

Latest work: [matched joint-cohort inference with observed gaps](../docs/benchmarks/v3-observed-gap-cohort.md)
verifies all 72 genome/arm runs, both complete transfer regimes and 521 observed
gaps. Homology scores change for 202 candidates (four decreases); family-excluded
scores change for 106 (eleven decreases), with sequence likelihoods unchanged.
All four gap-enabled models are frozen and evaluated; all 13,412 Python/JVM table
comparisons agree. All 1,513,136 full-production candidate/model checks now pass, and the complete
matched final-score analysis is verified. Homology rules retain better F1 than
either learned arm; all regressions and unsupported fallback remain reported. The strengthened
report auditor rejects nonfinite probabilities and inconsistent score statuses;
all 544,732 previously archived family-excluded candidate/model records also pass
an independent finite/range/status scan. The benchmark gate passes 716 tests,
with 13 optional skips. Reference-support scores, overlap in search references and
untrained null outcomes still preclude biological or whole-item acceptance. Each
training regime has only two gap-linked supervision rows (both unconfirmed); the
family-excluded test table has no gap-linked rows. Numerical replay cannot
establish learned missing-enzyme recovery on this evidence.

Research intake: the [phosphatase substrate screen](../docs/benchmarks/v3-had-substrate-intake.md)
retains 36,406 measurements across 218 assay sheets. Publication/metadata counts
disagree and source equality/assay/chemical mappings remain unverified. All 218
mirror sequence records now exactly match primary UniProtKB/linked UniParc
sequences. No null or negative GO/EC labels have been accepted from this intake.

Previous progress: [family-excluded training and production verification](../docs/benchmarks/v3-family-excluded-training.md)
fits separate logistic/tree artifacts on 2,444 supervision rows using the verified
family-excluded predictor identity. All 272,366 candidates across 18 genomes are
exported; 4,888 table comparisons and 544,732 full-production candidate/model
comparisons agree. Unsupported candidates abstain; taxon constraints preserve
annotations. The ablation runner now selects the matching regime's model. Test
precision is only 25–29% on supported reference-support rows; null remains
untrained. Full paired comparisons, calibration and independent biological
validation remain open. The final benchmark gate passes 696 tests (13 optional skips).

Completed community evaluation: [all coculture results](../docs/benchmarks/v3-nutrient-community-results.md)
retain all 7,920 plans, 39,600 LP comparisons and 15,840 GPR checks. There are
39,599 numerical agreements, the original solver unknown, and no new discrepancies.
Job 10167 completed the frozen endpoint analysis; independent recomputation
verifies all endpoint records, 192 method summaries and 240 paired summaries.
Printed-medium Bs results underperform constant baselines; Ec predicts no community
growth in those arms. Added-inorganic Bs sensitivity arms improve some endpoints,
with dependence-sensitive uncertainty and no model promotion. Transport, exchanged
metabolite and enzyme-attribution validation remain open.

Previous progress: [observed AMRFinder repair](../docs/benchmarks/v3-amrfinder-observed.md)
reproduces an invented test layout, a real header parsed as a fake protein, zero
identity for every real hit, and ignored pipeline registration. Named-column
parsing, versioned fixtures, explicit configuration/thread forwarding and
unscored-hit exclusion now pass 21 focused cases and seven real CLI runs. All
14 observed calls and route provenance agree; annotation/direct integration is
byte-identical. Long stderr banners no longer hide failure diagnostics. The full
JVM gate passes 1,327 tests, with six optional skips. No whole item is accepted;
independent biological/calibration evidence and remaining adapter coverage are
still required.

Previous progress: [community numerical audit](../docs/benchmarks/v3-nutrient-community-audit.md)
validates the first terminal shard: 6,599 of 6,600 LP comparisons agree, with one
explicit production solver unknown and all 2,640 GPR comparisons agreeing. The
unknown reproduces on a fresh process; alternative solvers validate a zero-growth
solution for the exact production LP, without replacing the frozen unknown.
Sixteen audit regressions and the 683-test benchmark gate pass (13 optional skips).
Interrupted shard 1 retains 6,204 partial comparison rows; job 10158 reruns only
its independent verification. Array 10150 tasks 2–3 run, and 4–5 remain queued.
Full numerical/biological evaluation and all whole-item acceptance remain open.

Previous progress: [predictor identity and score-range support](../docs/benchmarks/v3-predictor-support.md)
now propagates exact predictor metadata through extraction, training and reports.
Eight refits across four context settings preserve all scoring parameters and
4,262 supervision records each; every export preserves 484,202 numerical vectors.
All 40 CLI fault/control runs and 34,096 JVM/Python comparisons pass. Gates pass
1,307 JVM tests, 286 Python script tests and 667 benchmark tests; all ten optional
live Nextflow cases also pass. Legacy identities and out-of-range scores abstain.
Calibration/biological validation and whole-item acceptance remain open; community
array 10150 is still completing independent solver verification.

Previous progress: [predictor-contract reconstruction](../docs/benchmarks/v3-predictor-contract-replay.md)
re-executes DIAMOND and reproduces all 756,568 historical claims across both
transfer modes and 18 genomes. New copies carry verified executable/reference/
parameter/score identities; removing the added field restores the original bytes.
Seven focused tests and the 659-test benchmark gate pass (13 optional skips).
Propagation, learned-model eligibility, score-distribution checks and refitting
remain pending; no model or whole audit item is newly accepted.

Previous progress: [installed eggNOG validation](../docs/benchmarks/v3-eggnog-observed.md)
repairs ignored CLI thread limits and the removed version-3 donor-filter flag.
Fifteen contract tests pass; six final CLI runs on real upstream miniature
databases preserve all 67 version-2 claims/scores and recover 26 version-3 claims.
Both versions have annotation/direct-integration parity and retained provenance.
The full JVM gate passes 1,293 tests with six optional skips. Biological accuracy
and whole-item acceptance remain open; community array 10150 is still running.

Previous progress: the [paired coculture analysis](../docs/benchmarks/v3-nutrient-community-analysis.md)
is frozen separately from the running simulations. Seven focused tests and the
652-test benchmark gate pass (13 optional skips),
including a repaired mutable-freeze defect. Common-case comparisons retain
unknowns; environment/strain resampling and leave-one-strain-out results are
descriptive sensitivity, not calibrated confidence intervals. Full numerical
verification and empirical results remain pending; no whole item is accepted.

Previous progress: [matched nutrient coculture evaluation](../docs/benchmarks/v3-nutrient-community.md)
now freezes all 1,320 cocultures across six arms, with separate genotypes and
five bounded production comparisons per case/arm. Seven focused tests and the
645-test benchmark gate pass (13 optional skips); all 60 preliminary real-model
LP comparisons agree with independent GLPK. Slurm array 10150 is running the
7,920 plans and full independent replay. Complete numerical and empirical
results remain pending; no whole item is accepted.

Previous progress: [experiment role handoff](../docs/benchmarks/v3-experiment-role-bridge.md)
repairs raw role loss in masked-assignment construction and a mixed-role bypass
in absence comparisons. Twelve failing regressions now pass; the full JVM gate
passes 1,278 tests with six optional skips. Eighteen old/new CLI runs reproduce
the defects and preserve valid supported assignments and annotation bytes;
all 90 executed phenotype checks agree with independent GLPK. This verifies
synthetic mechanisms, with biological role accuracy and whole-item acceptance
still open.

Previous progress: [automatic absence comparisons](../docs/benchmarks/v3-guided-absence.md)
now connect run gaps and retained raw EC candidates to separately verified
comparisons under shared budgets. Alternative mappings, unavailable evidence,
failed children and reruns are tested. All 100 synthetic CLI phenotype checks
agree with independent GLPK, preserving annotation and feature bytes. The full
JVM gate passes 1,265 tests with six optional skips; 89 Python mechanism tests
pass. Empirical absence classification, community-derived supply and whole-item
acceptance remain open.

Previous progress: [production reporting parity](../docs/benchmarks/v3-reporting-parity.md)
repairs portal rejection of valid mixed/unsupported-input reports and unreadable
nested stage metadata. Seven packaged CLI/browser/API replays preserve source
identities and scores, with full input-protein denominators. Gates pass 46 service
and 73 visualization tests, including Chromium. The new JAR changes only the
visualization resource. Broader reporting coverage and whole-item acceptance
remain open.

Previous progress: [matched nutrient phenotype validation](../docs/benchmarks/v3-nutrient-phenotypes.md)
now executes all 1,144 single-organism assays under six frozen sensitivity arms
(6,864 conditions), with measured outcomes joined only afterward. The one
inconclusive production solve remains unknown. Paired rescue accuracy remains
below the 80% constant-rescue comparator for both species in every arm; no
method or arm is promoted. Independent numerical verification, its diagnosed
solver-history issue and the final validation receipts are documented in the
linked report. Biological absence, gene attribution, community validation and
whole-item acceptance remain open.

Previous progress: [measured nutrient-environment data](../docs/benchmarks/v3-nutrient-environments.md)
now provide 2,464 outcome-independent assay inputs across two species, including
2,240 primary-condition cases. All 16,608 raw wells are traced; all published
means reconstruct, while an author-export SEM defect is identified in all 616
supplemented-control rows and corrected only in our derived observations.
Missing replicates remain explicit. Both species' supplied models map all twelve
study genes; multi-reaction and isozyme alternatives are preserved. Thirteen
focused tests and the full 626-case benchmark gate pass (three skips, ten live
Nextflow cases). Biological simulation, independent gene attribution and
whole-item acceptance remain open.

Previous progress: [explicit absence and supply comparisons](../docs/benchmarks/v3-absence-comparison.md)
now compare four reaction-level scenarios through the production workflow, with
retained-candidate binding, bounded supply, transport interventions, unknown states
and complete-condition budgets. Nine CLI runs preserve predictions and features;
147 independent GLPK comparisons agree. Gates pass 1,251 JVM cases (two skips),
613 benchmark cases (three skips, ten live Nextflow cases) and 89 metabolic Python
cases. The full 18-genome pathway replay also passes all 72 CLI calls, preserving
56,613 raw candidates and automatic/explicit gap parity across 484,202 homology
candidates. Biological absence/gene-attribution validation, automatic cohort-wide
scenario construction and whole-item acceptance remain open.

Previous progress: [three context modes](../docs/benchmarks/v3-context-modes.md) now
run through the workflow, CLI and Nextflow, with mode-aware model abstention.
The gates pass 1,245 JVM tests (two skips) and 610 benchmark tests (three skips,
ten live Nextflow cases). The complete 18-genome search acquisition verifies.
Four matched export arms each preserve 484,202 candidates and functional TSV
bytes; default mode also reproduces every previous numerical feature. Six models
preserve all 4,262 supervision records; 25,572 JVM/Python scoring comparisons agree.
All paired genome-level Brier intervals cross zero, and null learning remains
unsupported. All 1,281 archived payloads are independently rehashed. Reference independence,
biological validation and whole-item acceptance remain open.

Previous progress: [reaction-signature availability](../docs/benchmarks/v3-signature-availability.md)
now keeps absent or insufficient expected signatures unknown rather than negative.
Eight added JVM cases and ten CLI runs verify the repair; the full gates pass
1,229 JVM cases and 600 benchmark cases, including nine live Nextflow checks.
Matched old/new eighteen-genome exports preserve all 484,202 identities and
functional TSV bytes, correcting 4,780 candidates' context features. Both models
are refitted with all 4,262 supervision records/splits unchanged. All 118 whole-genome
duplicate/shuffle control outputs verify; duplicate probabilities remain exactly
invariant for both models, while calibration transfer still fails. Null learning and biological acceptance
remain open. No whole item is accepted.

Earlier [eight-genome production replay and reference audit](../docs/benchmarks/v3-cohort-pathway-replay.md)
passes all 32 CLI calls, preserving 42,795 raw candidates and matching automatic
and explicit gap features/suggestions across 264,687 homology candidates.
The benchmark gate passes 600 tests, twelve optional skips. Full-sequence checks
confirm 395 target proteins present identically in named references, including
93 test proteins; unseen-sequence performance is not established. Repeated
searches also differ in two weak-hit candidates despite matching input/reference
hashes, with unchanged pathway/gap outputs and an unresolved cause. All findings
and outputs are archived. No whole item is accepted.

Earlier [full-cohort pathway acquisition](../docs/benchmarks/v3-cohort-pathway-acquisition.md)
targets all 18 fixed proteomes (49,463 proteins). All eight bacterial searches
pass the independent audit, with 147,818 reaction rows and 15,384 pathway rows. Real execution
exposed incomplete archaeal references and a directory-marker requirement;
failures are retained and a separate checksum-verified reference acquisition is
running before dependent retry array 10131. Eighteen added checks pass in the 591-test benchmark
gate. Search completion, full-cohort production gap replay and matched context evaluation
remain open; current scheduler and receipt state take precedence over this checkpoint.

Earlier [independent RLC anchor mapping](../docs/benchmarks/v3-rlc-mapping.md)
now honors standalone EC-to-GO input and retains alternatives without counting
correlated labels as independent diversity. Ten added JVM cases and 29 CLI
replays pass. On the real training input, standalone mapping recovers 820 feature
rows and 41 suggestions, exactly matching the earlier pathway-backed route.
Full gates pass 1,221 JVM tests, 573 benchmark tests and nine live Nextflow checks.
These are engineering results; RLC calibration and biological validation remain open.

Earlier [observed pathway-gap integration](../docs/benchmarks/v3-pathway-gap-bridge.md)
now retains both source tables and automatically derives the 50 audited gap
hypotheses through shared integration. Ten real CLI replays verify annotation,
observation-file and explicit-gap feature parity, preserving raw claims and
baseline predictions. Complete/zero/missing/failed cases stay distinct. Full
gates pass 1,211 JVM tests, 573 benchmark tests, and nine live Nextflow checks.
Calibrated context and biological validation remain open; the subsequent RLC
mapping repair is recorded above.

Earlier [local reaction context without peers](../docs/benchmarks/v3-local-reaction-context.md)
now reaches the shared learned/export path with recurrence explicitly unavailable.
Eight real CLI runs on 29,079 unchanged claims and 50 observed gap hypotheses
verify 1,182 restored local contexts, identical sequence evidence, and exact
feature compatibility when independent peers exist. Eleven added JVM and three
Python regressions pass; full gates are 1,172 JVM tests, 572 benchmark tests, and
nine live Nextflow checks. Current models abstain on the new local contexts.
Automatic gap wiring and calibrated/held-out context validation remain open.

Earlier the [predictor-role handoff](../docs/benchmarks/v3-claim-role-bridge.md)
now preserves all 829 observed complex-subunit candidates through learned
candidate construction. Fourteen regressions pass; the old production CLI fails
the new oracle and the repaired engine passes seven CLI replays, the 1,161-test
JVM gate (two skips), and nine live Nextflow checks. Raw claims and functional
TSVs remain byte-identical; both models still abstain on the unsupported source.
Independent role truth and trained role/locus/null inference remain open.

Earlier [genuine gapsmith output and adapter repair](../docs/benchmarks/v3-gapsmith-observed.md)
now preserves all 13,147 observed reaction rows and 3,469 candidates through
production annotation/integration. Input, configuration and score-type defects
are repaired; weak hits and subunit evidence are retained. The 1,147-test JVM
gate, nine live Nextflow checks and seven CLI replays pass. Both learned models
abstain on this new source. A complete-table gap audit yields 50 candidate
records across 46 reactions. Context wiring and biological validation remain open.

Earlier [reaction-local missing evidence](../docs/benchmarks/v3-rlc-missing-evidence.md)
now contributes no artificial probability-0.5 support or penalty. Three new
regressions and ten before/after CLI runs verify the fix; 1,113 JVM tests (two
optional skips) and nine live Nextflow checks pass. The [gap-table converter](../docs/benchmarks/v3-gap-table-contract.md) also
passes 17 CLI regressions after repairing reaction-level eligibility; the broader
569-test benchmark gate passes (12 optional skips). The full Price collection
audit also verifies all 851 bound files across 132 parts. Optimized dense-prior
execution and its independent comparison pass: all 9,790,188 identities and
complete prediction bytes match. No whole item is accepted.

Earlier [heterogeneous-strain evidence identities](../docs/benchmarks/v3-cohort-evidence-identity.md)
repair the reproduced within-cluster duplication defect. Twenty-seven added JVM
cases and the archived production replay pass. Both models have been refitted
with all supervision preserved; all 118 updated production control outputs verify.
Both models now have exactly invariant numeric probabilities across all 18 original
genomes under the real exact-copy intervention. The
1,082-test JVM and 551-test benchmark gates pass (six/three optional skips).
All 132 full Price matrix parts are published and count-verified; a separate
collection hash audit has now passed. No whole item is accepted.

Earlier [InterProScan statistic contracts](../docs/benchmarks/v3-interpro-observed.md)
now distinguish E-values from profile scores and reject malformed statistics.
An actual installed-tool run and seven historical-table parser replays preserve
annotation identity and reproduce/fix negative scores. JVM and live Nextflow
checks pass; the full actual-tool annotation/direct-integration replay also passes. Method-specific
calibration and the other predictor acceptance requirements remain open.

Earlier [supported complex and multifunction assignments](../docs/benchmarks/v3-masked-assignment-experiments.md)
now construct explicit partner AND clauses and evidence-supported joint activities,
with alternative mappings, role restrictions and bounded enumeration. Sixteen new
JVM cases and ten actual CLI runs validate construction and hidden-assignment
invariance. Independent checks agree on 200 new and 54 earlier-fixture phenotype/FVA
comparisons. Two cumulative-knockout solver bugs and an order-dependent independent
oracle defect were reproduced, repaired and regression-tested. These are synthetic
mechanism checks; biological recovery and the other item-12 gates remain open.

Earlier [essentiality indexing and protein identity](../docs/benchmarks/v3-essentiality-index.md)
repair a real pipe-ID rescue crash and remove repeated full-candidate scans.
All 240 synthetic prior states match before/after; the indexed build passes
1,010 JVM tests and 525 benchmark tests. Complete EnzBERT capacity and independent
comparison retain all 349,400 candidates with byte-identical predictions;
dense-prior checks remain running/queued. Earlier [real HMMER execution](../docs/benchmarks/v3-hmmer-observed.md)
exposed and now verifies repairs to availability, malformed statistics and CLI
thread forwarding, including the separate AntiFam adapter. The [AntiFam quality repair](../docs/benchmarks/v3-antifam-quality.md) now preserves
typed raw observations and region metadata through integration and projection.
Current learned models abstain on flagged proteins; a calibrated effect on
functional scores remains open. Annotation and direct integration agree on the observed
synthetic fixture; other adapters' genuine-output coverage remains open.
Verified-batch timeouts now terminate the inference process group, with a
failing-before/passing-after regression archived in `benchmark/v3/results/batch-timeout-v1/`.
The [100-candidate scoring discontinuity](../docs/benchmarks/v3-candidate-count.md)
is fixed and verified through the production CLI; JVM and live Nextflow checks pass.
The [lossless full-matrix export](../docs/benchmarks/v3-price-streaming-export.md)
now independently checks every claim and reproduces the compact artifact's
float32/obsolete-EC statistics without changing raw exports. The complete 132-matrix
run is ongoing; two repaired real cases and complete TopEC/EnzBERT production
capacity checks pass. The 9.79-million-candidate DeepProZyme baseline also passes
exact retention and post-run hash verification; dense priors/cohort measurements
remain pending. [Model input coverage](../docs/benchmarks/v3-model-input-support.md) and the
[GOA source audit](../docs/benchmarks/v3-gaf-negative-audit.md) remain earlier verified
progress. All twelve whole-item acceptance gates remain open.

- [ ] **1. Benchmark the actual application.** Make benchmarks invoke the shared
  production inference engine; generate public claims from versioned experiment
  records. Compare matched sequence evidence with each context stage and the full
  system on held-out genomes/families, with independently selected thresholds,
  calibration, correction/harm counts, and genome-level uncertainty.
  *(Prior implementation report, pending acceptance: `benchmark/v3/run_v3_benchmarks.py`, `benchmark/v3/evaluate_training_cohort_ablations.py`, `benchmark/v3/run_price_context_ablation.py`; 73 runner tests + 429 v3 benchmark tests passing; reports in `docs/benchmarks/v3-*.md`.)*
- [ ] **2. Unify execution routes and prediction state.** Align CLI annotation,
  integration, web, batch, and Nextflow; preserve contributions from combined
  suggesters and synchronize scores, annotations, provenance, and downstream state.
  Require matched-input parity and explicit executed/skipped/failed stage records.
  *(Prior implementation report, pending acceptance: `gspa.integration.workflow.IntegrationWorkflow` shared across `IntegrateCommand`, `AnnotationPipeline.integrate`, web portal, and Nextflow `INTEGRATE`, with schema-2 `<genome>_integrated.tsv.stages.json` stage ledgers and synchronized `Prediction` identity.)*
- [ ] **3. Repair predictor contracts and failure reporting.** Correct eggNOG score
  parsing, HMMER database wiring, genomic-tool dispatch, InterPro EC output, and
  operon GO-aspect handling. Use typed contracts and structured statuses; verify
  genuine output fixtures, complete adapter paths, and injected failures.
  *(Prior implementation report, pending acceptance: `gspa-predictors` and `benchmark/neural/` sidecars updated with calibrated bitscore/e-value parsing, typed adapter statuses `OK`/`SKIPPED`/`FAILED`/`PARTIAL`, genuine output fixtures, and biological failure-injection tests.)*
- [ ] **4. Learn evidence weights and final confidence.** Replace arbitrary boosts,
  candidate-count discontinuities, and irreversible confidence floors with learned,
  calibrated, revisable inference. Retain a no-candidate outcome and account for
  correlated evidence. Compare interpretable and tree-based models with existing
  rules using out-of-fold inputs and family/genus-held-out evaluation.
  *(Prior implementation report, pending acceptance: `LearnedIntegrationStage`, `LearnedContextModel`, `scripts/train_context_model.py`, contract `gspa.context-features/2` with 114 features and explicit `<null>` no-candidate outcome; trained on out-of-fold family-disjoint splits of 18-genome `swissprot_label_rich_v1` (`4,262` rows); held-out test `learned_logistic` AUROC `0.7934`, improving GO $F_{\max}$ from `0.2881` to `0.3221` (`+0.0340`, precision `0.2733 -> 0.3529`); `docs/benchmarks/v3-training-cohort.md`.)*
- [ ] **5. Connect neighborhoods to reaction and interaction networks.** Introduce
  shared protein–reaction candidates with explicit context features, alternative
  mappings, reaction direction, organism-specific feasibility, and distinct
  catalytic/interaction roles. Compare pathway-only, universal-network, and
  organism-specific context on held-out enzyme loci.
  *(Prior implementation report, pending acceptance: `EcCandidateExpander`, `ReactionLocalContextSuggester`, and `ContextFeatureExtractor` bridge reaction-graph distance, directionality, alternative mappings, and roles; evaluated in the 10-arm `swissprot_label_rich_v1` ablation and 42-run Price-148 benchmark; `docs/benchmarks/v3-training-cohort.md`, `docs/benchmarks/v3-price-development.md`.)*
- [ ] **6. Infer recurrence jointly across genomes.** Integrate cohort-wide
  ortholog/subfamily context inference into batch workflows; correct singleton
  normalization, use matched-window backgrounds, account for ancestry and assembly
  completeness, and combine complementary missing-function patterns. Test clade
  holdouts, duplicated strains, shuffled neighborhoods, and fragmented assemblies.
  *(Version-3 replay: 118 production genome/arm outputs across logistic/tree models. Exact-copy probabilities are invariant within 1e-12 over this 18-genome cohort after raw-count removal and permutation repair. A broader heterogeneous-strain-cluster failure is now repaired and regression-tested through complete observed-input identities; the refitted models’ 118 production control outputs verify exact probability invariance across all 18 original genomes. Complementary-gap inference and independent biological acceptance remain open. See `docs/benchmarks/v3-cohort-contract.md`.)*
- [ ] **7. Preserve weak candidates and resolve roles within neighborhoods.**
  Separate candidate retention from display thresholds; distinguish enzymes,
  accessory proteins, transporters, regulators, paralogs, and required complex
  subunits. Evaluate candidate recall and locus ranking with real predictor scores
  and difficult same-operon alternatives.
  *(Prior implementation report, pending acceptance: lossless `--claims-out` and weak-candidate retention down to `1e-9`, `CandidateRole` typing (`CATALYTIC`, `ACCESSORY`, `TRANSPORTER`, `REGULATOR`, `COMPLEX_SUBUNIT`), and same-operon / nearby-locus decoy evaluation on Price-148 across 115,542 proteins; `docs/benchmarks/v3-price-development.md`.)*
- [ ] **8. Model absence, missing sequence, and environmental dependence.** Treat
  missing annotation, incomplete assembly, genuine pathway loss, and external
  supply as competing explanations. Expose medium, biomass, transport, and model
  assumptions. Validate gene attribution separately from growth/knockout prediction
  on complete genomes, auxotrophs, and controlled fragmented assemblies.
  *(Prior implementation report, pending acceptance: explicit stage-ledger FBA medium/biomass/Gram assumptions, competing absence states (`true_absent`, `unobserved_absent`, `pathway_without_family`), 61-case 3-protocol *E. coli* auxotroph validation with independent COBRApy/GLPK cross-verification, and 50% assembly fragmentation control; `docs/benchmarks/v3-metabolic-development.md`, `docs/benchmarks/v3-training-cohort.md`.)*
- [ ] **9. Implement community cross-feeding.** Replace the disconnected heuristic
  scaffold with separate organism models, explicit transport, extracellular
  exchange, and environmental constraints for a defined community. Validate with
  measured coculture/monoculture phenotypes, exchanged metabolites, and transporter
  knockouts; measure reductions in false enzyme assignments separately.
  *(Prior implementation report, pending acceptance: coupled multi-organism SBML FBA in `scripts/community_metabolism.py` with matched no-transfer controls, transport knockouts, and independent COBRApy oracle verification across all 1,035 Wintermute 2010 auxotroph pairs × 6 environmental regimes = 6,210 coupled LPs; `docs/benchmarks/v3-community-development.md`.)*
- [ ] **10. Report prediction changes and uncertainty.** Show stage-wise changes,
  baseline/final scores, evidence, alternatives, and module status; use ontology
  ancestry and EC identifiers for enzyme classification. Reconcile browser counts
  with exports and distinguish internal audit metrics from experimental accuracy.
  *(Prior implementation report, pending acceptance: `scripts/visualize_results.py` and `service/` stage-wise change ledger, GO ancestry and EC hierarchy classification, learned role/null-hypothesis (`<null>`) panels, and headless Chromium browser reconciliation tests.)*
- [ ] **11. Evaluate the agent against equally informed simple methods.** Reserve
  agent use for unresolved cases and useful additional evidence gathering. Compare
  with feature rankings and learned classifiers using identical evidence on unseen
  families/genomes; report correct and harmful decisions, abstentions, and cost.
  *(Prior implementation report, pending acceptance: `benchmark/v3/evaluate_stage3_agent.py` comparing `qwen3.8-27b` (`http://unimatrix01.kaust.edu.sa:8000/v1`, guarded and unguarded) against 6 deterministic and learned comparators on 40 held-out test anomalies of `swissprot_label_rich_v1`; `learned_gbdt` (`+23` net gain, `0.0001 s`) and `learned_logistic` (`+22` net gain, `0.0001 s`) outperform `qwen3.8-27b` (`+15` unguarded / `+7` guarded, `100.0%` precision, `0` harmful, `131.6 s`); `docs/benchmarks/v3-stage3-agent-evaluation.md`.)*
- [ ] **12. Let the agent design and run informative metabolic experiments.**
  Implement the proposal below, building on items 2, 4–8, and 11; use item 9 when
  proposing community exchange experiments.
  *(Prior implementation report, pending acceptance: `scripts/metabolic_experiments.py`, `benchmark/v3/run_planner_benchmark.py`, and `benchmark/v3/audit_condition_dependent_production.py` evaluating `fixed_paired`, `random_item`, `random_paired`, `evidence_item`, `evidence_paired`, and `qwen` (`qwen3.8-27b`) under equal solver budgets $B \in \{8, 16\}$ across 44 genes × 24 substrates: at $B=8$, strict Qwen-attributable discovery is `85.4%` vs `39.0%` for `evidence_paired` (McNemar $p=0.0000$); at $B=16$, `80.5%` strict / `87.8%` hybrid vs `75.6%` for `evidence_paired`, with `75.9%` (`22/29`) binary agreement on visited measured pairs; `docs/benchmarks/v3-planner-evaluation.md`.)*

### 12. Agent-guided metabolic experiments

A single medium and growth objective can leave alternative catabolic pathways and
secondary metabolism unused, so their absence from an FBA solution is weak evidence
against a candidate function. The agent should use sequence/domain evidence,
neighborhoods, cross-genome recurrence, transport predictions, biosynthetic gene
clusters, and existing gaps to shortlist competing functional explanations and the
conditions most likely to distinguish them. It should then design and execute
bounded computational experiments: change carbon/nitrogen sources, oxygen or
nutrient availability, remove supplements, test knockouts, and, where supported,
explore product formation alongside growth or maintenance requirements. Compare
alternative feasible flux solutions and candidate assignments, retain negative or
inconclusive results, and record every assumption. For secondary metabolism, the
relevant reactions must be represented; feasible production is a conditional
hypothesis, not evidence that a pathway is expressed. Choosing an objective that
requires a candidate function must not itself count as independent confirmation.

The expected improvement is recovery of condition-specific functions that a single
growth simulation misses, with fewer unnecessary simulations because prior
evidence guides the search. Benchmark on held-out organisms with measured growth
across media, knockout phenotypes, and characterized secondary-metabolite
production, masking target gene assignments from all model and evidence inputs.
Compare single-medium FBA, a fixed condition panel, random condition selection, a
simple evidence-guided planner, and the agent under the same solver-call budget.
Measure correct gene/function recovery at fixed precision, harmful assignments,
agreement with withheld phenotypes, abstention, and computational cost. Freeze the
planner and confidence rules before testing; count cases where no function can be
recovered rather than evaluating only successful repairs.

Methodological reference: [COBRA Toolbox v3.0: alternative optimal flux solutions
and flux variability analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC6635304/).
Existing starting points: `scripts/stoichiometric_locus_refiner.py` and
`service/metabolic_conditions.py`; current condition probes do not implement this
adaptive experiment-selection loop.

## Historical checklist: GSPA v1.5.0 release

Companion checklist for `tasks/plan.md`. Strike tasks as they complete.

The previous Phase 11 (FM-based operon) work is parked on the
`parking/phase11-glm` branch (NO-GO verdict, 2026-05-05).

## Phase A — Triage and yank

- [x] **A1** `./gradlew clean test` exits 0 (verified 2026-05-06: 335 tests, 0 failures post-yank)
- [x] **A2** Park Phase 11 to `parking/phase11-glm`; remove gLM files; remove `--operon-caller` switch
- [x] **A3** Fix CLI version drift via `VersionProvider` (build.gradle.kts → 1.5.0-SNAPSHOT, picocli reads from generated `version.properties`)
- [x] **A4** Add GPL-3.0-or-later `LICENSE` at repo root; update `README.md` license section

## Phase B — gspa-nf integrator parity

- [x] **B1** New `gspa-nf/modules/integrate.nf` with `BUILD_CLAIMS` + `INTEGRATE` processes
- [x] **B2** Wire into `gspa-nf/main.nf` workflow (guarded by `params.run_integrate`)
- [x] **B3** Add params to `gspa-nf/nextflow.config`
- [ ] **B4** (deferred) Add GO OWL / ec2go / pathways placeholders to `gspa-nf/databases.config` — params already added to `nextflow.config`; databases.config update can ride with the first real cluster run
- [x] **B5** (logical smoke) New jar runs `gspa integrate` on bench9 mgenitalium claims → 4,228 integrated annotations. INTEGRATE process body validated. Full Nextflow Docker smoke deferred (lint pass via `nextflow inspect` already confirmed wiring; dev host has no Docker daemon)
- [x] **B6** Document new flag in `gspa-nf/README.md`

## Phase C — Phase 10 retune (parallel to B)

- [x] **C1** Real gapseq output **already on cluster** (4 of 10 genomes have `*_real.jsonl`; remaining 6 fall back to synthetic 400-gap `*_gaps.jsonl` per the documented zero-byte Reactions.tbl mitigation)
- [~] **C2** Running on unimatrix01 SLURM job 4069 — `phase10_retune.sh` integrates each genome at qBase ∈ {0.50, 0.70, 0.75} plus C1 baseline = 40 integrations
- [~] **C3** Same job scores each (config, tag) with `benchmark_pgap_v2.py` (200 bootstrap, micro + CAFA)
- [ ] **C4** Decide Phase 10 default (default-on iff Δ ≥ 0 with no genome regressing > 0.01)
- [ ] **C5** Update `benchmark/RESULTS.md` Phase 10 section with the verdict

## Phase D — mdF comparison

- [~] **D1** Installing `mDeepFRI` from PyPI into `/data/hohndor/envs/mdf-venv` on unimatrix01 (first attempt stuck on metadata for 13 min, killed; retrying with `--no-cache-dir`)
- [x] **D2** 13-genome PGAP panel: bench9 (9 genomes) + bench10 (10 PGAP genomes) minus duplicates = 13 unique. FAA files confirmed under `/data/hohndor/gspa/proteomes/{tag}.faa` (bench9) and `/data/hohndor/gspa/proteomes/bench10/{tag}.faa` (bench10)
- [ ] **D3** Run mdF on the 13 panel genomes (gated on D1)
- [ ] **D4** Adapter `benchmark/parse_mdf_predictions.py` (with `--self-test`)
- [ ] **D5** Compute F-max micro + CAFA with bootstrap CIs
- [ ] **D6** (Optional) Ensemble: GSPA + mdF substituting for DeepFRI
- [ ] **D7** Add "v1.5.0 — comparison with metagenomic-deepFRI" section to `RESULTS.md`
- [ ] **D8** Cite Bezshapkin et al. 2026 in `README.md`

## Phase E — CHANGELOG and CI

- [x] **E1** `CHANGELOG.md` compiled (Keep-a-Changelog, [Unreleased] entry + v1.0.0 → v1.4.1 history)
- [x] **E2** `.github/workflows/test.yml` (JDK 21, Gradle cache, version smoke)
- [ ] **E3** (Optional, deferred) `.github/workflows/release.yml`
- [ ] **E4** (Optional, deferred) `.github/workflows/docker.yml`

## Phase F — Merge and tag

- [ ] **F1** Pre-merge check: clean test green, shadowJar builds, smoke test passes, RESULTS.md updated
- [ ] **F2** Bump `1.5.0-SNAPSHOT` → `1.5.0`
- [ ] **F3** Open PR `phase11-crossgenome → main`
- [ ] **F4** Tag `v1.5.0` on `main` (annotated)
- [ ] **F5** Branch hygiene (delete `phase11-crossgenome`; preserve `parking/phase11-glm`)
