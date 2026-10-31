# Masked condition-dependent function recovery (WP5)

This benchmark asks whether choosing growth experiments helps recover a gene's
metabolic function once that function has been removed from every input. A
gene's reaction assignment is masked in the organism model, in the candidate
evidence and in every text an agent sees. Each arm then receives the same
target-blind candidate assignments (including null) and the same solver budget.
It chooses conditions, simulates the knockout phenotype of every candidate in
them, and the evaluation compares those simulations with measured phenotypes.
Function recovery, phenotype agreement and model-conditional capacity are
scored separately.

The result is negative. Under the calibration-frozen policy no experiment-design
arm recovers any test function at the required precision, because none of them
obtained a usable threshold on the calibration split. Sequence evidence alone
recovers a minority of test cases without harm. The test split is too small to
separate any two arms. Agreement between the production engine and the
independent COBRA replay is an execution check; it is not biological
confirmation, and a simulated phenotype that matches a measurement does not
confirm a function.

Design: `tasks/v3-wp5-masked-recovery-design-2026-10-29.md` (§1–§14, including
the §14 addendum on shared identifiers in masked models). Run log:
`tasks/v3-supervision-plan-2026-10-28.md`, entries from 2026-10-29 onward.
The engine itself is described in [masked assignment experiments](v3-masked-assignment-experiments.md);
the measured nutrient phenotypes come from [the nutrient dataset](v3-nutrient-environments.md)
and [the nutrient phenotype benchmark](v3-nutrient-phenotypes.md).

## Design

### Cases and splits

Three Tier A sources enter; they are gene-disjoint and frozen before any label is joined.

| Split | Source | Organism, model | Panel |
|---|---|---|---|
| Calibration | Keio single-gene auxotrophs (`manuscript_regulatory` protocol) | *E. coli* BW25113, iML1515 | Full 24-substrate grid plus glucose; the measured pairs were chosen by an iML1515 design that used the target GPR |
| Test, primary scope | Nutrient auxotrophs, `printed` medium; sulfonate transport controls (tauA, ssuA) | *E. coli* EcNR1 / iML1515; *B. subtilis* 168 / iYO844; *E. coli* MC4100 / iML1515 | Gene-independent: 44 environments × (none + six supplements); 20 sulfur sources |
| Test, sensitivity scope | Nutrient auxotrophs, `additional_inorganic` medium | as above | as above |

The case checklist (design §3) decides eligibility without reading any outcome.
Truth T(g) is the set of reactions whose original GPR contains the gene. A case
whose truth is mixed or spans more than two reactions is `ambiguous_truth` and is
reported only. Transport-only targets are negative controls: any catalytic call
on them counts as harm. Nutrient cases use a fixed oxygen uptake of
20 mmol gDW⁻¹ h⁻¹; a supplemented condition uses the unsupplemented wild type as
its control. Price RB-TnSeq (Tier B) stays gated on the four conditions of
design §9 and is not part of this run.

Construction produced 70 problems (`audit/protocol.json`):

| Protocol | Problems | Constructed | No admissible experiment |
|---|---:|---:|---:|
| Keio `manuscript_regulatory` (calibration) | 44 | 19 | 25 |
| Nutrient `printed` (test, primary) | 12 | 11 | 1 |
| Nutrient `additional_inorganic` (test, sensitivity) | 12 | 11 | 1 |
| Sulfonate `m63_sulfur_ceiling_1` (test, primary) | 2 | 0 | 2 |
| Total | 70 | 41 | 29 |

A problem with no admissible experiment stays in the denominator of every arm.

### Masking

Every channel in design §4 is masked: the SBML gene product and GPR clauses of
the target, sequence-database self-hits (DIAMOND against a frozen Swiss-Prot with
the target accession and same-genus entries removed), FASTA headers, case files,
condition modifications, product objectives (banned), agent text (opaque
`H<n>`/`C<n>` identifiers and `P<n>` gene aliases, no organism or strain), and
evidence free text. The candidate builder is the production
`MaskedAssignmentBuilder`; every arm, including the exhaustive reference, runs
through the production `run_assignment_experiments` loop on the masked public
SBML, and the engine's public-problem hash must equal the frozen hash.

`audit_wp5_masking.py` gates the run before any LP. It scans every public byte
for target aliases, rebuilds the bundle with the truth moved to decoy reactions,
permuted and with flipped labels, checks that orphan status is not read, checks
the candidate database for excluded hits, and checks that no product objective or
reaction-level condition exists. The masked-model channel scans the masking
boundary structurally (every `fbc:geneProduct` and `fbc:geneProductRef`,
case-sensitive). An alias that the original model already carries on another
gene's record is a shared identifier (operon and complex names on co-subunits,
shared InterPro families); its hits on that other gene are adjudicated and
counted but do not gate. Hits in the target's own record, in any
`geneProductRef`, or of aliases not shared in the original always gate
(design §14 addendum). The whole-file text scan is kept as an informational,
non-gating channel.

| Audit check (`audit/masking_audit.json`, `audit/transcript_audit.json`) | Result |
|---|---|
| Alias scan, gating hits | 0 |
| Alias scan, adjudicated shared-identifier hits (`masked_model_shared`) | 162 |
| Alias scan, informational whole-file hits (`masked_model_text`) | 4,614 |
| Hidden-truth and orphan invariance (decoy GPR, permuted truth, flipped labels, orphan) | 70 problems, 0 violations |
| Candidate database exclusion | pass, 0 violations |
| Objectives and condition modifications | pass, 0 violations |
| Qwen transcript scan (82 transcripts) | pass, 0 hits |

The first audit of this bundle failed (namespace prefix unreadable by the JVM
reader; whole-file alias collisions with non-gene content). That run was
treated as invalid, the writer and the scan were corrected as recorded in design
§14, and the audit was repeated; the failure is retained in the run log.

### Arms and budgets

All arms receive the identical public problem. One condition costs
c = Σ_h (2 + 2·|a_h|) LPs. Each case allows `max_solver_calls = k·c` and
`max_experiments = k` with k ∈ {3, 6}; unused calls are never reallocated.

| Arm | Definition |
|---|---|
| `minimal_medium` | Only the source's minimal medium |
| `fixed_panel` | Frozen source-table order of the gene-independent menu |
| `random` | Seeded shuffle of the menu, seeds 0–9, scored as the per-case mean over seeds |
| `evidence` | Production evidence rule (sum of per-channel maximum hypothesis evidence; ties by identifier) |
| `qwen` | `qwen3.8-27b` on the unimatrix01 endpoint, temperature 0, seed 0, thinking off, 600 s timeout; an invalid reply or endpoint failure stops the arm without fallback credit |
| `sequence_only` (reference) | Top evidence hypothesis, no LP |
| `exhaustive` (reference) | Every menu condition, outside any budget |

### Decisions, policy freeze and metrics

The case decision follows one frozen rule for all arms (design §6). C is the set
of non-null hypotheses that no measured phenotype contradicts. A call is the
unique top-evidence member of C with at least one supporting condition (an
agreeing condition on which the hypotheses' predictions differ). It is correct
if its reactions match T(g) and harmful otherwise. The remaining classes are
null calls (correct negative or false null), tied, unsupported, no consistent
hypothesis, inconclusive (solver or agent failure) and no admissible experiment.

The primary metric is function recovery at fixed precision over **all eligible
test cases**. Each arm's threshold on the decision score was selected on the
calibration split with the rule precision ≥ 0.90 on at least 10 decisions. An
arm without such a threshold abstains on the test split and scores zero.
`evaluate_wp5_recovery.py` joined calibration labels, froze
`evaluation/policy.json`, the hash was reviewed, and only then were test labels
joined; `evaluate` refuses a policy with a different hash. The reported result
is round 2 (policy SHA-256
`e683670a709e6c38cb60cb3a5398974cadb04b25490a9979eee9d94d478e3371`); round 1
(`483b3733e7f132d6deb40cbc2a792d5ecfdb78d45a53aa30f0bbeb88a98340e0`) is retained
in `evaluation-r1/`. Both rounds give every arm the same policy decision; see
[Qwen transport failures and rerun](#qwen-transport-failures-and-rerun). Harm, abstention, false null, correct negative, candidate
recall and conditional recovery are reported beside recovery. Proportions carry
Wilson 95% intervals. Paired tests use exact McNemar and a gene bootstrap with
10,000 resamples (seed 20261029); the random arm enters only the bootstrap.

## Execution and provenance

| Stage | Where | Jobs |
|---|---|---|
| Container smoke test (WP5 unit tests, COBRA fixture) | Shaheen III, k10018 | 15690694 |
| prepare → reference → query → diamond → claims → construct | Shaheen III | 15690699, 15691814, 15691830, 15691909 (failed attempts and their fixes in the run log) |
| construct, masking audit (pass), protocol, exhaustive, deterministic arms | Shaheen III | 15692056 |
| COBRA replay, deterministic arms: one process | Shaheen III | 15692326 |
| COBRA replay, deterministic arms: 256 shards and merge | Shaheen III | 15694094, 15694128–15694142; merge 15696670 |
| Round 1: Qwen arm and transcript audit (8 concurrent workers) | unimatrix01 node004 | 10697 |
| Round 1: COBRA replay, Qwen arms: 82 shards and merge | Shaheen III | 15697662–15697667; merge 15697668 |
| Round 1: join + calibrate (policy frozen) | Shaheen III | 15697684 |
| Round 1: evaluate + report | Shaheen III | 15697736 |
| Round 2: rerun of the six transport-failed Qwen runs (1 worker) and transcript audit | unimatrix01, package `package-r2` | 10700 |
| Round 2: COBRA replay, Qwen arms: shards and merge | Shaheen III | 15701081 |
| Round 2: join + calibrate (policy frozen) | Shaheen III | 15701082 |
| Round 2: evaluate + report | Shaheen III | 15701254 |

The deterministic stages ran from a frozen source snapshot and input package
verified against `SHA256SUMS` in every job. The Qwen endpoint is not reachable
from Shaheen compute nodes, so the Qwen and transcript stages ran on unimatrix01
on the same frozen bundle and deterministic runs, copied from Shaheen with hash
checks; the Qwen runs and transcripts were imported back with hash checks and
no overwrite. The protocol (`audit/protocol.json`,
SHA-256 `b45b3789ceb8d74e100d09061d3d168dfbcf0ec1b99e901065a46acd8ba7bc6c`) binds
the implementation, engine configuration, arms, budgets, seeds and every masked
model hash; the calibration join receipt (`evaluation/join.receipt.json`) binds
the 1,271 run files by hash.

The single-process COBRA replay and the merge of its 256 shards
(`merge_wp5_cobra.py`) produce byte-identical outputs
(SHA-256 `efb68f1e8920a38e73317f99f50d12ba0cabb88a76cd82a997d8a78cac35bd0c`).

Small outputs are mirrored in `benchmark/v3/results/wp5-masked-recovery-v1/`
(`evaluation/`, `oracle/`, `audit/` for round 2; `evaluation-r1/`,
`oracle/cobra-qwen-r1.*`, `audit/transcript_audit-r1.json` and
`qwen-transport-attempt1/` for round 1; all covered by `SHA256SUMS`). The bundle, run
files, transcripts, joined truth and label files stay on Shaheen under
`/scratch/project/k10018/gspa-v3/out/wp5-recovery-v1/`; the unimatrix packages are
`/data/hohndor/gspa-v3/wp5-recovery-v1/package` (round 1) and
`/data/hohndor/gspa-v3/wp5-recovery-v1/package-r2` (round 2).

## Results

All tables are round 2. Relative to round 1 only the Qwen rows change; the
round-1 tables are in `evaluation-r1/report.md`.

### Function recovery at fixed precision (test split)

From `evaluation/report.md`. Denominators include every eligible test case,
including candidate-generation failures; random seeds are averaged per case.

| Arm | Scope | Policy | Threshold | N | Recovery | 95% CI | Harm | Abstention | False null | Correct negative | Candidate recall | Conditional recovery |
|---|---|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|
| evidence:k3 | primary | no_supported_threshold_abstain | – | 13 | 0.000 | 0.000–0.228 | 0.000 | 0.769 | 0.000 | 0.000 | 0.692 | 0.000 |
| evidence:k3 | sensitivity | no_supported_threshold_abstain | – | 11 | 0.000 | 0.000–0.259 | 0.000 | 0.909 | 0.000 | 0.000 | 0.818 | 0.000 |
| evidence:k6 | primary | no_supported_threshold_abstain | – | 13 | 0.000 | 0.000–0.228 | 0.000 | 0.769 | 0.000 | 0.000 | 0.692 | 0.000 |
| evidence:k6 | sensitivity | no_supported_threshold_abstain | – | 11 | 0.000 | 0.000–0.259 | 0.000 | 0.909 | 0.000 | 0.000 | 0.818 | 0.000 |
| exhaustive | primary | no_supported_threshold_abstain | – | 13 | 0.000 | 0.000–0.228 | 0.000 | 0.769 | 0.000 | 0.000 | 0.692 | 0.000 |
| exhaustive | sensitivity | no_supported_threshold_abstain | – | 11 | 0.000 | 0.000–0.259 | 0.000 | 0.909 | 0.000 | 0.000 | 0.818 | 0.000 |
| fixed_panel:k3 | primary | no_supported_threshold_abstain | – | 13 | 0.000 | 0.000–0.228 | 0.000 | 0.769 | 0.000 | 0.000 | 0.692 | 0.000 |
| fixed_panel:k3 | sensitivity | no_supported_threshold_abstain | – | 11 | 0.000 | 0.000–0.259 | 0.000 | 0.909 | 0.000 | 0.000 | 0.818 | 0.000 |
| fixed_panel:k6 | primary | no_supported_threshold_abstain | – | 13 | 0.000 | 0.000–0.228 | 0.000 | 0.769 | 0.000 | 0.000 | 0.692 | 0.000 |
| fixed_panel:k6 | sensitivity | no_supported_threshold_abstain | – | 11 | 0.000 | 0.000–0.259 | 0.000 | 0.909 | 0.000 | 0.000 | 0.818 | 0.000 |
| minimal_medium:k3 | primary | no_supported_threshold_abstain | – | 13 | 0.000 | 0.000–0.228 | 0.000 | 0.769 | 0.000 | 0.000 | 0.692 | 0.000 |
| minimal_medium:k3 | sensitivity | no_supported_threshold_abstain | – | 11 | 0.000 | 0.000–0.259 | 0.000 | 0.909 | 0.000 | 0.000 | 0.818 | 0.000 |
| minimal_medium:k6 | primary | no_supported_threshold_abstain | – | 13 | 0.000 | 0.000–0.228 | 0.000 | 0.769 | 0.000 | 0.000 | 0.692 | 0.000 |
| minimal_medium:k6 | sensitivity | no_supported_threshold_abstain | – | 11 | 0.000 | 0.000–0.259 | 0.000 | 0.909 | 0.000 | 0.000 | 0.818 | 0.000 |
| qwen:k3 | primary | no_supported_threshold_abstain | – | 13 | 0.000 | 0.000–0.228 | 0.000 | 0.615 | 0.000 | 0.000 | 0.692 | 0.000 |
| qwen:k3 | sensitivity | no_supported_threshold_abstain | – | 11 | 0.000 | 0.000–0.259 | 0.000 | 0.545 | 0.000 | 0.000 | 0.818 | 0.000 |
| qwen:k6 | primary | no_supported_threshold_abstain | – | 13 | 0.000 | 0.000–0.228 | 0.000 | 0.615 | 0.000 | 0.000 | 0.692 | 0.000 |
| qwen:k6 | sensitivity | no_supported_threshold_abstain | – | 11 | 0.000 | 0.000–0.259 | 0.000 | 0.000 | 0.000 | 0.000 | 0.818 | 0.000 |
| random:k3 | primary | no_supported_threshold_abstain | – | 13 | 0.000 | 0.000–0.228 | 0.000 | 0.769 | 0.000 | 0.000 | 0.692 | 0.000 |
| random:k3 | sensitivity | no_supported_threshold_abstain | – | 11 | 0.000 | 0.000–0.259 | 0.000 | 0.909 | 0.000 | 0.000 | 0.818 | 0.000 |
| random:k6 | primary | no_supported_threshold_abstain | – | 13 | 0.000 | 0.000–0.228 | 0.000 | 0.769 | 0.000 | 0.000 | 0.692 | 0.000 |
| random:k6 | sensitivity | no_supported_threshold_abstain | – | 11 | 0.000 | 0.000–0.259 | 0.000 | 0.900 | 0.009 | 0.000 | 0.818 | 0.000 |
| sequence_only | primary | selected_on_calibration | 0.864 | 13 | 0.231 | 0.082–0.503 | 0.000 | 0.538 | 0.000 | 0.000 | 0.692 | 0.333 |
| sequence_only | sensitivity | selected_on_calibration | 0.864 | 11 | 0.273 | 0.097–0.566 | 0.000 | 0.636 | 0.000 | 0.000 | 0.818 | 0.333 |

Only `sequence_only` obtained a calibration threshold; its calibration counts
are in the policy table of `evaluation/report.md`. Every experiment-design arm,
the exhaustive reference included, had no threshold reaching 0.90 precision on
at least 10 calibration decisions and therefore abstains on the test split by
rule. The zero recovery of those arms is a consequence of the frozen policy, not
of wrong test calls: no arm made a harmful call, and no catalytic call was made
on a transport control. Candidate recall is identical across arms by
construction, so roughly a third of the primary test cases are lost before any
experiment is chosen.

### Decision classes (test split)

Counts per case from `evaluation/evaluation.json`; random arms are means over ten
seeds. `abstain_policy` marks a call that the frozen threshold withheld. No arm
produced a `wrong_call`, `null_call` or `correct_negative`; those columns are
omitted.

| Arm | Scope | N | correct_call | false_null | abstain_policy | abstain_tied | abstain_unsupported | no_consistent_hypothesis | inconclusive | no_admissible_experiment |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| evidence:k3 | primary | 13 | 0 | 0 | 3 | 3 | 3 | 1 | 0 | 3 |
| evidence:k3 | sensitivity | 11 | 0 | 0 | 7 | 1 | 0 | 2 | 0 | 1 |
| evidence:k6 | primary | 13 | 0 | 0 | 3 | 3 | 3 | 1 | 0 | 3 |
| evidence:k6 | sensitivity | 11 | 0 | 0 | 7 | 1 | 0 | 2 | 0 | 1 |
| exhaustive | primary | 13 | 0 | 0 | 1 | 3 | 3 | 3 | 0 | 3 |
| exhaustive | sensitivity | 11 | 0 | 0 | 2 | 1 | 0 | 7 | 0 | 1 |
| fixed_panel:k3 | primary | 13 | 0 | 0 | 3 | 3 | 3 | 1 | 0 | 3 |
| fixed_panel:k3 | sensitivity | 11 | 0 | 0 | 7 | 1 | 0 | 2 | 0 | 1 |
| fixed_panel:k6 | primary | 13 | 0 | 0 | 3 | 3 | 3 | 1 | 0 | 3 |
| fixed_panel:k6 | sensitivity | 11 | 0 | 0 | 7 | 1 | 0 | 2 | 0 | 1 |
| minimal_medium:k3 | primary | 13 | 0 | 0 | 3 | 3 | 3 | 1 | 0 | 3 |
| minimal_medium:k3 | sensitivity | 11 | 0 | 0 | 7 | 1 | 0 | 2 | 0 | 1 |
| minimal_medium:k6 | primary | 13 | 0 | 0 | 3 | 3 | 3 | 1 | 0 | 3 |
| minimal_medium:k6 | sensitivity | 11 | 0 | 0 | 7 | 1 | 0 | 2 | 0 | 1 |
| qwen:k3 | primary | 13 | 0 | 0 | 0 | 3 | 5 | 0 | 2 | 3 |
| qwen:k3 | sensitivity | 11 | 0 | 0 | 1 | 1 | 4 | 0 | 4 | 1 |
| qwen:k6 | primary | 13 | 0 | 0 | 1 | 3 | 4 | 0 | 2 | 3 |
| qwen:k6 | sensitivity | 11 | 0 | 0 | 0 | 0 | 0 | 0 | 10 | 1 |
| random:k3 | primary | 13 | 0 | 0 | 0 | 3 | 6.9 | 0.1 | 0 | 3 |
| random:k3 | sensitivity | 11 | 0 | 0 | 1.7 | 2.4 | 5.2 | 0.7 | 0 | 1 |
| random:k6 | primary | 13 | 0 | 0 | 0 | 3 | 6.8 | 0.2 | 0 | 3 |
| random:k6 | sensitivity | 11 | 0 | 0.1 | 2.8 | 2 | 3.8 | 1.3 | 0 | 1 |
| sequence_only | primary | 13 | 3 | 0 | 4 | 3 | 0 | 0 | 0 | 3 |
| sequence_only | sensitivity | 11 | 3 | 0 | 4 | 3 | 0 | 0 | 0 | 1 |

Visiting more conditions does not help under this decision rule: the exhaustive
reference ends more often than the budgeted arms with every hypothesis
contradicted by a measured phenotype, most clearly in the sensitivity medium.
The Qwen arm loses most of the sensitivity-scope cases at k = 6 to inconclusive
runs. After the round-2 rerun every inconclusive Qwen run is a model failure
(see below).

### Descriptive test precision–recovery (primary scope)

These curves are computed on the test split and do not select or replace the
frozen operating point. The table gives, per arm, the curve point with the most
decisions; the full curves are in `evaluation/report.md`. Arms without a row made
no call at any threshold.

| Arm | Threshold | Decisions | Precision | Recovery | Harm |
|---|---:|---:|---:|---:|---:|
| evidence:k3, evidence:k6 | 0.430 | 3 | 1.000 | 0.231 | 0.000 |
| fixed_panel:k3, fixed_panel:k6 | 0.430 | 3 | 1.000 | 0.231 | 0.000 |
| minimal_medium:k3, minimal_medium:k6 | 0.430 | 3 | 1.000 | 0.231 | 0.000 |
| exhaustive | 0.762 | 1 | 1.000 | 0.077 | 0.000 |
| qwen:k6 | 0.430 | 1 | 1.000 | 0.077 | 0.000 |
| qwen:k3, random:k3, random:k6 | – | 0 | – | – | – |
| sequence_only | 0.430 | 7 | 1.000 | 0.538 | 0.000 |

Where the experiment arms call on the test split they are correct, but they call
on at most three of 13 cases, and sequence evidence alone calls on more cases at
the same observed precision. With these denominators the curves cannot
distinguish the arms.

### Paired comparisons (primary scope)

| Arm(s) | Comparator | McNemar only-arm | only-comparator | p | Bootstrap 95% recovery difference |
|---|---|---:|---:|---:|---|
| evidence, exhaustive, fixed_panel, minimal_medium, qwen (each budget) | sequence_only | 0 | 3 | 0.250 | -0.462 to 0.000 |
| evidence, fixed_panel, qwen (each budget) | minimal_medium (same budget) | 0 | 0 | 1.000 | 0.000 to 0.000 |
| qwen (each budget) | fixed_panel, evidence (same budget) | 0 | 0 | 1.000 | 0.000 to 0.000 |
| random (each budget) | sequence_only | – | – | – | -0.462 to 0.000 |
| random (each budget) | minimal_medium (same budget) | – | – | – | 0.000 to 0.000 |
| qwen (each budget) | random (same budget) | – | – | – | 0.000 to 0.000 |

No difference is detectable. The three cases recovered only by `sequence_only`
give p = 0.250 against every experiment arm.

### Phenotype agreement

Called-pair agreement compares the simulated phenotype of the called hypothesis
with the measurement, on binary-evaluable pairs; partial growth is excluded.
Because the experiment arms abstain under the frozen policy, almost no arm has a
called pair. The table lists every row of `evaluation/report.md` with called
pairs; all other rows have none.

| Arm | Scope | Set | Binary pairs | Called pairs | Called agree | Called agreement | Coverage | Constant growth | Constant no-growth |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| random:k6 | sensitivity | visited | 150 | 2 | 2 | 1.000 | 0.013 | 0.600 | 0.400 |
| random:k6 | sensitivity | common exhaustive | 7600 | 80 | 44 | 0.550 | 0.011 | 0.533 | 0.467 |
| sequence_only | sensitivity | common exhaustive | 760 | 240 | 235 | 0.979 | 0.316 | 0.533 | 0.467 |

The unmasked-original phenotype comparator of design §7 is not available: the
public-only runner has no stage for it, and no private model was introduced into
the evaluation.

### Model-conditional capacity and cost

C1 is the fraction of eligible cases with at least two hypotheses whose visited
conditions include one that discriminates between them in the masked model. It
measures what the model can distinguish, not recovery. Primary scope, plus the
Qwen sensitivity rows; the full table is in `evaluation/report.md`.

| Arm | Scope | C1 cases | C1 discriminating | C1 fraction | LPs used | LP cap | Unused LPs | Runs | Failed runs | Qwen requests | Qwen failures |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| evidence:k3 | primary | 10 | 4.000 | 0.400 | 182 | 270 | 88 | 14 | 0 | 0 | 0 |
| evidence:k6 | primary | 10 | 4.000 | 0.400 | 364 | 540 | 176 | 14 | 0 | 0 | 0 |
| exhaustive | primary | 10 | 4.000 | 0.400 | 17766 | 26656 | 8890 | 14 | 0 | 0 | 0 |
| fixed_panel:k3 | primary | 10 | 4.000 | 0.400 | 196 | 270 | 74 | 14 | 0 | 0 | 0 |
| fixed_panel:k6 | primary | 10 | 4.000 | 0.400 | 378 | 540 | 162 | 14 | 0 | 0 | 0 |
| minimal_medium:k3 | primary | 10 | 4.000 | 0.400 | 70 | 270 | 200 | 14 | 0 | 0 | 0 |
| minimal_medium:k6 | primary | 10 | 4.000 | 0.400 | 70 | 540 | 470 | 14 | 0 | 0 | 0 |
| qwen:k3 | primary | 10 | 0.000 | 0.000 | 160 | 270 | 110 | 14 | 2 | 33 | 0 |
| qwen:k3 | sensitivity | 10 | 1.000 | 0.100 | 166 | 270 | 104 | 12 | 5 | 33 | 0 |
| qwen:k6 | primary | 10 | 2.000 | 0.200 | 310 | 540 | 230 | 14 | 2 | 60 | 0 |
| qwen:k6 | sensitivity | 10 | 0.000 | 0.000 | 186 | 540 | 354 | 12 | 11 | 41 | 0 |
| random:k3 | primary | 10 | 1.900 | 0.190 | 1806 | 2700 | 894 | 140 | 0 | 0 | 0 |
| random:k6 | primary | 10 | 2.300 | 0.230 | 3584 | 5400 | 1816 | 140 | 0 | 0 | 0 |
| sequence_only | primary | 10 | 0.000 | 0.000 | 0 | 0 | 0 | 14 | 0 | 0 | 0 |

The deterministic budgeted arms reach the same discrimination as the exhaustive
sweep at a small fraction of its LPs; Qwen and random reach less. "Failed runs"
counts runs whose case decision is inconclusive.

### Independent COBRA/GLPK replay

`check_wp5_cobra.py` re-executes every arm-visited cell with COBRApy and GLPK on
the masked public SBML (no sample of unvisited exhaustive cells) and compares
cell contrasts and growth/flux endpoints with the production engine
(`math.isclose`, relative 1e-6, absolute 1e-7). This is execution agreement,
not function or phenotype accuracy.

| Arms | Runs | Shards | Cells agree | Endpoints agree | Merged output SHA-256 |
|---|---:|---:|---:|---:|---|
| Deterministic (minimal_medium, fixed_panel, random, evidence) | 1,066 | 256 | 12,320 / 12,320 | 24,554 / 24,640 | `efb68f1e8920a38e73317f99f50d12ba0cabb88a76cd82a997d8a78cac35bd0c` |
| Qwen, round 2 (reported) | 82 | 82 | 742 / 742 | 1,484 / 1,484 | `f8818f3a542580000fc5315d137fe8be47005636f3f6ce4290f5e323242859e2` |
| Qwen, round 1 (`oracle/cobra-qwen-r1.*`) | 82 | 82 | 703 / 703 | 1,406 / 1,406 | `8ccc88c5d5f394e6b0607c51e25dcac67a7a4182a758f6dd895356f3eaf8ff71` |

All 86 disagreeing endpoints lie in random-arm runs (28 at k = 3, 58 at k = 6)
of the six `additional_inorganic` problems W051–W056, under three condition
identifiers (C242: 56, C298: 28, C284: 2). In every one the COBRA/GLPK growth
is 0.1024002. The run log records the production value as 0.10240007, a
relative difference of 1.3e-6 against the frozen 1e-6 tolerance. No cell
contrast or status differs. The tolerance was not changed after the replay, so
the strict replay exits 1 and is not reported as a complete pass.

### Qwen transport failures and rerun

In round 1, 35 of the 82 Qwen runs ended `agent_failed`. Two causes are
distinguished. A model failure is a reply that re-selects a condition already
visited, which the engine rejects with `ValueError`; under the frozen rule this
stops the arm without credit, and these runs are kept. A transport failure is an
`HTTPError` from the endpoint before any reply could be judged. The six
transport-failed runs (`qwen-transport-attempt1/`) each record one failed HTTP
call with no status code. Five failed on the first request after about 10 s;
the sixth (W035, k = 6) failed on its fifth request after four completed calls.
The same requests succeeded when sent alone. Round 1 used 8 concurrent Qwen
workers instead of the design default of 4, which is the most likely cause. A
transport failure says nothing about the model's choice of experiment, so these
six runs were rerun.

The six runs were archived unchanged (runs, transcripts, the round-1 transcript
audit and its detail, `http-failed.list`) and rerun with one worker on
unimatrix01 (job 10700); every other run file was verified unchanged by hash.
The rerun was decided on the transport evidence before any round-2 outcome was
seen. It does not change the frozen decision rule, the budgets or the agent
failure rule. Calibration and evaluation were then rerun from scratch. This
happened after the round-1 test results had been seen. The policy decision of
every arm is unchanged (only `sequence_only` has a threshold, 0.864); the policy
hash differs only because run files changed. The transcript audit of round 2
passes with 0 hits over 82 transcripts.

| Round | Budget | Qwen runs | Completed | Agent failed | Model (`ValueError`, repeated condition) | Transport (`HTTPError`) |
|---|---|---:|---:|---:|---:|---:|
| 1 | k = 3 | 41 | 28 | 13 | 11 | 2 |
| 1 | k = 6 | 41 | 19 | 22 | 18 | 4 |
| 1 | total | 82 | 47 | 35 | 29 | 6 |
| 2 | k = 3 | 41 | 30 | 11 | 11 | 0 |
| 2 | k = 6 | 41 | 20 | 21 | 21 | 0 |
| 2 | total | 82 | 50 | 32 | 32 | 0 |

Sources: completed and failed totals per round and the round-1 per-budget
totals come from the run files on Shaheen, tabulated by the supervisor; the
transport counts per budget come from `qwen-transport-attempt1/`; the remaining
cells follow from these and the outcomes of the six reruns. Of the six reruns,
three completed (W035 k = 6, W049 k = 3, W050 k = 3) and three failed with
`ValueError` (W045, W047 and W050, all k = 6). The model-failure count therefore
rises in round 2. The mirrored evaluation counts inconclusive cases per arm and
scope, not failure causes; the runs per arm and scope do not sum to 82 here
because `evaluation.json` scores only test cases, while the table counts every
Qwen run including calibration.

In the evaluation, round 2 moves one qwen:k3 case in each scope from
inconclusive to unsupported; no Qwen call is made in either round, so recovery,
harm and the paired comparisons are unchanged.

## Negative results

Against the criteria of design §10:

- **Experiments add no recovery.** No experiment-design arm exceeds
  `minimal_medium` or `sequence_only`; every one recovers nothing at fixed
  precision.
- **No calibration threshold for any experiment arm.** None of them reaches 0.90
  precision on at least 10 calibration decisions, so under this protocol none can
  recover at fixed precision.
- **No agent benefit.** Qwen does not exceed `evidence` or any other comparator
  and has the most inconclusive runs; after the transport rerun all of them are
  model failures (repeated conditions).
- **Candidate generation limits recovery.** Truth is absent from the candidate
  set in roughly a third of the primary test cases, for every arm.
- **Harm.** No arm made a wrong call or a catalytic call on a transport control;
  one random seed in the sensitivity scope produced a false null.

Sequence evidence alone recovers a minority of the test cases with zero harm.
That is the only non-zero fixed-precision result, and its interval is wide.

## Limitations

- The test split has 13 eligible cases (11 in the sensitivity scope). Every
  interval is wide and no paired difference is detectable; a larger labelled test
  set is required before any arm can be ranked.
- The calibration split is small and design-biased: Keio measured pairs were
  chosen by an iML1515 design that used the target GPR (the `note` in
  `evaluation/policy.json`). A policy calibrated there may not transfer to the
  gene-independent test panels.
- Calibration and test differ in source, organism coverage and medium protocol;
  the abstention of the experiment arms on test is driven entirely by calibration.
- Nutrient cases use one oxygen ceiling (20 mmol gDW⁻¹ h⁻¹); no oxygen sweep,
  product objective, reaction-level condition knockout, learned predictor or
  Tier B (Price) case is admitted.
- The phenotype agreement of a called hypothesis is agreement between a
  simulation and a measurement. It does not confirm a function, and the
  unmasked-original comparator is not available.
- The engine–COBRA disagreement on 86 endpoints is a solver-precision effect on
  one growth value; it changes no cell outcome but is retained as a failed
  strict comparison.
- Shared identifiers (operon and complex names, InterPro families on other
  genes) remain in the masked models by design; 162 such hits were adjudicated
  under the §14 ruling, which is recorded as pending Robert's confirmation.

## Reproduction

All stages run on the cluster; local runs are limited to the unit tests
(`benchmark/v3/tests/test_wp5_*.py`).

```bash
# Shaheen III: deterministic stages through the deterministic arms (audit gates every LP)
sbatch benchmark/v3/shaheen/wp5_recovery.sbatch
#   STAGES="prepare reference query diamond claims construct audit protocol exhaustive deterministic cobra" (default);
#   rerun from a later stage with e.g. STAGES="construct audit protocol exhaustive deterministic"

# Shaheen III: COBRA replay fanned out (16 nodes x 16 shards), then merged
for first in $(seq 0 16 240); do
  SHARDS=256 FIRST=$first COUNT=16 ARM_STAGE=deterministic sbatch benchmark/v3/shaheen/wp5-cobra-shards.sbatch
done
SHARDS=256 ARM_STAGE=deterministic sbatch benchmark/v3/shaheen/wp5-cobra-merge.sbatch

# unimatrix01: Qwen arm and transcript audit, resuming from the Shaheen bundle and runs
# (<package>/resume/{bundle,runs} hash-listed in resume/SHA256SUMS)
# round 1 used QWEN_SHARDS=8; the design default is 4
STAGES="qwen transcripts" QWEN_SHARDS=8 \
  sbatch benchmark/v3/wp5_recovery.sbatch /data/hohndor/gspa-v3/wp5-recovery-v1/package
# round 2: the six transport-failed runs only, one worker; all other run files verified unchanged by hash
STAGES="qwen transcripts" QWEN_SHARDS=1 \
  sbatch benchmark/v3/wp5_recovery.sbatch /data/hohndor/gspa-v3/wp5-recovery-v1/package-r2

# Shaheen III: COBRA replay of the Qwen runs after import
SHARDS=82 FIRST=0 COUNT=82 ARM_STAGE=qwen sbatch benchmark/v3/shaheen/wp5-cobra-shards.sbatch
SHARDS=82 ARM_STAGE=qwen sbatch benchmark/v3/shaheen/wp5-cobra-merge.sbatch

# Shaheen III: freeze the policy, review its SHA-256, then join test labels
PHASES="join calibrate" sbatch benchmark/v3/shaheen/wp5-evaluate.sbatch
PHASES="evaluate report" POLICY_SHA256=e683670a709e6c38cb60cb3a5398974cadb04b25490a9979eee9d94d478e3371 \
  sbatch benchmark/v3/shaheen/wp5-evaluate.sbatch
```

The original Qwen-arm COBRA replay ran as six node jobs over the 82 shards; the
`FIRST`/`COUNT` split does not change the merged output. Submit Shaheen jobs with
`--exclude=nid00022` (the wrappers already do where package extraction is
involved). `prepare` refuses to overwrite an existing bundle, and `evaluate`
refuses a policy whose hash differs from `POLICY_SHA256`.

## Archived evidence

`benchmark/v3/results/wp5-masked-recovery-v1/` (11 MB, `SHA256SUMS`):

- `evaluation/`: round 2 `report.md` (all generated tables), `evaluation.json`,
  `report.json`, `policy.json`, `join.receipt.json`; `evaluation-r1/`: the same for round 1;
- `oracle/`: merged COBRA reports for the deterministic arms and the round-2 Qwen
  arms, each with its merge receipt binding every shard hash; `cobra-qwen-r1.*` for round 1;
- `audit/`: `protocol.json`, `masking_audit.json`, `transcript_audit.json` (round 2),
  `transcript_audit-r1.json`;
- `qwen-transport-attempt1/`: the six transport-failed round-1 runs and transcripts,
  `http-failed.list`, and the round-1 transcript audit with its detail.

The archive verifies file integrity and provenance. It is not biological
acceptance of any function call.
