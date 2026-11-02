# Frozen plan: closing rows 7 (retention and roles) and 11 (agent comparison), 2026-11-01

Status: frozen before any check is run. The sha256 of this file as frozen is recorded in
`benchmark/v3/results/roles-agent-closure-v1/PLAN_SHA256` before the first evaluation. Deviations are
appended, dated, under "Deviations" and never edit the frozen sections.

Scope: the MISSING elements of acceptance rows 7 and 11 in `tasks/v3-acceptance-audit-2026-11-01.md`.
Row 7: "difficult same-operon alternatives" (registry cases carry no operon membership; no same-operon
case is confirmed) and intervals/n for the role benchmark. Row 11: "correctness/harm" is not
identifiable because the unresolved cases of the production anomaly extractor had no matched
experimental truth; the bar also asks for unseen genomes, calibration-frozen decisions, abstention and
complete cost. Established data only (Price-148 reannotated enzymes, Swiss-Prot experimental ECs, the
Swiss-Prot label-rich development cohort, RefSeq model genomes); AllTheBacteria is not used.

Compute: Shaheen III's ControlMaster accepts no new sessions since 07:17 (keepalive log "master down",
`ssh` session open refused); the plan therefore runs on ws (56 cores, 125 GB). Dense Price predictors
(CLEAN 116 GB, GraphEC 113 GB, DeepProZyme 62 GB full-matrix parts; 96–168 GB JVM heaps) do not fit and
are out of scope; the three sparse predictors (TopEC, MAPred, EnzBERT) do. If Shaheen becomes reachable
the dense predictors may be added as a dated deviation; nothing in the criteria depends on them.

## Part A. Row 7: operon membership and same-operon alternatives

### A1. Operon tagging of every registry case

Rule (the production distance rule, `gspa-core/src/main/resources/config/bacteria.yaml`
`predictors.operons.maxIntergenicDistance: 300`, `OperonEnsemble` "distance" predictor): two
consecutive loci on the same contig and strand belong to one operon when the intergenic gap
(start of the downstream gene minus end of the upstream gene minus 1) is at most 300 bp; membership
is transitive along the chain. A strict variant (gap at most 50 bp, `OperonEnsemble` "strict") is
reported as a sensitivity column, never used for the primary tag.

Layouts: the RefSeq layouts already built for the three model genomes
(`/data/hohndor/gspa-v3/role-truth-v1/layouts/{mtb,ecoli,bsubtilis}.layout.tsv`, copied read-only to
`/mnt/data2/leechuck/gspa-v3-data/roles-agent-closure-v1/unimatrix/role-truth-layouts/`) and the
cohort layouts `~/.cache/gspa-v3-data/training-cohort/layout/<assembly>.layout.tsv` for cohort cases.
Each case (and each negative-only locus) receives an `operon_context` record: rule, layout file and its
sha256, the operon group and the strand/contig of every member, and a relation:
`same_operon` (the active member and at least one evidence negative share an operon group),
`different_operon` (all members placed, none shares a group with the active member),
`unplaced` (a member is absent from the layout). For negative-only loci the relation is between the
negatives. The tag is stored as a new optional field on the registry line (the validator ignores
unknown fields; checked before writing) and duplicated in
`benchmark/v3/results/roles-agent-closure-v1/role_truth_operon_tags.jsonl`.

### A2. New same-operon cases

Target: at least three new cases in which the active member and an evidence negative lie in the same
operon under the rule in A1, of the two kinds named in the bar: non-catalytic subunit versus catalytic
subunit, and paralog in the same operon. Candidate loci (chosen from textbook biology before any
score is looked at; the list may shrink but not grow after scores are seen):

| candidate case | genome | function | active | negative | expected operon |
|---|---|---|---|---|---|
| ATCase catalytic vs regulatory chain | E. coli K-12 | GO:0004070 / EC 2.1.3.2 | PyrB (P0A786) | PyrI (P0A7F3) | pyrBI |
| tagatose-bisphosphate aldolase GatY vs GatZ | E. coli K-12 | EC 4.1.2.40 | GatY (P0C8J6) | GatZ (P0C8J8) | gatYZABCDR |
| tagatose-bisphosphate aldolase KbaY vs KbaZ | E. coli K-12 | EC 4.1.2.40 | KbaY (P0AB71) | KbaZ (P42909) | kba/aga cluster |
| biotin carboxylase AccC vs carrier AccB | E. coli K-12 | GO:0004075 / EC 6.3.4.14 | AccC (P24182) | AccB (P0ABD8) | accBC |
| biotin carboxylase AccC vs carrier AccB | B. subtilis 168 | GO:0004075 / EC 6.3.4.14 | AccC (P49787) | AccB (P49786) | accBC |
| IMP dehydrogenase GuaB2 vs GuaB3 (already in registry) | M. tuberculosis H37Rv | GO:0003938 | GuaB2 (P9WKI7) | GuaB3 (P9WKI5) | Rv3411c/Rv3410c |

Label standard (as set by the task): the active member's function and the negative statement must
come from a reviewed Swiss-Prot entry with an experimental evidence code (ECO:0000269 or a descendant)
on the relevant comment or EC, citing a PubMed ID, taken from the pinned Swiss-Prot flat file
(`/home/leechuck/deepgoplus/data/uniprot_sprot.dat.gz`, 687,123,790 bytes, sha256 recorded in the
results folder) and cross-checked against the live UniProt REST record when reachable. Provenance
(accession, comment text, ECO code, PubMed IDs, flat-file release) is recorded per member. A member
whose negative statement is only "by similarity" or absent is labelled `unknown` (label null), never
a negative; a case without an experimentally labelled negative is not added. The proteins' accessions
are confirmed in the RefSeq layout (operon placement) and the claims files (candidate universe). The
candidate universe of each new case is every protein the claims propose for the function at any
score plus the declared members, as for the existing drafted cases. Cases are added with status
`accepted` only when both the positive and at least one negative meet the standard; otherwise
`pending_review`.

### A3. Benchmark rerun with real predictor scores

Arms (all on ws, JAR `gspa-cli/build/libs/gspa-3.0.0.jar` built 2026-11-01 from this worktree, sha256
recorded):

| arm | claims | integration |
|---|---|---|
| L0 | legacy bench9 claims (DIAMOND identity + Pfam; `/data/hohndor/gspa/proteomes/bench9/{mtb,ecoli,bsubtilis}_claims.jsonl`) | `integrate --batch-manifest`, priors none, genome layout, operons from A1 rule |
| L1 | same | default priors, genome layout, operons from A1 rule |
| F0 | role-fenced Swiss-Prot DIAMOND claims with predictor contract (`role-genome-arm-v1/local/arm-inputs/claims`) | priors none, layout, operons |
| F1 | same | default priors, layout, operons |
| S-L, S-F | the same two claim sets | simple method: best single claim score per candidate, no integration (`evaluate_role_truth.py --format claims_jsonl`) |

Learned arms are not run: the only learned models available on ws (`predictor-support-v1/models/default`)
carry the cohort-fenced predictor contract and abstain on both claim sets by design; the role-fenced
models exist only on Shaheen. This is recorded as a limit, not measured as a result.

Evaluator: `benchmark/v3/evaluate_role_truth.py`, threshold 0.5, statuses accepted and pending_review,
registry and negative-only sha256 recorded. For each arm: cases scored of total, strict and lenient
top-1 with Wilson 95% intervals over scored cases, MRR, mean margin with a 2,000-resample bootstrap
percentile interval (seed 20261101) over scored cases, evidence negatives called, negative-only members
called; the same stratified by operon relation (same_operon vs different_operon) and by case kind
(subunit vs paralog). "No-candidate" cases (no member proposed by any claim) are listed per arm with
the count of negatives wrongly proposed (always 0 by construction) and are not counted as scored.

### A4. Pass/fail criteria for Part A

- C7.1 PASS if every registry case and negative-only locus carries an `operon_context` with a relation
  other than `unplaced` for all members present in the layouts (members absent from a layout are listed).
- C7.2 PASS if at least three new cases with relation `same_operon` and status `accepted` are added
  under the A2 standard (at least one of each kind: subunit, paralog). Three `pending_review` cases do
  not pass.
- C7.3 PASS if every arm in A3 completes on the enlarged registry and the report gives, per arm, n
  scored, strict and lenient top-1 with intervals, MRR, mean margin with interval, negatives called,
  negative-only called, and the per-case table includes every same-operon case (scored or no-candidate).
- Reported, not a criterion: whether same-operon cases are ranked correctly (strict top-1 on same-operon
  cases) in each arm; whether default priors with operons change any same-operon ranking relative to
  priors none.

## Part B. Row 11: production anomalies with experimental truth, frozen policies, live agent

### B1. Production anomaly extraction

Inputs: the full-matrix Price claims of the three sparse predictors (`price-full-matrix-v2/parts/{topec,
mapred,enzbert}/<assembly>/claims.jsonl`, sha256 as in each `part.json`), the 22 Price layouts
(`price-full-matrix-v2/layout`, equal to `~/.cache/gspa-v3-data/price-layout-v4`; equality checked),
reference files identical to the Price release-candidate plan (`go.owl`, `ec2go.txt`,
`kegg_pathways.tsv`, `kegg_modules.tsv`; the local copies' sha256 compared with
`/data/hohndor/gspa-v3/price-production-v1/reference`), taxon `NCBITaxon_<taxid>` from the layout report,
essential profile `bacteria`. Configuration = the release-candidate stage1 cell (default priors,
`--genome-layout`), plus `--anomalies-out`. Per predictor three arms exactly as `run_agent_evidence.py`
requires: `stage1`, `learned_logistic`, `learned_gbdt` (the contract-3 models of
`predictor-support-v1/models/default`; whether they score or abstain on these predictors is measured).
Runs go through `v3_runner.run_arm` so that each arm has a content-verified receipt and schema-2 stage
ledgers. Output root `/mnt/data2/leechuck/gspa-v3-data/roles-agent-closure-v1/agent/<predictor>/`.

Dense-label genomes: `run_agent_evidence.py` is rerun on ws for the seven calibration/test genomes of
the development cohort (Swiss-Prot GO labels) with the same JAR; this is the second source and is
evaluated separately (GO truth), not pooled with Price.

### B2. Split

Genus-blocked, fixed before any extraction: NCBI genus from `price-cohort/cohort_metadata.json`
(`clade`), genera sorted, shuffled with `random.Random(20261101)`, genera assigned to calibration in
that order until at least 7 genomes are in calibration; the rest are test. Result:

| split | genus | assemblies |
|---|---|---|
| calibration | ncbi_genus_13687 | GCF_002797435.1 |
| calibration | ncbi_genus_286 | GCF_000698265.1, GCF_001307155.1, GCF_001307275.1, GCF_001440325.1, GCF_001623525.1, GCF_001625455.1 |
| test | ncbi_genus_146937 | GCF_000236665.1 |
| test | ncbi_genus_12916 | GCF_001633105.1 |
| test | ncbi_genus_2901164 | GCF_000327065.1 |
| test | ncbi_genus_1822464 | GCA_004402975.1, GCF_000020125.1 |
| test | ncbi_genus_22 | GCF_000015245.1, GCF_000016065.1, GCF_000146165.2, GCF_000203935.1 |
| test | ncbi_genus_2742 | GCF_000166295.1 |
| test | ncbi_genus_106589 | GCF_000832305.1 |
| test | ncbi_genus_191 | GCF_000237365.1 |
| test | ncbi_genus_309512 | GCF_000018145.1 |
| test | ncbi_genus_106591 | GCA_000006965.1 |
| test | ncbi_genus_302485 | GCA_003443575.1 |

No genus appears in both splits. The six predictors' own training data may contain Price enzymes
(declared control of the Price manifests); the learned integration models were trained on the
development cohort, which shares no genome with Price. The cohort source keeps its own family-blocked
calibration/test split (3 + 4 genomes).

### B3. Truth

A labels table in the `prepare_agent_truth.py` format (genome_id, protein_id, function_id,
swissprot_accession, evidence, match_basis, tier=experimental) built from
(i) Price-148: `price-full-inputs/labels/<assembly>.ec.tsv` (complete ECs only, as exported;
evidence `Price2018-RBTnSeq`, match basis `price_label_audit`), and
(ii) `benchmark/v3/offgrowth/data/swissprot_experimental_ecs_price22.tsv` (ECO:0000269 and
descendants, exact-sequence matched). EC identifiers are normalised to the surface form the anomaly
records use (`EC:` prefix present or absent), nothing else. Truth statuses:
`positive` for (protein, labelled EC); `negative` for (protein, f) when f is a complete EC different
from every labelled EC of that protein and the protein is in the labelled set (the off-growth harm
definition, design §6); every other pair is `unknown`. Partial-EC targets are never matched. The join
is done by `prepare_agent_truth.py` logic (positives) plus the negative rule in a thin wrapper
`benchmark/v3/roles_closure/prepare_price_agent_truth.py`; the wrapper is unit-tested on a fixture.
Truth is read only after the public case bundles are frozen (their sha256 recorded first).

### B4. Methods, calibration and evaluation

The six simple methods of `evaluate_stage3_agent.py suggest` (abstain_only, rank_by_raw_seq,
rank_by_stage1_posterior, rank_by_context_support, learned_logistic, learned_gbdt) and the live agent
(`evaluate_stage3_agent.py agent`: model `qwen3.8-27b`, endpoint `GSPA_LLM_URL` default
`http://unimatrix01.kaust.edu.sa:8000/v1`, key from `~/.config/borg-llm/api-key`, two turns,
temperature 0, seed 42, four workers; the configuration of the existing Stage 3 comparison). Shared
guardrail `shared-reject-v2`. Calibration rule frozen as in the existing comparison: at least ten active
calibration decisions and at least 90% demonstrated precision with unknowns counted against precision;
otherwise abstain on every test case. Policies are written once and never edited.

Agent budget: if the Price case total exceeds 2,000, the agent runs on a seeded (20261101) random
subsample of 2,000 cases stratified by genome, calibration and test alike; every method is then also
evaluated on that subsample so that comparisons are paired. The cohort source (at most 7 × 150 cases)
runs in full.

Reported per method and source: cases, abstained, correct, harmful, unknown, unchanged, error, known
decision precision, demonstrated precision lower bound, net known gain; cost = shared production wall
time, decision wall time, HTTP requests, prompt and completion tokens (agent). Abstention after
calibration is reported separately from pre-calibration actions. Intervals: Wilson 95% for correct
among known decisions when known > 0.

### B5. Pass/fail criteria for Part B

- C11.1 PASS if the extractor runs on all 22 Price genomes for all three sparse predictors with verified
  receipts (66 stage1 cells plus the two learned arms) and on the 7 cohort genomes; the number of cases
  and the 150-cap genomes are reported.
- C11.2 PASS if the truth intersection is measured and reported: number of (genome, candidate protein,
  target function) pairs in the public cases with positive, negative and unknown truth, by split and
  source. The value may be zero; zero is a measured negative result for the bar element
  "correctness/harm" (still not identifiable), not a pass of that element.
- C11.3 PASS if all seven methods produce proposals on every case of their scope, every policy is
  frozen from calibration truth only, and every evaluation reads test truth only.
- C11.4 PASS if correctness, harm, abstention and complete cost are reported per method; the element
  "correctness/harm" counts as evidenced only if at least one method makes at least one known
  (correct or harmful) test decision under its frozen policy.

## Files

- Plan: this file. Code: `benchmark/v3/roles_closure/` (operon tagging, case drafting support, role
  arm runner, interval summary, Price truth wrapper, Price agent source builder) with tests in
  `benchmark/v3/tests/test_roles_closure.py`. Registry edits: `benchmark/v3/role_truth_cases.jsonl`,
  `benchmark/v3/role_truth_negative_only.jsonl` (new field and new lines only; existing lines' fields
  unchanged).
- Results: `benchmark/v3/results/roles-agent-closure-v1/` with `PLAN_SHA256`, `SHA256SUMS`.
- Report: `docs/benchmarks/v3-roles-agent-closure.md`. One progress line in
  `tasks/v3-supervision-plan-2026-10-28.md`. `tasks/v3-acceptance.md` and `README.md` are not edited.

## Deviations

(appended, dated)

- 2026-11-02 (A2): the hand-listed candidate table was replaced by a systematic, score-blind miner
  (`benchmark/v3/roles_closure/mine_operon_role_cases.py`) over the three model genomes and the 18
  development-cohort genomes: an experimental EC on one member, an experimental sentence denying that
  activity on another, operon relation from the layout. Every pair that met the label standard was added;
  none was chosen or dropped by score. Disclosure: during feasibility checks the claim files' candidate
  universes (and hence raw scores) for the candidate functions were looked at before the registry was
  frozen; the selection rule does not use them. PyrB/PyrI (same operon) was drafted and then withdrawn
  because the flat file carries no experimental negative for PyrI; it is listed in the report, not in the
  registry. KbaY/KbaZ was folded into the GatY/GatZ case as a second active/negative pair.
- 2026-11-02 (A1): the strict 50 bp variant is reported per active/negative pair
  (`same_strict_operon`) rather than as a case column.
- 2026-11-02 (B1): the Price arms run on a private frozen copy of the JAR
  (`roles-agent-closure-v1/jar/gspa-3.0.0.jar`, sha256 b9dd0b41...) because the worktree JAR was rebuilt
  by other work during the first attempt and the runner's identity check voided that arm. Learned arms
  whose stage ledger shows `learned_integration` not executed on every genome are recorded as
  "abstained" and the run continues, instead of aborting the collection.
- 2026-11-02 (B): every Price integration carries EC rows only and every anomaly type keys on GO
  terms (`AnomalyExtractor.groovy`: pathway/complex/essential targets are GO; `BORDERLINE_CONFLICT`
  skips non-GO keys), so the extractor emits no case for any Price cell. A diagnostic
  (`price_label_candidate_status.py`) was added to measure where the experimental EC labels sit in the
  integrated posteriors (called / sub-threshold / absent) for what an EC-aware extractor could have seen.
- 2026-11-02 (B, cohort source): the second source (7 cohort genomes with Swiss-Prot GO truth) was not
  rerun: its inputs and code are unchanged from `stage3-agent-evaluation-v2` (196 cases, zero
  experimental pairs), and the split there is already family-blocked; rerunning adds no information.
- 2026-11-02 (B4): the live agent was not exercised because the Price source yields zero cases; the
  endpoint and model remain those of the v2 comparison. Cost of record stays the v2 cost.
