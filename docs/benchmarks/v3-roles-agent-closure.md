# Roles and agent comparison: closing rows 7 and 11 (2026-11-02)

Frozen plan: `tasks/v3-row7-row11-roles-agent-2026-11-01.md` (sha256 of the frozen text
`fbe76850df62cd36…`, copy in `benchmark/v3/results/roles-agent-closure-v1/plan_frozen.md`; dated
deviations are appended to the plan). Results folder: `benchmark/v3/results/roles-agent-closure-v1/`.
Raw runs, inputs and receipts: `/mnt/data2/leechuck/gspa-v3-data/roles-agent-closure-v1/` on ws.
Code: `benchmark/v3/roles_closure/`, tests `benchmark/v3/tests/test_roles_closure.py` (11 tests).

Compute: ws only. Shaheen accepted no session during the work (ControlMaster down since
2026-11-01 07:17; `squeue` unreachable), and the worktree JAR was rebuilt by other work during the
first attempt, so every arm ran on a private frozen copy of the JAR
(`sha256 b9dd0b41b92d78763203d39afa383f73906b8dd14829cc06ff0e33acc328cd3f`, git HEAD at copy time
in `jar/git_head.txt`).

## Verdict against the frozen criteria

| criterion | result |
|---|---|
| C7.1 operon context on every case and negative-only locus | PASS: 15/15 cases and 3/3 negative-only loci tagged; one member unplaced (P9WP43 in the cutinase negative-only group, absent from the RefSeq layout) |
| C7.2 at least three new accepted same-operon cases, both kinds | FAIL (measured): two accepted same-operon cases were added (both non-catalytic subunit vs catalytic subunit); the systematic search over 21 genomes found no third, and no same-operon paralog pair with an experimental negative. The existing GuaB2/GuaB3 case turned out to be same-operon (paralog kind), so the enlarged registry holds three same-operon cases, one of them pre-existing. |
| C7.3 every arm completes, intervals, per-case table | PASS: six arms, n scored, Wilson and bootstrap intervals, strata by operon relation and kind, every same-operon case listed |
| C11.1 extractor on 22 Price genomes × 3 sparse predictors with verified receipts | PASS for the stage-1 arms (66 cells, schema-2 ledgers and content-verified receipts); both learned arms abstained on every cell (engine stage `learned_integration` skipped: "no candidates supported by the model input-coverage manifest"), recorded, not executed |
| C11.2 truth intersection measured | PASS as a measurement: 0 cases, therefore 0 positive / 0 negative / 0 unknown pairs; the cause is structural (below) |
| C11.3 all seven methods under frozen policies | NOT APPLICABLE: no case to decide on; no policy was fitted or applied |
| C11.4 correctness, harm, abstention, cost per method | NOT EVIDENCED: no known test decision exists; the bar element "correctness/harm" remains unmet, now with a mechanism and a quantified bound |

## Part A. Operons and the role benchmark

### A1. Operon tags

Rule: production distance rule (same contig and strand, intergenic gap at most 300 bp, transitive;
`predictors.operons.maxIntergenicDistance` in the bundled presets), computed from the RefSeq layouts
already built for E. coli K-12, B. subtilis 168 and M. tuberculosis H37Rv and from the cohort layouts.
Tags live on every registry line (`operon_context`) and in `role_truth_operon_tags.jsonl` /
`role_truth_negative_only_operon_tags.jsonl`, with the layout sha256, each member's coordinates and
operon group, and a per-pair relation plus the strict (50 bp) variant.

| relation | cases (accepted + pending) | negative-only loci |
|---|---:|---:|
| same_operon | 3 (guaB2/guaB3 mtb; gatY/gatZ ecoli; addA/addB bsub) | 1 (cutinase group mtb, one member unplaced) |
| different_operon | 12 | 0 |
| unplaced (no labelled pair) | 0 | 2 (members not in the cohort layouts' id space) |

### A2. New cases

Search: `mine_operon_role_cases.py` over the pinned Swiss-Prot flat file (release copy
sha256 `a6907c4d…`, the same file the label-rich cohort was built from) for the three model genomes
(4,530 + 6,522 reviewed entries) and the 18 cohort genomes (`cohort_sprot.dat`): an experimental
four-level EC (ECO:0000269) on one protein, an experimental sentence denying that activity on another,
operon relation from the layout. Candidate universes come from the claim files; no score chose a case.

| case | kind | relation | active | negative (experimental statement) |
|---|---|---|---|---|
| gaty-gatz-tagatose-bisphosphate-aldolase-ecoli (GO:0009025) | non-catalytic subunit | same_operon (gatYZ) | GatY P0C8J6 (EC 4.1.2.40, PMID 11976750); KbaY P0AB74 (PMID 10712619), multi-positive | GatZ P0C8J8 "When expressed alone, GatZ does not show any aldolase activity" (PMID 11976750, 8955298); KbaZ P0C8K0 same statement (different operon) |
| adda-addb-3prime5prime-helicase-bsubtilis (GO:0043138) | non-catalytic subunit | same_operon (addAB, gap 13 bp) | AddA P23478 (EC 5.6.2.4, PMID 21071401) | AddB P23477 "does not have helicase activity" (PMID 22307084) |
| thie-teni-thiamine-phosphate-synthase-bsubtilis (GO:0004789) | inactive paralog | different_operon | ThiE P39594 (EC 2.5.1.3, PMID 9139923) | TenI P25053 "Has no thiamine phosphate synthase activity, despite a high sequence similarity with ThiE" (PMID 21534620) |
| hemn-hemw-coproporphyrinogen-dehydrogenase-ecoli (GO:0051989) | inactive paralog | different_operon | HemN P32131 (EC 1.3.98.3, PMID 12114526) | HemW P52062 "Does not have coproporphyrinogen III dehydrogenase activity in vitro" (PMID 29282292) |

Withdrawn: PyrB/PyrI (same operon, pyrBI). PyrI is the regulatory chain, but the flat file carries no
experimental negative for the catalytic activity, so under the label standard it is unknown, and a
case without an experimental negative is not added. SufE/SufS (same operon) was flagged by the miner
and rejected on reading: "does not affect the selenocysteine lyase activity of SufS" is not a
negative about SufE. The cohort genomes yielded no same-operon pair that met the standard.

Registry after the additions: 15 cases (12 accepted, 3 pending), 3 negative-only loci;
`REGISTRY_SHA256` holds the frozen digests
(`f55d00a7…` cases, `d433260a…` negative-only), recorded before any arm was evaluated.

### A3. Arms

Legacy claims (bench9 DIAMOND identity + Pfam) and role-fenced Swiss-Prot DIAMOND claims (predictor
contract, genus-fenced reference), each through `integrate` with the genome layout and the derived
operons (869 / 726 / 850 operons), priors off (`0`) or default (`1`), plus the best single claim per
candidate (S-L, S-F). Threshold 0.5; statuses accepted and pending_review; 15 cases.
The priors arms differ from the no-priors arms (7,596 and 16,714 integrated rows carry fired priors
in L1 and F1), yet no member score of any role locus changed: default priors with operons did not
touch a single reviewed locus.

| arm | scored / 15 | top-1 strict (Wilson 95%) | top-1 lenient | MRR | mean margin (bootstrap 95%) | negatives called | negative-only called |
|---|---:|---|---:|---:|---|---:|---:|
| L0 legacy, no priors | 9 | 0.667 (0.354–0.879) | 0.889 | 0.815 | 0.271 (0.125–0.395) | 4/18 | 0/7 |
| L1 legacy, default priors | 9 | 0.667 (0.354–0.879) | 0.889 | 0.815 | 0.271 (0.125–0.395) | 4/18 | 0/7 |
| F0 fenced, no priors | 6 | 1.000 (0.610–1.000) | 1.000 | 1.000 | 0.202 (0.089–0.317) | 0/18 | 0/7 |
| F1 fenced, default priors | 6 | 1.000 (0.610–1.000) | 1.000 | 1.000 | 0.202 (0.089–0.317) | 0/18 | 0/7 |
| S-L legacy, best claim | 9 | 0.667 (0.354–0.879) | 1.000 | 0.815 | 0.364 (0.177–0.529) | 4/18 | 0/7 |
| S-F fenced, best claim | 6 | 1.000 (0.610–1.000) | 1.000 | 1.000 | 0.329 (0.146–0.513) | 0/18 | 0/7 |

Strata (strict top-1, scored/cases): same_operon L 0.5 (2/3), F 1.0 (1/3); different_operon
L 0.714 (7/12), F 1.0 (5/12). Intervals on one to seven cases are wide; the table reports them, it does
not rank the arms.

Same-operon cases one by one:

| case | legacy (L0 = L1) | fenced (F0 = F1) |
|---|---|---|
| guaB2 vs guaB3 (and guaB1) | three-way tie at 0.667, rank 3 | correct, 0.263 vs 0.132 (margin 0.056) |
| addA vs addB | correct, 0.613 vs 0.262 | AddB not proposed; AddA 0.095 (no-candidate, negative never wrongly proposed) |
| gatY vs gatZ | GatZ and KbaZ never proposed (unscored); GatY and KbaY 0.667 | none of the four proposed |

Reading: when a non-catalytic subunit is proposed at all (AddB by the legacy claims), both claim sets
rank the catalytic subunit above it; the fenced homology claims simply do not propose the non-catalytic
subunits (GatZ, KbaZ, AddB), which is the behaviour the design wants from Stage 0 but leaves the
within-operon ranking untested there. Legacy claims saturate at posterior 0.667 (raw 1.0 identity
hits), producing ties (guaB, folP, coaA/coaX) that the fenced contract does not.

Six cases remain unscored in every arm because the claims never propose their members
(the three archaeal/cyanobacterial cohort cases, and gatY/gatZ in the fenced arms); they are listed
in `role_arms.json` (`unscored_cases`) and counted as no-candidate cases, never as correct.

## Part B. Production anomalies with experimental truth

### B1. What ran

Three sparse Price predictors (TopEC, MAPred, EnzBERT) from the frozen full-matrix parts (66 parts,
every claim digest re-verified against `part.json`), the 22 Price layouts (identical to the Price
release-candidate layouts, sha256 checked), the Price reference files (go.owl, ec2go, KEGG pathways
and modules; sha256 equal to `price-production-v1/reference`), taxon from the layout report,
essential profile bacteria, default priors, `--anomalies-out`; arms stage1, learned_logistic,
learned_gbdt through `v3_runner.run_arm` (receipts and schema-2 ledgers). Genus-blocked split fixed
before extraction (7 calibration genomes in 2 genera, 15 test genomes in 11 genera; no genus in both).

| predictor | stage1 | learned_logistic | learned_gbdt | anomalies (22 genomes) |
|---|---|---|---|---:|
| TopEC | executed, 22/22 verified | abstained 22/22 | abstained 22/22 | 0 |
| MAPred | executed, 22/22 verified | abstained 22/22 | abstained 22/22 | 0 |
| EnzBERT | executed, 22/22 verified | abstained 22/22 | abstained 22/22 | 0 |

### B2. Why zero, and what it bounds

The integrated tables of these cells carry EC rows only (TopEC, GCF_000015245.1: 3,635 rows, all
`EC`, 1,435 called at 0.5). The production extractor (`gspa-core/.../agent/AnomalyExtractor.groovy`)
derives pathway, complex and essential anomalies from the set of annotated GO terms and explicitly
skips non-GO keys for borderline conflicts, so a single-predictor EC collection cannot yield an
unresolved case whatever its scores. The Price release-candidate collection, the only production
collection with experimental EC truth on new genomes, is therefore outside the extractor's reach by
construction. Supplying operons would not change this (operons only decorate candidates of GO
targets).

How far away the truth sits even for an EC-aware extractor (`price_label_candidate_status.py`,
candidate floor 0.05, annotated 0.5):

| source of label | labels | called | sub-threshold (0.05–0.5) | absent from the integration |
|---|---:|---:|---:|---:|
| Price-148 complete ECs, TopEC stage1 | 151 | 2 | 3 | 146 |
| Swiss-Prot experimental ECs, TopEC stage1 | 35 | 3 | 2 | 30 |
| Price-148 complete ECs, MAPred stage1 | 151 | 61 | 53 | 37 |
| Swiss-Prot experimental ECs, MAPred stage1 | 35 | 24 | 2 | 9 |
| Price-148 complete ECs, EnzBERT stage1 | 151 | 41 | 55 | 55 |
| Swiss-Prot experimental ECs, EnzBERT stage1 | 35 | 21 | 2 | 12 |

At most a handful of experimentally labelled (protein, EC) pairs would have been visible as candidates
at all; the remaining labels are not proposed by the predictor at any score. The missing bar element
therefore needs (i) a GO-bearing production claim set on genomes with experimental truth and (ii) a
predictor that proposes the labelled functions below threshold, neither of which the current
collections provide.

### B3. Methods, policies, agent

No case exists, so no simple method and no agent produced a proposal, no calibration policy was
fitted, and no test truth was read. The truth wrapper (`prepare_price_agent_truth.py`, unit-tested)
is in place for a future GO-bearing collection. The cohort source (7 label-rich genomes, Swiss-Prot
GO truth) was not rerun: inputs and code are unchanged from `stage3-agent-evaluation-v2`, whose
196 cases met no experimental assertion. Cost of record for the agent remains that run's.

## Files

- `benchmark/v3/results/roles-agent-closure-v1/`: `PLAN_SHA256`, `plan_frozen.md`, `REGISTRY_SHA256`,
  `role_truth_operon_tags.jsonl`, `role_truth_negative_only_operon_tags.jsonl`, `role_arms.json`,
  `price_label_status_*.json`, `price_sources_*.json`, `SHA256SUMS`.
- Registry: `benchmark/v3/role_truth_cases.jsonl` (+4 cases, `operon_context` on every line),
  `benchmark/v3/role_truth_negative_only.jsonl` (`operon_context` added).
- Code: `benchmark/v3/roles_closure/{operon_context,mine_operon_role_cases,add_operon_cases,
  run_role_arms,run_price_anomalies,prepare_price_agent_truth,price_label_candidate_status}.py`.
