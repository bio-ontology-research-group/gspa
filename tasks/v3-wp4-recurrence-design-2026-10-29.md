# WP4 design: recurrence and complementary gaps on independent clade holdouts (item 6)

Status: **design only, nothing registered or run.** Requirement: handoff `tasks/v3-resume-2026-10-28.md`
§8 item 3; acceptance row 6 of `tasks/v3-acceptance.md` stays open whatever the outcome. Builds on
`benchmark/v3/results/peer-selection-v1/` (111 peer rows, 105 unique accessions), the 10406
`diamond_homology` plans/models, `prepare_context_modes.py`, `prepare_role_genome_context.py`,
`prepare_sequence_loss.py` and `prepare_cohort_controls.py`. The 2026-10-29 role-genome arm lessons are
binding: new rows need their real layout path and `unplacedProteins = FASTA − layout`, and claims must
carry the **same predictor contract** as the model's training claims, or the models abstain.

## 1. Question and what can be claimed

The 18 cohort genomes are 18 genera, so every earlier recurrence feature came from distant genomes. WP4
asks whether, for the same sequence candidates, adding unlabelled close relatives to the joint cohort
(a) concentrates the context evidence of a homolog on its experimentally labelled function, and (b) does so
through complementary pathway pieces carried by different relatives. It also asks whether that changes the
recovery of frozen learned models. Claims are limited to these genomes, this engine configuration and
DIAMOND/Swiss-Prot candidates. Context concentration is not evidence that a function is correct.

**Interpretation (decision D1).** The three context relations are the arms *pathway-only* (`PATHWAY`),
*universal-reaction* (`UNIVERSAL`) and *organism-model* (`ORGANISM`); `SEQONLY` is the candidate anchor.
The *joint-recurrence* arm is the peer regime crossed with each relation: R0 = cohort only (the 10406
setting), R2 = cohort + the target family's peers. R2 − R0 within a relation is the joint-recurrence effect.

## 2. Holdout unit and genomes

Holdout unit = **NCBI family block**: target + all its `peers.tsv` rows. Tiers: ND = near-duplicate (same
species), C = close (same genus), F = same family/other genus, D = the 17 other cohort genomes. Peer sets of
different families are disjoint, except that three accessions are shared within Sulfolobaceae.
`GCF_900079115.1` is ND for 273057 and C for the other two. Model fitting (frozen 10406) used only the 11
train targets with R0 features. No test-family genome (target or peer) entered fitting or calibration.

| Target taxid | Assembly | Genus / family | 10406 split | WP4 role | C | F | ND |
|---|---|---|---|---|---:|---:|---:|
| 192952 | GCA_000007065.1 | Methanosarcina / Methanosarcinaceae | test | holdout: model + model-free | 6 | 0 | 1 |
| 272562 | GCA_000008765.1 | Clostridium / Clostridiaceae | test | holdout: model + model-free | 6 | 0 | 1 |
| 198628 | GCA_000147055.1 | Dickeya / Pectobacteriaceae | test | holdout: model + model-free | 6 | 0 | 1 |
| 342108 | GCA_000009985.1 | Paramagnetospirillum / Magnetospirillaceae | test | no-peer control (R0 only) | 0 | 0 | 0 |
| 71421 | GCA_000027305.1 | Haemophilus / Pasteurellaceae | calibration | calibration-used holdout: model-free | 6 | 0 | 1 |
| 192222 | GCA_000009085.1 | Campylobacter / Campylobacteraceae | calibration | calibration-used holdout: model-free | 6 | 0 | 1 |
| 79929 | GCA_000145295.1 | Methanothermobacter / Methanobacteriaceae | calibration | calibration-used holdout: model-free | 6 | 0 | 1 |
| 273057 | GCA_000007005.1 | Saccharolobus / Sulfolobaceae | train | model-free only (dev-examined) | 5 | 1 | 1 |
| 330779 | GCA_000012285.1 | Sulfolobus / Sulfolobaceae | train | model-free only | 3 | 3 | 1 |
| 273063 | GCA_000011205.1 | Sulfurisphaera / Sulfolobaceae | train | model-free only | 2 | 4 | 0 |
| 69014 | GCA_000009965.1 | Thermococcus / Thermococcaceae | train | model-free only | 6 | 0 | 1 |
| 309800 | GCA_000025685.1 | Haloferax / Haloferacaceae | train | model-free only | 6 | 0 | 1 |
| 1140 | GCA_000012525.1 | Synechococcus / Synechococcaceae | train | model-free only | 6 | 0 | 1 |
| 267377 | GCA_000011585.1 | Methanococcus / Methanococcaceae | train | model-free only | 3 | 2 | 1 |
| 62977 | GCA_000046845.1 | Acinetobacter / Moraxellaceae | train | model-free only | 6 | 0 | 1 |
| 272557 | GCA_000011125.1 | Aeropyrum / Acidilobaceae | train | model-free only | 1 | 5 | 0 |
| 273075 | GCA_000195915.1 | Thermoplasma / Thermoplasmataceae | train | model-free only | 1 | 1 | 0 |
| 272624 | GCA_000008485.1 | Legionella / Legionellaceae | train | model-free only | 6 | 0 | 1 |

Peer accessions are exactly the `peers.tsv` rows with that `target_taxid`. The batch for family block B
contains all 18 cohort genomes plus the unique peers of B. This gives 15 batches (Sulfolobaceae is one batch
of 18 + 14 genomes); Paramagnetospirillum has none.

**Engine metadata and holdout.** Cohort convention, kept unchanged: `clade = ncbi_genus_<taxid>`,
`strainCluster = ncbi_species_<taxid>` from NCBI Taxonomy (the cached peer-selection batches). R0 uses
`--cohort-holdout CLADE` as in 10406. R2 uses `STRAIN_CLUSTER`, because CLADE would drop every C peer from
its target. This removes the ND (same species) and keeps the C peers. Under the default `cladeCap = 1`, the
six C clusters share one genus weight. With `jointMinEffective = 0.5`, a reaction is jointly covered by
congeners only if at least three of the six carry it. These are engine defaults and are not tuned (D5).
Because no two cohort genomes share a genus, R0 under STRAIN_CLUSTER equals R0 under CLADE. That is checked
directly (R0-SC). Peers: `COMPLETE`, contigs and `closed` taken from the NCBI report,
`familySearchComplete = true` only when every peer protein has a family row.

## 3. Identical sequence candidates

- Target claims = the 10406 checkpoint `training/diamond_homology` claim files. They are byte-identical in
  every arm, and their sha256 values are listed in each plan.
- Peer claims: DIAMOND 2.1.11 with the 10406 parameters against the **same** cohort-fence reference
  (`~/.cache/gspa-v3-data/training-cohort/baseline/`, faa sha256 `1d8c4ebb…`, labels `db8a8ac9…`).
  `prepare_role_genome_arm.make_claims` builds the contract, and `contract_differs_from_cohort` must be
  empty. Otherwise stop: a mismatch makes the models abstain (2026-10-29 defect 2). D2 covers the case where
  WP2e replaces the models.
- Families: the 10406 cohort orthogroup rows are kept byte-identical. Peer proteins are assigned to the
  existing cluster representatives with MMseqs2 (30 % identity / 80 % coverage, as for the cohort), and the
  rest are clustered among themselves. This assignment differs from joint re-clustering, which is disclosed.
- Checks, for target genomes: identical `candidate_key` sets in `candidate_features.jsonl` across all
  arms; R2 SEQONLY outputs byte-identical to R0 SEQONLY. In `command.json` the only allowed differences
  from the 10406 template are manifest, metadata, orthogroups, pathways, holdout, learned-model and output
  paths. JAR 10402, as for the 10406 models.

## 4. Arms and controls

| Block | Definition | Families | Contexts × kinds |
|---|---|---|---|
| R0 | 18 cohort genomes, CLADE holdout; plus R0-SC (STRAIN_CLUSTER) check | all targets | SEQONLY + 3 relations; none/logistic/gbdt |
| R2 | batch B, STRAIN_CLUSTER holdout | 15 batches (none); 3 test (learned) | 3 relations |
| Add-one | 18 cohort + one C peer j (G0+j) | 6 test/calibration | 3 relations, none |
| R1 check | batch B with CLADE holdout: target features must equal R0 | 3 test | UNIVERSAL |
| S1 orthogroup shuffle | in **peer genomes only**, permute family assignments among proteins of one stratum (contig, gene count within ±`halfWidthBp`, contig-end state). Family-size multiset and window density are preserved. Claims stay on their proteins. 5 seeds | 6 | 3 relations, none |
| S2 pathway shuffle | degree-preserving swap of `reaction_id` across `pathway-groups.tsv` groups (group sizes and per-reaction group counts kept). 3 seeds | 6 | PATHWAY |
| Duplicates | (i) exact copy of top C peer, same cluster; (ii) same copy as a new species cluster; (iii) ND relabelled as its own cluster; (iv) exact copy of the target relabelled as a congeneric species; (v) R2 without ND | 6 | UNIVERSAL none; (i) also logistic |
| True removal | one copy of the target per removed protein: protein, claims, layout row and family row removed; gapsmith re-run on the reduced proteome. Removed: every strict-EC-labelled protein, plus an equal number of random unlabelled EC-claimed proteins per genome (seed 42) as the base rate. All copies share the target's cluster, so holdout keeps them out of each other's peers | 3 test (R0, R2) + 342108 (R0) | 3 relations; none + logistic |
| Fragmentation | target genome cut at contig lengths drawn i.i.d. from fixed log-normals (median 100 kb / 25 kb / 5 kb, σ = 1, min 500 bp, seeds 1–2, the same draws for every target). CDSs crossing a cut are removed together with their claims/memberships (`prepare_sequence_loss` rule); `DRAFT`, `closed = false`; gapsmith re-run. Copies share the target's cluster | as removal | as removal |

The copy layout is checked: in removal/fragmentation batches, the intact target's features must equal its
features in the plain batch.

## 5. "Complementary pieces narrow the same homolog's function": operational definition

For target protein p in orthogroup f, E(p) is the set of complete ECs in p's raw claims; it is the same in
every arm. For genome set G and relation m:

- **Per-homolog, per-function quantity** Q_G(p,e) = `cohort.w.joint_coverage` of candidate (p, e). This is
  the share of e's signature reactions carried next to members of f by peer units of total weight ≥ 0.5.
  When the value is null, the candidate is *unavailable*, never 0.
- Share s_G(p,e) = Q_G(p,e) / Σ Q_G(p,e') over the available e'. It is undefined when the sum is 0.
  Concentration H_G(p) = normalized entropy of s_G when |available| ≥ 2.
- Genome sets: G0 = R0 peers; G0+j = add-one (**before/after adding each genome j**: Δ_j Q(p,e) =
  Q_{G0+j} − Q_{G0}); Gfull = R2 peers.
- **Narrowing** of p: H_Gfull − H_G0 < −δ. For a labelled function F: **toward** if
  s_Gfull(p,F) − s_G0(p,F) > δ, **away** if < −δ, otherwise unchanged. δ = 0.05, fixed now.
- **Complementary narrowing**: toward F **and** `cohort.w.complementarity`(p,F) (joint − mean single
  coverage) is higher under Gfull than under G0 by more than δ. Different relatives then carry different
  signature pieces, not the same pieces again.
- Secondary: rank of F within E(p) by `cohort.w.support_fraction` (ties pessimistic) in G0 versus Gfull.
- The same quantities are computed label-free for **all** multi-candidate proteins of each target. Real
  values are compared with S1/S2 seeds.

## 6. Leakage controls (no holdout annotation builds a network, group or candidate)

| Input | Risk | Mask | Verification |
|---|---|---|---|
| DIAMOND reference (claims) | target self-hits | cohort fence (cohort taxids + genera + Price genera; all test/calibration peers are congeners, so they are fenced as well) | exact-sequence join: 0 reference sequences identical to any target protein. Count target labelled proteins with ≥95 % identity / ≥80 % coverage hits (other genera, disclosed). Count F-tier peer proteins with such hits (better peer anchors, disclosed) |
| gapsmith reference → pathway predictions, gaps, `plausibleReactions` | target (and congener-peer) hits to the target's own characterized Swiss-Prot sequences | new fenced reference: drop every sequence whose UniProt accession falls in the claim fence and every sequence with ≥95 % / ≥80 % to **any** protein of the 18 targets (label-blind). **All** 123 genomes and all perturbation proteomes are searched with it | 0 reference sequences ≥95 % to target proteins; hash list; report how R0 PATHWAY/ORGANISM inputs differ from the 10406 inputs (D4) |
| Orthogroups | none (sequence only) | — | receipt inputs = FASTA + 10406 families only |
| Pathway groups | DB curation only | — | union check of `prepare_context_modes`; diff against 10406 groups |
| Peer proteomes / layouts | RefSeq product names | FASTA headers rewritten to bare ids before any tool; 5-column layouts | header regex; layout header equality |
| Models / thresholds | test labels | frozen 10406 models and display thresholds; no WP4 re-selection | model sha256 equal to 10406 checkpoint |
| Analysis | label-shaped choices | `analyze_wp4_recurrence.py seal` writes the label-blind table (hashed) before `score` joins labels | `score` refuses a table whose hash differs from the seal |

`audit_wp4_leakage.py` fails the run if the declared input of any preparation receipt is a label file
(`labels.tsv`, `swissprot_matches.tsv`, `uniprot_sprot.dat*`, `label_policy.json`, `role_truth*`). The
documented exception is `prepare_wp4_perturbations.py removal`, which reads labels only to choose the
proteins it removes. Static check: no WP4 preparation module imports `match_swissprot_labels`.

## 7. Labels

Strict Swiss-Prot experimental labels only: `training-cohort/labels/labels.tsv` under the experimental
policy (`label_policy.json` `10f3f664…`): complete four-level EC with ECO:0000269 (primary) and
experimental GO (secondary, learned arms only; is_a/part_of propagation with the pinned go-basic). The
curated-inferential tier is excluded. EC identity follows the 10406 `EC:`-prefix normalization; an ENZYME
transfer is used only when it is unambiguous. **Unconfirmed is never negative:** competitor ECs on
labelled proteins and calls on unlabelled proteins are counted as "unconfirmed", never as errors. Without
negatives, no precision, Brier or ECE is reported.

## 8. Metrics (each table partitions its complete denominator)

| Table | Denominator | Categories / statistic |
|---|---|---|
| T1 ledger | all strict EC assertions per target (test: 39 on 37 proteins) | no candidate · F candidate only · multi-candidate with F unavailable G0/Gfull · f present in ≥1 / ≥3 C peers |
| T2 narrowing | same, per relation, pooled over 17 targets and by split | toward · away · unchanged · F unavailable · no alternative · no candidate. One-sided exact sign test (toward vs away), Holm over 3 relations; family-block bootstrap (15 blocks) for mean Δs |
| T3 complementarity | toward events | complementary vs not; real vs each S1/S2 seed |
| T4 label-free | all multi-candidate proteins of the targets | ΔH distribution, share narrowed; real vs shuffles |
| T5 add-one | (assertion, peer j) pairs, 6 families | Δ_j Q(p,F) and Δ_j Q(p,best competitor), by peer tier |
| T6 learned (test) | 39 EC assertions + all target candidates | recovered at frozen threshold, strict top-1, MRR; paired gains/harms R2 vs R0 and vs SEQONLY (exact McNemar); scored / abstained (`unchanged_unsupported`) counts; unconfirmed calls counted |
| T7 controls | target candidates | max abs Δ of features and probabilities, changed-candidate counts for R0-SC, R1, dup (i)/(v), SEQONLY parity (expected 0); inflation for (ii)–(iv) |
| T8 removal | all removals (labelled, random) | F's reaction in observed gaps / `cohort.w.complementary_unannotated` (labelled rate vs random base rate); F reassigned to another protein above threshold (count; unconfirmed) |
| T9 fragmentation | assertions × tier × seed | sequence lost (kept in denominator) · window FULL→TRUNCATED/UNKNOWN · Δ narrowing / Δ recovery |

Power gate (fixed now): fewer than 20 pooled labelled multi-candidate assertions with F available in G0
and Gfull ⇒ T2/T3 are reported as **not testable**, not as a positive or a negative result. Test-only T6 is
expected to be underpowered; it is reported per genome and supports no positive claim on its own.

## 9. Compute budget on unimatrix01

Known costs: gapsmith search ≈ 8 min/genome (≤ 8 G, peak RSS 209 MB); 21-genome integration 5–13 min per
arm on 4 CPUs (6.2–7.8 G). Larger batches are extrapolated linearly in genome count. The pilot measures
the real cost.

| Work | Units | Est. per unit | Task-hours |
|---|---|---|---:|
| Fenced gapsmith searches: 123 base + 74 removal + 24 fragmentation | 221 | 8 min | 29.5 |
| DIAMOND peers (4 threads), measured in pilot | 105 | ≤ 10 min | ≤ 17.5 |
| R0 (5 unscored + 8 learned) | 13 | 8–13 min | 2.4 |
| R2 unscored 45 + SEQONLY parity 3 + R1 3 | 51 | 15 min | 12.8 |
| R2 learned (3 test × 3 × 2) | 18 | 20 min | 6.0 |
| Add-one (6 × 6 × 3) | 108 | 10 min | 18.0 |
| S1 (5 × 6 × 3) + S2 (3 × 6) | 108 | 15 min | 27.0 |
| Duplicates (5 × 6 + 6 logistic) | 36 | 15 min | 9.0 |
| Removal (≈45-genome batches, 24 G) | 42 | 30 min | 21.0 |
| Fragmentation (≈31-genome batches) | 42 | 20 min | 14.0 |
| **Total** | 418 integrations, 221 gapsmith, 105 DIAMOND | | **≈ 157** |

With 8 concurrent searches and 2–4 concurrent 4-CPU integrations (Price tasks leave about 16 G per node;
nodes 003/007 avoided), this is about 2–3 days of wall time. Storage is about 170 MB per run (10406 ratio),
about 70 GB under `/data/hohndor/gspa-v3/wp4-recurrence-v1/`. Only target-genome features, learned JSON,
`cohort.json` and receipts are retrieved. Pre-registered reduction, applied in this order if the pilot run
time exceeds twice the estimate: S1 3 seeds → add-one on the 3 test families only → removal/fragmentation
unscored only.

## 10. Scripts, tests and order of work

| Script (`benchmark/v3/`) | Inputs → outputs | Tests (`benchmark/v3/tests/`) |
|---|---|---|
| `acquire_peer_genomes.py` | `peers.tsv`, `peer_selection.json`, NCBI cache → per-accession protein FASTA (bare ids), GFF-derived 5-column layout, contigs, clade/strain from taxonomy, `provenance.json` with sha256; reuses `acquire_training_cohort` | `test_acquire_peer_genomes.py`: synthetic package, header stripping, twin exclusion, layout header, unplaced = FASTA − layout |
| `prepare_wp4_sequence_inputs.py claims\|families` | peer FASTA + cohort baseline reference → contracted claims (contract gate), MMseqs2 assignment → orthogroups (cohort rows byte-identical) | `test_wp4_sequence_inputs.py`: contract mismatch refused, cohort rows unchanged, every peer protein assigned exactly once |
| `build_fenced_gapsmith_reference.py` | gapsmith reference + claim fence + all target proteins → fenced reference, removal list, hashes (no label input) | `test_fenced_gapsmith_reference.py`: accession and identity rules on synthetic FASTA |
| `prepare_wp4_perturbations.py removal\|fragment\|shuffle\|duplicate` | genomes, claims, families, layouts (+ labels for `removal` only) → reduced/fragmented proteomes, claims, layouts, memberships, metadata; S1/S2 files; duplicate genomes | `test_wp4_perturbations.py`: shuffle preserves per-genome family multiset and strata, never touches targets; fragments use the fixed draws; boundary CDSs dropped everywhere; copies share cluster/clade |
| `observe_gapsmith_cohort.py` (existing) | new manifest naming the fenced reference; scoped flag only if the reference is hard-coded | existing tests + flag test |
| `prepare_wp4_batches.py` | fenced searches → `prepare_context_modes.prepare` (metadata, pathway groups), `prepare_role_genome_context.build_gap_records` (gaps), 10406 `diamond_homology/<CTX>/command.json` → `plans/<block>/<arm>/{manifest.tsv,command.json}`, `arms.tsv`, `plan_report.json` | `test_wp4_batches.py`: every `genome_layout_path` is a layout with the 5-column header (defect 1), `unplacedProteins` correct (defect 2), target claim sha256 identical across arms, holdout flag per regime, SEQONLY has no context flags, allowed command diff only |
| `audit_wp4_leakage.py` | all receipts and references → `leakage.json` (gate) | `test_wp4_leakage.py`: planted label file in a receipt fails; identical-sequence detection |
| `analyze_wp4_recurrence.py seal\|score\|report` | arm outputs (+ labels at `score`) → sealed table, T1–T9 JSON, `docs/benchmarks/v3-wp4-recurrence.md` tables | `test_wp4_analysis.py`: category sums equal denominators, null ≠ 0, pessimistic ties, unconfirmed never FP, seal hash enforced, δ/power gate |

Order, with gates:
1. **Pilot:** 3 Methanosarcina peers → claims (contract gate) → one R2 UNIVERSAL logistic run plus R0-SC.
   Stop and report if the models abstain under STRAIN_CLUSTER or with peers, or if gapsmith reference headers
   cannot be mapped to accessions (then use the identity-only fence). Measure time and memory.
2. Acquire 105 peers → claims + families → fenced reference → perturbation proteomes → search array.
3. Batches + plans → leakage audit → freeze source/plan hashes and `analysis-plan.json` (§5, §8, δ, power
   gate). Model-based arms wait for the verified 10406 `diamond_homology` replay (Task G/G2 resume).
4. Runs in waves: R0, parity and invariance checks first (a failure stops the next wave) → R2/learned →
   add-one → shuffles/duplicates → removal/fragmentation.
5. `seal` → `score` → report; archive `benchmark/v3/results/wp4-recurrence-v1/` with hashes. The supervisor
   updates row 6.

## 11. What counts as a negative result (reported plainly)

- T2: toward ≤ away, or the sign test is not significant after Holm, in a relation → no narrowing toward
  labelled function from close relatives in that relation.
- T2/T3: real toward counts or complementary events do not exceed the maximum over S1 (or, for PATHWAY,
  S2) seeds → the signal cannot be told apart from density-matched background.
- T3: zero complementary narrowing events → no evidence that complementary pieces combine.
- T4: label-free narrowing no more frequent than under shuffles.
- T6: gains − harms ≤ 0 for R2 versus R0 (or versus SEQONLY) at frozen thresholds; or the frozen models
  abstain on R2 features → close-relative recurrence does not improve, or cannot be used by, the current
  models.
- Relations do not differ on paired assertions → no relation advantage.
- T7: any target feature or probability changes beyond 1e-12 under R0-SC, R1, duplicate (i) or (v), or in
  SEQONLY parity → engine/holdout defect (stop and fix before interpreting anything).
- T8: the labelled-removal gap rate is not above the random-removal base rate → recurrence does not
  localize genuinely missing functions; learned reassignment of F more often than SEQONLY → harm.
- T9: R2 gains vanish at the 25 kb tier → the benefit does not hold on draft assemblies.

## 12. Decisions for the supervisor

- **D1:** whether "joint-recurrence arm" means the R2 − R0 contrast inside each relation (this design) or
  a separate fourth arm.
- **D2:** if WP2e retrains `diamond_homology` on the role fence, the peer claims must use that reference
  instead, and WP4 waits for those models.
- **D3:** refitting models on R2 features of the train targets (same splits) is excluded from the budget
  (about +1 day).
- **D4:** the fenced gapsmith reference makes R0 PATHWAY/ORGANISM/gap inputs differ from 10406; the
  frozen models see shifted inputs.
- **D5:** genus-level clade weighting dilutes congeners (at least 3 of 6 needed). Keep the engine
  defaults; no tuning on test.

## 13. Supervisor decisions (2026-10-29)

- **D1:** the joint-recurrence arm is the R2 − R0 contrast inside each relation, as designed.
- **D2:** model-free arms (unscored context) start now with the cohort-fence reference; the learned arms wait for
  the WP2e outcome and use whichever reference the models they score were trained with (one contract per run).
- **D3:** no refitting on R2 features in this work package.
- **D4:** accepted; R0 is re-run under the fenced gapsmith reference so that R2 − R0 is internally consistent, and
  the difference from the 10406 inputs is reported, not interpreted.
- **D5:** engine defaults kept; no tuning on test.
- Order confirmed: pilot (§10 step 1) first; the power gate and the pre-registered reductions stand.

## 14. Pilot outcome (step 1 of §10, Shaheen job 15691175, 2026-10-29)

Methanosarcina block, three close peers, UNIVERSAL × {none, logistic}, regimes R0 / R0-SC / R1 / R2 (plans with
per-arm organism metadata restricted to each manifest; the engine rejects metadata naming genomes outside the
manifest). All eight arms ran (10–12 min each, 21 genomes, 4 threads, ≤17 GB). On the target's 27,062 candidate
rows: the unscored arms are identical across all four regimes; logistic R1 equals R0 exactly (the same-genus peers
are held out under CLADE, as designed); R0-SC differs from R0 in 194 rows, none by more than 0.05; R2 differs from R0
in 1,000 rows, 257 by more than 0.05, so the R2 − R0 contrast of D1 is non-empty; the frozen 10406 logistic model
rescored 19,137 target rows and 21,061 rows of a peer, i.e. the models do not abstain with peers or under
STRAIN_CLUSTER. The gates of step 1 therefore pass; steps 2–5 proceed with the fenced gapsmith reference.
