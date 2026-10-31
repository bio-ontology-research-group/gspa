# Joint recurrence from close relatives (WP4)

This benchmark asks whether adding unlabelled close relatives to the joint
cohort helps assign function to the same sequence candidates. It tests whether
the context evidence of a homolog then concentrates on its experimentally
labelled function, and whether different relatives contribute complementary
pathway pieces. It also asks whether frozen learned models recover more labelled
functions from the resulting features. The 18 cohort genomes are 18 genera, so
every earlier recurrence feature came from distant genomes. WP4 is the first run
with congeners in the cohort.

The result is negative. Too few labelled assertions keep their function
available with and without the relatives, so the narrowing test (T2) fails the
pre-registered power gate in all three context relations. It is therefore not
testable. Among the assertions that are available, none moves toward or away
from the labelled function, and no complementary narrowing occurs. Label-free
narrowing is rare and exceeds the shuffles only in the universal-network
relation, on very few proteins. The frozen models gain no recovery from the
close-relative features. Removing a labelled gene does not put its function
into the observed gaps more often than removing a random gene. Every invariance
control passes on the corrected engine. The first full run (v1) failed one of
these controls; it was invalidated, the engine was fixed and the run was
repeated in full.

Design: `tasks/v3-wp4-recurrence-design-2026-10-29.md` (§1–§14; §13 records the
supervisor decisions D1–D5, §14 the pilot). Run log:
`tasks/v3-supervision-plan-2026-10-28.md`, entries 2026-10-29 to 2026-10-30.
All tables below are copied or condensed from the generated report
`benchmark/v3/results/wp4-recurrence-v2/report/v3-wp4-recurrence-tables.md` and
from `scores/metrics.json` in the same directory. Claims are limited to these
genomes, this engine configuration and DIAMOND/Swiss-Prot candidates. A
concentration of context evidence does not show that a function is correct.

## Outcome against the pre-registered criteria

The criteria are those of design §11.

| Test | Negative result if (design §11) | Outcome |
|---|---|---|
| T2 narrowing toward the labelled function | toward ≤ away, or not significant after Holm; or not above shuffles | **Not testable** in all three relations (power gate); no toward or away event |
| T3 complementary narrowing | zero complementary events, or not above shuffles | No complementary event in any relation |
| T4 label-free narrowing | not more frequent than under shuffles | Negative for PATHWAY and ORGANISM; UNIVERSAL above every S1 seed on very few proteins (see [below](#what-the-t4-universal-observation-means)) |
| T5 add-one congener | (descriptive) | No pair in which the labelled function gains above its best competitor |
| T6 frozen learned models | gains − harms ≤ 0 for R2 vs R0, or abstention on R2 | No gain in any relation or model; models do not abstain on R2 |
| Relation comparison | relations do not differ on paired assertions | No relation advantage |
| T7 invariance controls | any change beyond 1e-12 in a zero control | No defect in v2 (v1 failed R0-SC; see [v1 defect](#first-run-v1-invalid-and-the-seed-fix)) |
| T8 true removal | labelled gap rate not above the random base rate | Negative in every regime, relation and model; harm versus SEQONLY not evaluable |
| T9 fragmentation | R2 gains vanish at the 25 kb tier | No R2 narrowing exists to lose; recovery mostly drops with fragmentation |

## Design

### Targets, splits and peers

The holdout unit is the NCBI family block: a target plus all its peer rows in
`benchmark/v3/results/peer-selection-v1/peers.tsv` (111 rows, 105 unique
accessions). Tiers: C = same genus, F = same family, other genus, ND = near
duplicate (same species). The frozen 10406 `diamond_homology` models were fitted
on the eleven train targets with distant-cohort features only, so no test-family
genome entered fitting or calibration.

| Target taxid | Assembly | Genus / family | 10406 split | WP4 role | C | F | ND |
|---|---|---|---|---|---:|---:|---:|
| 192952 | GCA_000007065.1 | Methanosarcina / Methanosarcinaceae | test | holdout: model + model-free | 6 | 0 | 1 |
| 272562 | GCA_000008765.1 | Clostridium / Clostridiaceae | test | holdout: model + model-free | 6 | 0 | 1 |
| 198628 | GCA_000147055.1 | Dickeya / Pectobacteriaceae | test | holdout: model + model-free | 6 | 0 | 1 |
| 342108 | GCA_000009985.1 | Paramagnetospirillum / Magnetospirillaceae | test | no-peer control (R0 only) | 0 | 0 | 0 |
| 71421 | GCA_000027305.1 | Haemophilus / Pasteurellaceae | calibration | model-free | 6¹ | 0¹ | 1 |
| 192222 | GCA_000009085.1 | Campylobacter / Campylobacteraceae | calibration | model-free | 6 | 0 | 1 |
| 79929 | GCA_000145295.1 | Methanothermobacter / Methanobacteriaceae | calibration | model-free | 6 | 0 | 1 |
| 273057 | GCA_000007005.1 | Saccharolobus / Sulfolobaceae | train | model-free (dev-examined) | 5 | 1 | 1 |
| 330779 | GCA_000012285.1 | Sulfolobus / Sulfolobaceae | train | model-free | 3 | 3 | 1 |
| 273063 | GCA_000011205.1 | Sulfurisphaera / Sulfolobaceae | train | model-free | 2 | 4 | 0 |
| 69014 | GCA_000009965.1 | Thermococcus / Thermococcaceae | train | model-free | 6 | 0 | 1 |
| 309800 | GCA_000025685.1 | Haloferax / Haloferacaceae | train | model-free | 6 | 0 | 1 |
| 1140 | GCA_000012525.1 | Synechococcus / Synechococcaceae | train | model-free | 6 | 0 | 1 |
| 267377 | GCA_000011585.1 | Methanococcus / Methanococcaceae | train | model-free | 3 | 2 | 1 |
| 62977 | GCA_000046845.1 | Acinetobacter / Moraxellaceae | train | model-free | 6 | 0 | 1 |
| 272557 | GCA_000011125.1 | Aeropyrum / Acidilobaceae | train | model-free | 1 | 5 | 0 |
| 273075 | GCA_000195915.1 | Thermoplasma / Thermoplasmataceae | train | model-free | 1 | 1 | 0 |
| 272624 | GCA_000008485.1 | Legionella / Legionellaceae | train | model-free | 6 | 0 | 1 |

¹ Design §2 values. At planning, one Pasteurellaceae peer (`GCF_031191205.1`) is
placed in *Parahaemophilus* by current NCBI taxonomy, i.e. tier F. The add-one
block therefore has one peer fewer for this family (see [Execution](#execution)).

Sulfolobaceae share three peer accessions and form one batch, so there are 15
peer batches; Paramagnetospirillum has no peers.

### Regimes, relations and model kinds

| Regime | Genomes in the batch | Cohort holdout | Role |
|---|---|---|---|
| R0 | the 18 cohort genomes | `CLADE` | baseline (10406 setting, fenced gapsmith reference, D4) |
| R0-SC | the 18 cohort genomes | `STRAIN_CLUSTER` | invariance control: no two cohort genomes share a genus, so it must equal R0 |
| R1 | cohort + the family's peers | `CLADE` | invariance control: congeners are held out, so target features must equal R0 (3 test families, UNIVERSAL) |
| R2 | cohort + the family's peers | `STRAIN_CLUSTER` | joint-recurrence arm: the near duplicate is held out, congeners stay |
| ref10406 | the 18 cohort genomes, 10406 inputs | `CLADE` | report-only comparison of fenced versus original inputs (D4) |

The joint-recurrence effect is the R2 − R0 contrast inside each relation
(decision D1). The relations are the context arms PATHWAY (pathway-only
groups), UNIVERSAL (universal reaction network) and ORGANISM (organism model);
SEQONLY is the candidate anchor without context. Each arm runs as `none`
(unscored context features) and, for the learned tables, with the frozen 10406
`logistic` and `gbdt` models (D3: no refitting). Engine defaults are unchanged
(D5): `clade = ncbi_genus_<taxid>`, `strainCluster = ncbi_species_<taxid>`,
`cladeCap = 1`, `jointMinEffective = 0.5`, so congeners share one genus weight.

Target claims are the byte-identical 10406 `diamond_homology` claim files in
every arm. Peer claims come from DIAMOND against the same cohort-fence reference
under the same predictor contract. The cohort orthogroup rows are kept
byte-identical and peer proteins are assigned to them with MMseqs2.

### Perturbations and controls

| Block | Definition (design §4) | Families |
|---|---|---|
| Add-one | 18 cohort genomes + one C peer | 6 test/calibration |
| S1 orthogroup shuffle | family labels permuted within matched strata, in peer genomes only; 5 seeds | 6 test/calibration |
| S2 pathway shuffle | degree-preserving swap of reactions between pathway groups; 3 seeds | 6, PATHWAY only |
| Duplicates (i)–(v) | (i) exact copy of the top C peer in its cluster; (ii) the same copy as a new species; (iii) the ND relabelled as its own species; (iv) an exact copy of the target relabelled as a congeneric species; (v) R2 without the ND | 6 |
| True removal | one target copy per removed protein (protein, claims, layout row, family row removed; gapsmith re-run). Removed: every strict-EC-labelled protein, plus as base rate an equal number of random unlabelled proteins with a complete EC claim (seed 42) | 3 test (R0, R2) + 342108 (R0) |
| Fragmentation | target cut at i.i.d. log-normal contig lengths (medians 100 / 25 / 5 kb, σ = 1, seeds 1–2), CDSs crossing a cut removed with their claims and memberships; gapsmith re-run | as removal |

### Quantities and gates

For target protein p, the per-function quantity Q(p,e) is
`cohort.w.joint_coverage`. A null value means *unavailable*, never 0. The share
s(p,e) normalizes Q over the available candidate functions, and the
concentration H(p) is the normalized entropy of s. G0 is the R0 peer set and
Gfull the R2 peer set. An assertion moves *toward* its labelled function F if
s(p,F) rises by more than δ from G0 to Gfull. It moves *away* if s(p,F) falls by
more than δ. It is *complementary* if it moves toward F and
`cohort.w.complementarity` also rises by more than δ. A protein *narrows* if H
falls by more than δ. The following values are frozen in `analysis_plan` of
`seal/seal.json`, identical in v1 and v2.

| Parameter | Value |
|---|---|
| δ (toward / away / narrowing) | 0.05 |
| Power gate (T2/T3): labelled multi-candidate assertions with F available in G0 and Gfull, pooled | ≥ 20, otherwise *not testable* |
| Significance | one-sided exact sign test, α = 0.05, Holm over the 3 relations; exact McNemar for paired contrasts |
| Bootstrap of mean Δs | family-block, 2,000 replicates, seed 20261029 |
| Invariance tolerance (T7) | 1e-12 |
| Draft tier (T9) | 25 kb |
| Labels | strict Swiss-Prot experimental EC (ECO:0000269), four-level; experimental GO secondary (learned arms only) |
| Unconfirmed calls | counted, never false positives; no precision, Brier or ECE |

T6 on the test split was expected to be underpowered. It is reported per genome
and supports no positive claim on its own (design §8).

### Leakage controls

`audit_wp4_leakage.py` ran per domain on the re-planned v2 inputs and gates the
runs. The analysis seals a label-blind table (hashed) before `score` joins
labels; `score` refuses a table whose hash differs from the seal.

| Check (`leakage-archaea.json` / `leakage-bacteria.json`) | Archaea | Bacteria |
|---|---:|---:|
| Status | passed | passed |
| Preparation receipts scanned for label inputs | 1,069 | 1,069 |
| Receipts with a label-file hit | 0 | 0 |
| Fenced gapsmith reference sequence files | 46,464 | 41,265 |
| Reference sequences identical to any target protein | 0 | 0 |
| Distinct target sequences checked | 49,204 | 49,204 |
| Accessions on the reference removal list | 6,514 | 36,375 |
| Static import check (`match_swissprot_labels`), modules clean | 9 of 9 | 9 of 9 |

## First run (v1), invalid, and the seed fix

The first full analysis (Shaheen job 15692636) failed the R0-SC control.
STRAIN_CLUSTER and CLADE holdouts select identical peer pools in R0, yet
candidate feature records differed. The sealed numeric columns were identical
and only `features_sha256` differed. In the diagnosed case (Saccharolobus,
ORGANISM, run log 2026-10-30 02:40) the record diff was in
`cohort.w.perm_excess`. `CohortInference` derived the seed of the permutation
null from family, reaction and the holdout key (`c:<clade>` or `s:<cluster>`).
The same peers therefore received different Monte Carlo draws under the two
holdouts. Because R2 uses STRAIN_CLUSTER and R0 uses CLADE, this noise would also
enter the R2 − R0 contrast through the learned permutation features. Under
design §11 a failed invariance control stops interpretation. The v1 results are
kept only as a record of the failure (`benchmark/v3/results/wp4-recurrence-v1/`).

The engine now seeds the permutation from content (`permutationSeed`
`content:v2`). The seed is a SHA-256 over the scheme, the configured seed, the
family, the reaction and a canonical signature of the peer evidence units and
their background pools; see
[the permutation note](v3-cohort-permutation.md#holdout-independent-permutation-seed-contentv2).
`CohortPermutationSeedSpec` adds seven cases that fail on the old code. The JAR
bundle `ec-prefix-identity-v4` was built from the 10402 snapshot plus this patch
(unimatrix01 job 10699). Every WP4 run was re-planned onto it and rerun, then
sealed, scored and reported with the unchanged analyzer and analysis plan.

| Control R0-SC versus R0 (`metrics.json`, T7) | Rows (target × relation) | Rows with changed candidates | Changed candidates per row | Max abs Δ, sealed numerics |
|---|---:|---:|---|---:|
| v1, hold-out-keyed seed | 54 | 54 | 117–561 | 0 |
| v2, `content:v2` seed | 54 | 0 | 0 | 0 |

The two seals carry the same analyzer SHA-256 and the same analysis plan.
Between the v1 and v2 generated tables, only two things differ: the R0-SC rows
and a single unconfirmed-call count of one R0 learned arm (Methanosarcina,
PATHWAY, logistic). The new seed changes the permutation draws in every regime,
not only in R0-SC. In these tables that change is visible nowhere else.

## Results (v2)

### T1: ledger of labelled assertions

All strict EC assertions per target. "No candidate": the labelled function is
not among the protein's DIAMOND candidates. "F candidate only": it is the only
candidate. "Multi-candidate": competitors exist. The C-peer family columns are
NA when the C-peer count is unknown for at least one assertion protein of that
target.

| Assembly | Split | Assertions | Proteins | No candidate | F candidate only | Multi-candidate | Family in ≥1 C peer | Family in ≥3 C peers |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| GCA_000007005.1 | train | 32 | 25 | 11 | 9 | 12 | NA | NA |
| GCA_000007065.1 | test | 21 | 20 | 16 | 4 | 1 | NA | NA |
| GCA_000008485.1 | train | 4 | 4 | 2 | 1 | 1 | NA | NA |
| GCA_000008765.1 | test | 13 | 12 | 7 | 2 | 4 | NA | NA |
| GCA_000009085.1 | calibration | 15 | 14 | 9 | 3 | 3 | NA | NA |
| GCA_000009965.1 | train | 42 | 36 | 12 | 13 | 17 | NA | NA |
| GCA_000009985.1 | test | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| GCA_000011125.1 | train | 20 | 17 | 11 | 2 | 7 | NA | NA |
| GCA_000011205.1 | train | 28 | 19 | 15 | 4 | 9 | 27 | 0 |
| GCA_000011585.1 | train | 15 | 10 | 11 | 2 | 2 | NA | NA |
| GCA_000012285.1 | train | 28 | 26 | 15 | 5 | 8 | NA | NA |
| GCA_000012525.1 | train | 8 | 8 | 2 | 4 | 2 | 6 | 5 |
| GCA_000025685.1 | train | 23 | 22 | 11 | 5 | 7 | NA | NA |
| GCA_000027305.1 | calibration | 21 | 20 | 2 | 9 | 10 | 13 | 12 |
| GCA_000046845.1 | train | 13 | 12 | 2 | 4 | 7 | 13 | 13 |
| GCA_000145295.1 | calibration | 28 | 27 | 21 | 4 | 3 | NA | NA |
| GCA_000147055.1 | test | 5 | 5 | 1 | 3 | 1 | 5 | 5 |
| GCA_000195915.1 | train | 19 | 14 | 10 | 5 | 4 | NA | NA |
| all | | 335 | 291 | 158 | 79 | 98 | NA | NA |

Multi-candidate assertions by relation, summed over the per-target rows of the
generated table:

| Relation | Multi-candidate | F unavailable in G0 | F unavailable in Gfull only | F available in both | Arm missing |
|---|---:|---:|---:|---:|---:|
| PATHWAY | 98 | 87 | 2 | 9 | 0 |
| UNIVERSAL | 98 | 90 | 0 | 8 | 0 |
| ORGANISM | 98 | 93 | 0 | 5 | 0 |

The binding constraint sits before the relatives are added. Most labelled
functions have no candidate at all, and for most multi-candidate assertions the
labelled function has no joint coverage even in the distant cohort.

### T2: narrowing toward the labelled function (R2 − R0)

| Relation | Assertions | Toward | Away | Unchanged | F unavailable | No alternative | No candidate | Arm missing | F available both / gate | Sign test p | Mean Δs (block bootstrap 95%) | Bootstrap blocks | Verdict |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|---:|---|
| PATHWAY | 335 | 0 | 0 | 9 | 89 | 79 | 158 | 0 | 9 / 20 | NA | 0 (0–0) | 6 | not testable |
| UNIVERSAL | 335 | 0 | 0 | 8 | 90 | 79 | 158 | 0 | 8 / 20 | NA | 0 (0–0) | 4 | not testable |
| ORGANISM | 335 | 0 | 0 | 5 | 93 | 79 | 158 | 0 | 5 / 20 | NA | 0 (0–0) | 3 | not testable |

All three relations fail the power gate, so T2 is **not testable**. This is
neither a positive nor a negative finding about narrowing. Descriptively, no
assertion with F available in both genome sets changed its share beyond δ. The
rank of F by support fraction is also unchanged for every such assertion (T2
secondary in the generated tables). The bootstrap resamples only the family
blocks that contribute assertions; with every Δs equal to zero its interval is
degenerate. The split breakdown (`T2 by split`) contains no toward or away
event in any split.

Relation comparison on paired assertions (exact McNemar):

| Relations | Paired | Toward only first | Toward only second | p |
|---|---:|---:|---:|---|
| PATHWAY vs ORGANISM | 3 | 0 | 0 | NA |
| PATHWAY vs UNIVERSAL | 5 | 0 | 0 | NA |
| UNIVERSAL vs ORGANISM | 5 | 0 | 0 | NA |

### T3: complementary narrowing

| Relation | Toward events | Complementary | Not complementary | Real toward / complementary above every seed | Shuffle seeds (all nulls zero) |
|---|---:|---:|---:|---|---|
| PATHWAY | 0 | 0 | 0 | no / no | S1 × 5, S2 × 3 |
| UNIVERSAL | 0 | 0 | 0 | no / no | S1 × 5 |
| ORGANISM | 0 | 0 | 0 | no / no | S1 × 5 |

With no toward event, there is nothing for complementary pieces to explain. The
pre-registered statement applies: no evidence that complementary pieces
combine. It carries the T2 power caveat.

### T4: label-free concentration (all multi-candidate target proteins)

| Relation | Proteins | Narrowed | Broadened | Unchanged | H undefined | Arm missing | Share narrowed | ΔH median (q25–q75) | Narrowed above every seed |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|
| PATHWAY | 7613 | 17 | 20 | 698 | 6239 | 639 | 0.02313 | 0 (0–0) | no |
| UNIVERSAL | 7613 | 3 | 22 | 669 | 6280 | 639 | 0.004323 | 0 (0–0) | yes |
| ORGANISM | 7613 | 6 | 17 | 455 | 6496 | 639 | 0.01255 | 0 (0–0) | no |

"Share narrowed" is relative to proteins with H defined in both genome sets.
The shuffle comparison uses the six families with shuffles (six target genomes):

| Relation | Shuffle | Multi-candidate proteins | Real narrowed | Null narrowed per seed |
|---|---|---:|---:|---|
| PATHWAY | S1 (5 seeds) | 2554 | 1 | 0, 0, 0, 1, 0 |
| PATHWAY | S2 (3 seeds) | 2554 | 1 | 13, 20, 9 |
| UNIVERSAL | S1 (5 seeds) | 2554 | 3 | 0, 0, 0, 0, 0 |
| ORGANISM | S1 (5 seeds) | 2554 | 0 | 0, 0, 0, 0, 0 |

For PATHWAY and ORGANISM, label-free narrowing is no more frequent than under
shuffles. Every PATHWAY S2 seed narrows more proteins than the real pathway
groups. For the UNIVERSAL observation, see [below](#what-the-t4-universal-observation-means).

### T5: adding one congener at a time

| Relation | Tier | Pairs (assertion × peer) | F gains above best competitor | Not above (ties incl.) | F unavailable | Competitor unavailable | No alternative | No candidate | Mean ΔQ(F) | Mean ΔQ(best competitor) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| PATHWAY | C | 597 | 0 | 43 | 62 | 17 | 141 | 334 | 0 | 0 |
| UNIVERSAL | C | 597 | 0 | 61 | 50 | 11 | 141 | 334 | 0 | 0 |
| ORGANISM | C | 597 | 0 | 32 | 62 | 28 | 141 | 334 | 0 | 0 |

No single congener raises the joint coverage of a labelled function above its
best competitor. This is consistent with the genus-level weighting of D5, under
which one congener alone carries limited joint weight. The design did not test
that mechanism directly.

### T6: frozen learned models on the test targets

Recovery at the frozen display threshold over the strict EC assertions of the
test targets with peers, summed over the per-target rows of the generated
table. Paramagnetospirillum has no EC assertion and no R2 arm. The per-target
rows, with candidate and unconfirmed-call counts, are in the generated tables.

| Relation | Kind | R0 recovered | R2 recovered | R0 strict top-1 | R2 strict top-1 | Below threshold R0 / R2 | Abstained R0 / R2 | Candidate absent |
|---|---|---:|---:|---:|---:|---|---|---:|
| SEQONLY | logistic | 0 | – | 0 | – | 0 / – | 15 / – | 24 |
| SEQONLY | gbdt | 0 | – | 0 | – | 0 / – | 15 / – | 24 |
| PATHWAY | logistic | 11 | 9 | 14 | 14 | 4 / 6 | 0 / 0 | 24 |
| PATHWAY | gbdt | 12 | 12 | 11 | 12 | 3 / 3 | 0 / 0 | 24 |
| UNIVERSAL | logistic | 11 | 10 | 14 | 14 | 4 / 5 | 0 / 0 | 24 |
| UNIVERSAL | gbdt | 10 | 10 | 13 | 13 | 5 / 5 | 0 / 0 | 24 |
| ORGANISM | logistic | 7 | 3 | 12 | 14 | 7 / 11 | 1 / 1 | 24 |
| ORGANISM | gbdt | 10 | 10 | 11 | 12 | 4 / 4 | 1 / 1 | 24 |

Paired gains and harms, pooled over the test targets (exact McNemar):

| Contrast | Relation | Kind | EC: assertions | EC: gains | EC: harms | EC: p | GO: assertions | GO: gains | GO: harms | GO: p |
|---|---|---|---:|---:|---:|---|---:|---:|---:|---|
| R2 vs R0 | PATHWAY | logistic | 39 | 0 | 2 | 0.5 | 114 | 0 | 0 | NA |
| R2 vs R0 | PATHWAY | gbdt | 39 | 0 | 0 | NA | 114 | 0 | 0 | NA |
| R2 vs R0 | UNIVERSAL | logistic | 39 | 0 | 1 | 1 | 114 | 0 | 0 | NA |
| R2 vs R0 | UNIVERSAL | gbdt | 39 | 0 | 0 | NA | 114 | 0 | 0 | NA |
| R2 vs R0 | ORGANISM | logistic | 39 | 0 | 4 | 0.125 | 114 | 0 | 0 | NA |
| R2 vs R0 | ORGANISM | gbdt | 39 | 0 | 0 | NA | 114 | 0 | 0 | NA |
| R2 vs SEQONLY | PATHWAY | logistic | 39 | 9 | 0 | 0.003906 | 114 | 3 | 0 | 0.25 |
| R2 vs SEQONLY | PATHWAY | gbdt | 39 | 12 | 0 | 0.0004883 | 114 | 11 | 0 | 0.0009766 |
| R2 vs SEQONLY | UNIVERSAL | logistic | 39 | 10 | 0 | 0.001953 | 114 | 3 | 0 | 0.25 |
| R2 vs SEQONLY | UNIVERSAL | gbdt | 39 | 10 | 0 | 0.001953 | 114 | 10 | 0 | 0.001953 |
| R2 vs SEQONLY | ORGANISM | logistic | 39 | 3 | 0 | 0.25 | 114 | 2 | 0 | 0.5 |
| R2 vs SEQONLY | ORGANISM | gbdt | 39 | 10 | 0 | 0.001953 | 114 | 10 | 0 | 0.001953 |

Close-relative recurrence does not improve the current models in any relation.
R2 never gains over R0. The logistic models lose recoveries under R2, but no
loss is significant and the test split is small. The models score R2 features;
none of the six relation/model combinations abstains on R2. The R2-versus-SEQONLY
gains do not show a benefit of relatives. The SEQONLY models abstain on every
test candidate (`unchanged_unsupported`), so any context arm that scores at all
gains over them. R0 alone reaches the same or a higher recovery than R2. The same
pattern holds for the secondary GO labels.

### T7: invariance and inflation controls

Condensed from 441 control rows. A row is one target × relation × model kind
(× regime for copy-layout checks).

| Control | Expectation | Targets | Rows | Rows with changed candidates | Changed candidates per row | Max abs Δ | Defect rows |
|---|---|---:|---:|---:|---|---:|---:|
| R0-SC vs R0 | zero | 18 | 54 | 0 | 0 | 0 | 0 |
| R1 vs R0 | zero | 3 | 3 | 0 | 0 | 0 | 0 |
| Duplicate (i): exact copy of a C peer, same cluster (none + logistic) | zero | 6 | 9 | 0 | 0 | 0 | 0 |
| Duplicate (v): R2 without the near duplicate | zero | 6 | 6 | 0 | 0 | 0 | 0 |
| SEQONLY parity, R2 vs R0 | zero | 3 | 3 | 0 | 0 | 0 | 0 |
| Copy layout: intact target in removal/fragment batches vs plain batch | zero | 4 | 294 | 0 | 0 | 0 | 0 |
| Duplicate (ii): copy as a new species | inflation | 6 | 6 | 6 | 520–2013 | 0.2069 | – |
| Duplicate (iii): ND as its own species | inflation | 6 | 6 | 6 | 535–2014 | 0.25 | – |
| Duplicate (iv): target copy as a congeneric species | inflation | 6 | 6 | 6 | 626–2018 | 0.6667 | – |
| Fenced vs 10406 inputs (D4) | report only | 18 | 54 | 54 | 117–422 | 1 | – |

The candidate-key sets match across arms (0 mismatches). Every control that must
be invariant is exactly invariant, so no engine or holdout defect remains under
design §11. Duplicates labelled as independent species change features as
expected. That is the inflation the design asks to report, not a defect.
Replacing the 10406 gapsmith inputs with the fenced reference changes R0
features for every target. Per D4 this is reported, not interpreted.

### T8: true removal

Gap rates are identical across relations and model kinds within each class,
because the observed gaps come from the gapsmith search on the reduced
proteome. Labelled removals count one row per strict EC label of a removed
labelled protein. Random removals count one row per complete EC claim of a
removed random protein (`removal_selection.tsv` convention), hence the different
denominators.

| Class | Removals | F in observed gaps | F not in gaps | Gaps unavailable | Gap rate |
|---|---:|---:|---:|---:|---:|
| labelled | 39 | 2 | 37 | 0 | 0.05128 |
| random (base rate) | 69 | 4 | 65 | 0 | 0.05797 |

| Regime | Relation | Mean complementary unannotated, labelled | Mean complementary unannotated, random | F reassigned above threshold, logistic, labelled (unconfirmed) | Same, random |
|---|---|---:|---:|---:|---:|
| R0 | PATHWAY | NA | 0.01077 | 6 | 16 |
| R0 | UNIVERSAL | 0.02646 | 0.004705 | 6 | 17 |
| R0 | ORGANISM | 0.2083 | 0.0433 | 2 | 14 |
| R2 | PATHWAY | 0 | 0.01801 | 6 | 13 |
| R2 | UNIVERSAL | 0.03968 | 0.009337 | 7 | 14 |
| R2 | ORGANISM | 0.4167 | 0.05812 | 1 | 8 |

The labelled gap rate is not above the random base rate in any of the twelve
regime × relation × kind cells (`T8 labelled rate vs random base rate`).
Recurrence therefore does not localize genuinely missing functions in these
genomes. No removal arm was run without context, so the pre-registered harm
comparison (learned reassignment versus SEQONLY) has no anchor and is not
evaluated. The reassignment counts above are unconfirmed calls, not errors.

### T9: fragmentation

The sequence and window columns compare the R0 features of the intact target
with those of its fragmented copy. The toward columns apply the T2
classification to the intact target and to the copy. The recovery columns
compare R2 logistic recoveries (from `metrics.json`; there are no gbdt fragment
arms). Denominator:
the EC assertions of the test targets.

| Relation | Tier | Seed | Sequence lost | Window FULL → TRUNCATED/UNKNOWN | Window unchanged | Window unobserved in plain | No candidate | Toward plain / fragment | Recovered plain / fragment (R2 logistic) |
|---|---|---:|---:|---:|---:|---:|---:|---|---|
| PATHWAY | 100k | 1 | 0 | 1 | 7 | 7 | 24 | 0 / 0 | 9 / 8 |
| PATHWAY | 100k | 2 | 0 | 0 | 8 | 7 | 24 | 0 / 0 | 9 / 9 |
| PATHWAY | 25k | 1 | 1 | 2 | 6 | 6 | 24 | 0 / 0 | 9 / 5 |
| PATHWAY | 25k | 2 | 0 | 3 | 5 | 7 | 24 | 0 / 0 | 9 / 8 |
| PATHWAY | 5k | 1 | 2 | 2 | 4 | 7 | 24 | 0 / 0 | 9 / 6 |
| PATHWAY | 5k | 2 | 3 | 6 | 0 | 6 | 24 | 0 / 0 | 9 / 5 |
| UNIVERSAL | 100k | 1 | 0 | 1 | 7 | 7 | 24 | 0 / 0 | 10 / 10 |
| UNIVERSAL | 100k | 2 | 0 | 0 | 8 | 7 | 24 | 0 / 0 | 10 / 10 |
| UNIVERSAL | 25k | 1 | 1 | 2 | 6 | 6 | 24 | 0 / 0 | 10 / 7 |
| UNIVERSAL | 25k | 2 | 0 | 2 | 6 | 7 | 24 | 0 / 0 | 10 / 10 |
| UNIVERSAL | 5k | 1 | 2 | 3 | 4 | 6 | 24 | 0 / 0 | 10 / 7 |
| UNIVERSAL | 5k | 2 | 3 | 5 | 1 | 6 | 24 | 0 / 0 | 10 / 5 |
| ORGANISM | 100k | 1 | 0 | 1 | 5 | 9 | 24 | 0 / 0 | 3 / 3 |
| ORGANISM | 100k | 2 | 0 | 0 | 6 | 9 | 24 | 0 / 0 | 3 / 3 |
| ORGANISM | 25k | 1 | 1 | 2 | 4 | 8 | 24 | 0 / 0 | 3 / 1 |
| ORGANISM | 25k | 2 | 0 | 2 | 4 | 9 | 24 | 0 / 0 | 3 / 2 |
| ORGANISM | 5k | 1 | 2 | 3 | 3 | 7 | 24 | 0 / 0 | 3 / 1 |
| ORGANISM | 5k | 2 | 3 | 3 | 1 | 8 | 24 | 0 / 0 | 3 / 0 |

R2 has no narrowing gain, so there is none to vanish on draft assemblies. The
§11 draft criterion has nothing to act on. Fragmentation removes sequences and
truncates windows more often at shorter contig medians. R2 logistic recovery
stays the same or falls under every fragmentation; it never rises. Two seeds per
tier give no estimate of variability.

## Negative results

Stated against design §11 (the analyzer's selection is
`report/v3-wp4-negative-results.md`):

- **Narrowing is not testable** in the PATHWAY, UNIVERSAL and ORGANISM
  relations. Too few labelled multi-candidate assertions keep F available in
  both genome sets for the power gate. This is neither a positive nor a negative
  result, and no descriptive toward or away event occurred.
- **No complementary narrowing**: no event in any relation, so there is no
  evidence that complementary pieces from different relatives combine.
- **Label-free narrowing is no more frequent than under shuffles** for PATHWAY
  and ORGANISM. The PATHWAY S2 shuffles narrow more proteins than the real
  groups.
- **No single-congener effect**: adding one congener never lifts a labelled
  function above its best competitor.
- **No learned-model gain**: R2 does not improve either frozen model over R0 in
  any relation. The logistic models lose recoveries (not significant).
- **No relation advantage** on paired assertions.
- **No localization of missing functions**: removing a labelled gene does not
  put its function into the observed gaps more often than removing a random gene,
  in either regime or any relation.
- **No draft-assembly benefit to lose**: fragmentation lowers or keeps R2
  logistic recovery; there was no narrowing gain to begin with.

## What the T4 UNIVERSAL observation means

The only "yes" in the generated tables is T4 "Narrowed above every seed" for
UNIVERSAL. The rule compares real narrowing with the maximum over seeds. In the
six shuffled families, the real R2 − R0 contrast narrowed the context
concentration of a few multi-candidate target proteins. No S1 seed narrowed any.

What it does mean: on the universal reaction network, adding the real
relatives' family structure narrowed these few proteins. Permuting family
labels among matched windows of the peer genomes removes that narrowing.

What it does not mean:

- **It is not a frequent effect.** It concerns only the proteins counted in the
  T4 table, among thousands of multi-candidate proteins. In the same relation
  more proteins broadened than narrowed, and the median ΔH is zero.
- **It is label-free.** It says nothing about whether the concentration moved
  toward a correct function. On the labelled assertions, no UNIVERSAL
  assertion moved toward or away (T2), and the test is not powered.
- **It is not a statistical test.** T4 has no pre-registered significance level.
  With five seeds whose nulls are all zero, a real count above every seed has a
  rank probability of at best one in six under exchangeability. Only the S1
  null applies to UNIVERSAL; S2 is defined for PATHWAY only.
- **It is not evidence of correctness.** By design, context concentration is not
  evidence that a function is correct, and the narrowed proteins were not
  checked against any label.
- **It does not generalize across relations.** Neither PATHWAY nor ORGANISM
  shows the same pattern.

It is recorded as an observation for future work. It is not a finding of WP4.

## Limitations and threats to validity

- **Power.** The labelled ledger is dominated by assertions without a DIAMOND
  candidate or without joint coverage for the labelled function in the distant
  cohort. The narrowing tests therefore lack the assertions they need, whatever
  the relatives contribute. A larger strict-label set, or candidates beyond
  DIAMOND/Swiss-Prot, would be required.
- **Small test split.** Only three test targets have peers; T6 supports no
  positive claim on its own (design §8), and the McNemar tests have few
  discordant pairs.
- **Engine defaults (D5).** With `cladeCap = 1` and `jointMinEffective = 0.5`,
  congeners share one genus weight. Joint coverage from congeners needs about
  half of them to carry a reaction. Other weightings were not tested, by design;
  the result applies to the defaults only.
- **Frozen models (D3).** The models never saw close-relative features in
  training. A refit on R2 features could behave differently and is outside WP4.
- **Shifted inputs (D4).** The fenced gapsmith reference changes R0 inputs
  relative to 10406 for every target (T7). The frozen models see shifted inputs
  in R0 and R2 alike.
- **Family assignment.** Peer proteins were assigned to the fixed cohort
  orthogroups and the rest clustered among themselves, not jointly
  re-clustered (design §3).
- **Coverage gaps in the generated tables.** The C-peer family columns of T1 are
  NA for most targets. T8 harm versus SEQONLY is not evaluable without a SEQONLY
  removal arm. T9 recovery covers logistic only. `metrics.json` records the
  label policy hash as null; the labels file hash is recorded.
- **Deviations from the design counts.** The add-one block has 35 peers instead
  of 36 because one Pasteurellaceae peer is tier F under current NCBI taxonomy.
  The T2 block bootstrap uses only the contributing family blocks.
- **Seed fix scope.** The `content:v2` seed changes permutation draws in every
  regime. WP4 v2 is internally consistent (all arms on one JAR). Earlier cohort
  results produced with hold-out-keyed seeds were not re-run here (see
  [the permutation note](v3-cohort-permutation.md#holdout-independent-permutation-seed-contentv2)).

## Execution

The planner wrote 630 arm rows, packed into 420 integration runs. Base arms use
fenced inputs; `ref10406` arms use the 10406 inputs. Fragment copies are packed
per target, and duplicates reuse their source's fenced gaps. The add-one block
has 105 rows (35 peers × 3 relations). Each run left a `verified_run.json`
receipt. `seal` requires a verified receipt for every arm, and `seal.json` binds
420 distinct receipts, one per run directory, all with status `verified`.

| Stage | Where | Job |
|---|---|---|
| Pilot (Methanosarcina, R0 / R0-SC / R1 / R2) | Shaheen III | 15691175 |
| v1 seal → score → report (invalid: R0-SC defect) | Shaheen III | 15692636 |
| Seed fix: JAR bundle `ec-prefix-identity-v4`, focused spec 7/7, full suite 1,480 tests, 0 failures, 6 skipped | unimatrix01 | 10699 |
| Re-plan onto the v4 JAR + leakage audit (420 commands on v4, 0 on the old JAR) | Shaheen III | 15699717 |
| 420 integration runs | Shaheen III | 15700359 onward |
| v2 seal → score → report | Shaheen III | 15700443 |

The v1 trees on the cluster were renamed `*-v1-seedbug`. The v2 sealed table
(68,048,314 rows) stays on Shaheen at the path in `PROVENANCE.txt`.

## Provenance

`benchmark/v3/results/wp4-recurrence-v2/` (`SHA256SUMS` verified locally):

| Object | SHA-256 |
|---|---|
| `seal/seal.json` | `7f9f97424c667474c81b464319640b0e298f525529214d58c6ebad77327207ad` |
| sealed table `sealed_table.tsv` (Shaheen only) | `5bb5dd9d1b178d261335d0f367cc2b90933568673b8533c4ea6421e51a7a6be1` |
| `seal/sealed_gaps.tsv` | `0e3de1af7c849578c87bd71936d747cde98ed737d433571e5c930886772b8905` |
| `scores/metrics.json` | `8501cfaf80d8bb94296c8bae3c42dce6351f70ed7de4518ad8a7f33b3b895608` |
| `scores/score_receipt.json` | `3a6e06a54592ce6059d16e9d83c3893fa35f68a5217345ad8ddae37d37a07ece` |
| `report/v3-wp4-recurrence-tables.md` | `cda9cc4dd908f7e425448dde6b1aa105e5b569f87a2d456f606cbc9bcee0b0ef` |
| `report/v3-wp4-negative-results.md` | `deaadc50162255cf1387cfe8d8575a516faf64cb5134f2de538e5ea925a69c45` |
| `leakage-archaea.json` | `a65975ddb3c7120a95b09ae44b4b10f72a07a97bb36272dc98f662ccfd1a14fb` |
| `leakage-bacteria.json` | `3683a3a54fc01bb8a25a528cb7f62a19f77f28003eee03bbbc0e833153a3b6e7` |
| `PROVENANCE.txt` | `f244d1ffd7054003f99c078ec3774e4e9063d7b065b3e269df14118b9d809d5f` |
| Analyzer `benchmark/v3/analyze_wp4_recurrence.py` (v1 and v2 seals, current file) | `ee6a25183758db235e05097cbca3f01d2893bbd5b6d1270cb0ff9360f5c26960` |
| Leakage auditor `audit_wp4_leakage.py` | `db0e3b4e037d78b9eea3bba725c0b3f2a13759168d0d273e9e6789c68812a6d0` |
| Labels `training-cohort/labels/labels.tsv` | `7be9a66904d214be03abcc580cffb8430c443a6dc6d4b2fdd1ea80dd7b1afdc3` |
| JAR `ec-prefix-identity-v4/jobs/10699/gspa-3.0.0.jar` | `56fe479bd64a2acf4a3c29c3734893e3ffb52d42bf163ddeba986456393c58ac` |

The invalid v1 run is archived in `benchmark/v3/results/wp4-recurrence-v1/`
(seal `9e960057…`, sealed table `da649811…`, own `SHA256SUMS`). Its results must
not be interpreted.

The archive verifies file integrity and provenance. It is not biological
acceptance of any function call.
