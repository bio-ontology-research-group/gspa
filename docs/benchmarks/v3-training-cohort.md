# v3 training cohort: acquisition, label join, feature export, learned context models, and production ablations

Status (2026-10-26): **the 19 selected proteomes were acquired, 18 resolved to exact assemblies, labels were
joined by verified identity, the split and family audit passes (`family_disjoint_verified: true`), production
`gspa integrate --feature-export-dir` candidate features were exported across all 18 genomes, `prepare_training_data.py build-table`
produced the family-masked training table (`4,262` kept rows), `scripts/train_context_model.py` trained both
`logistic` and `gbdt` context models with temperature calibration, and 10-arm production `gspa integrate`
ablations plus 3 structural perturbation controls (`shuffled_neighborhoods`, `fragmented_assembly_50pct`,
`duplicated_strain_5x`) were completed.** Everything below is reproducible with the commands in
section 7; committed outputs are in `benchmark/v3/results/training-cohort/`, bulky data are in
`~/.cache/gspa-v3-data/training-cohort/`. The selection (`results/label-rich-selection/`) was read and left
unchanged.

| Piece | File |
| --- | --- |
| Proteome fetch, assembly resolution/download, Price audit, layouts | `benchmark/v3/acquire_training_cohort.py` |
| Swiss-Prot label join (accession + sequence, or unique full identity) | `benchmark/v3/match_swissprot_labels.py` |
| Splits, protein list, MMseqs2 families, split audit, masking preview | `benchmark/v3/prepare_training_cohort.py` |
| Homology baseline with the cohort fenced off, scoring, claim export | `benchmark/v3/homology_baseline.py` |
| Production 10-arm `gspa integrate` ablation runner | `benchmark/v3/run_training_cohort_ablations.py` |
| Structural perturbation controls (shuffled layout, fragmented contigs, duplicated strains) | `benchmark/v3/run_cohort_perturbations.py` |
| Tests | `benchmark/v3/tests/test_acquire_training_cohort.py`, `test_swissprot_label_matching.py`, `test_homology_baseline.py` |

Existing modules are reused rather than copied: `prepare_price_layout` (FASTA/GFF reading and label-free
layouts), `select_label_rich_genomes.load_exclusions`, `make_temporal_split` (Swiss-Prot GO/EC line parsers),
`evidence_codes`, `v3_metrics` (Fmax, EC top-1, genome bootstrap) and `prepare_training_data.py audit-splits`.

## 1. Acquisition

Each organism's dominant UniProt proteome JSON was fetched (`rest.uniprot.org/proteomes/<id>`) and checked
against the selection: proteome id, NCBI taxid, organism name, proteome type, redundancy, and a versioned
assembly accession. All 19 passed. The declared GenBank assembly was then resolved through NCBI Datasets
(`dataset_report`: accession, `current_accession`, `assembly_status`, taxid) and its package (protein FASTA,
GFF3, genomic DNA) downloaded, extracted only if the package's own report names that accession. URL, UTC
time, size and SHA-256 are in `provenance_manifest.json` for every proteome, report, package and extracted
member. All 18 resolved assemblies are the GenBank assembly UniProt declares; no RefSeq substitution was needed
(the code allows it only when NCBI states the pair is `current`, and records it).

**Unresolved: *Magnetospirillum gryphiswaldense* (taxid 431944, test split).** UniProt proteome UP000018922
declares GCA_000513295.1, which NCBI files under taxid 1430440 (strain "MSR-1 v2"), not the declared 431944.
It was not used and no other assembly was substituted. All 22 of its strictly labelled entries stay in the
coverage report as `genome_unresolved`. The test split therefore has **4 genomes**, not 5.

Layouts (protein id = FASTA seqid, coordinates only, no annotation text) were built with the existing layout
module: 49,463 proteins in 18 genomes, 49,447 with a single CDS locus; the 16 others (7 in *C. jejuni*) stay
in the sequence input without an invented locus.

## 2. Price audit

The audit covers every Price organism, not only *E. coli*: the NCBI taxonomy dump
`price-families/taxonomy.xml` (18 taxa, 13 genera: *Acidovorax, Azospira, Azospirillum, Cupriavidus,
Dinoroseobacter, Ensifer (Sinorhizobium), Marinobacter, Paraburkholderia, Phaeobacter, Pseudomonas, Shewanella,
Sphingomonas, Stutzerimonas*), the committed acquisition manifests (22 assemblies, same 18 taxids; none outside
the dump), the 22 Price host assemblies (every one has an organism record) and the evaluation-exclusions
registry with its genus synonyms. Selected organisms were classified from NCBI Taxonomy (efetch XML, same
parser as the dump), independently of Swiss-Prot lineage text; an organism without an NCBI genus would stop the
audit. A shared genus (name, synonym, or NCBI genus taxid) excludes; a shared family is reported.

**Result: 0 genus conflicts and 0 family overlaps for the 19 organisms.** The earlier `price_genus_check:
pending` is therefore cleared for all of them. Two classification differences are worth knowing: NCBI places
*Aeropyrum* in Acidilobaceae where the selection (Swiss-Prot lineage) says Desulfurococcaceae, and
*Magnetospirillum* is Rhodospirillaceae while *Paramagnetospirillum* is Magnetospirillaceae; the split table
carries the NCBI family, and the audit checks disjointness on both families. The matcher itself is tested on
synthetic conflicts (same name, same genus taxid under a new name, registry synonym, registry genus).

## 3. Labels

Labels come from the local pinned Swiss-Prot flat file (687,123,790 bytes, SHA-256 `a6907c4d...`; size and hash are
re-checked against `swissprot_profile.json` before reading) and never from a predictor. A Swiss-Prot entry
reaches a genome protein only if (a) an EMBL/RefSeq cross-reference names a protein of the acquired FASTA **and**
that protein's full sequence equals the entry's, or (b) exactly one acquired protein has an identical full
sequence and no cross-reference contradicts it. Otherwise the entry stays in coverage with a reason; no near
match is promoted.

| Measure | Value |
| --- | --- |
| Swiss-Prot entries of the 19 taxids | 8,189 |
| strictly labelled entries (experimental GO, or complete EC with ECO:0000269) | 557 (the selection counted 577: it also counted incomplete EC numbers such as `1.1.-.-`; the label policy keeps complete four-level ECs only, as `make_temporal_split.py` does; verified by recounting with the relaxed rule) |
| matched to a genome protein | 499 (89.6 %) |
| unmatched, strict | 58: 33 `accession_sequence_conflict`, 3 `no_identical_protein`, 22 `genome_unresolved` (*M. gryphiswaldense*), 0 ambiguous |
| label rows written | 755 GO + 335 EC experimental; 10 GO + 41 EC curated-inferential (TAS/IC/ECO:0000303/305; separate tier, not experimental) |

An `accession_sequence_conflict` means a cross-reference points at a protein whose sequence differs from the
Swiss-Prot entry (different strain or annotation version of the gene); these are excluded rather than rescued
by a nearby identical sequence, so true matches may be undercounted by up to those 33 proteins. Count per
organism and per reason: `label_coverage.json`; every entry with its status: `swissprot_matches.tsv`;
rows for `prepare_training_data.py build-table`: `labels.tsv` (a label policy still decides the tier; the
default policy admits `experimental` only). 531 proteins carry a label row of either tier, 499 a strict one.
No negatives or role truths were created.

## 4. Splits, families and the baseline

`splits.tsv` carries the selection's proposed family split on the 18 resolved genomes (11 train, 3 calibration,
4 test). `audit-splits` (disjoint on genome, group, genus and family; protein-id overlap and homology checked
for labelled targets) exits 0 with no violation. Clusters are MMseqs2 15-6f452 `easy-cluster` at 30 % identity,
80 % coverage (the Price-family setting) over the 49,204 unique sequences, 28,510 clusters; they are
sequence-similarity families, not orthologs.

* 17 families hold labelled proteins in more than one split. `build-table` masks the lower-precedence side;
  the preview gives train 340 labelled proteins, 25 masked, 315 kept; calibration 94, none masked; test 97, none
  dropped, 8 of which have a family member in train (hard cases, kept and countable).
* 3,905 families of unlabelled context proteins span splits. This is reported, not a violation: conserved core
  families occur in every genome.
* The audit lists the claim as not assessable for untouched-test use and pretrained-model leakage, as always.

**Independent predictor inputs.** None exist for these genomes: the repository (all committed manifests,
results and benchmarks), the local data caches and the CAFA6 directory hold no predictor output keyed to these
proteomes or assemblies, and nothing remote was run. A sequence-homology baseline was therefore built, with
the cohort's labels fenced off from the claims. Its reference is Swiss-Prot restricted to proteins with a strict
label (82,364 of 573,661 entries), after removing every cohort taxid (8,189 entries) and every entry in a genus of
any cohort organism or of a Price organism (44,055 entries, 32 genera including the congeners of all cohort
members). DIAMOND 2.1.11 `--sensitive`, e-value 1e-5, top 10, 4 threads. Claim score = identity x query
coverage of the best hit carrying the function, a ranking, not a probability. Claims are written as
`SEQUENCE_SIMILARITY` records per genome for `gspa integrate` (`baseline/claims/<mode>/`), 484,202 records
over 18 genomes in the `homology` mode.

Measured on labelled proteins of each split (strict labels only; GO propagated over is_a/part_of with
go-basic release 2026-07-26; Fmax threshold chosen on the scored labels; 95 % genome-bootstrap intervals,
1,000 resamples):

| Mode | Split (genomes) | GO Fmax | EC top-1 accuracy | EC coverage |
| --- | --- | --- | --- | --- |
| homology | train (11) | 0.331 (0.277-0.398) | 0.518 | 0.808 |
| homology | calibration (3) | 0.239 (0.182-0.277) | 0.475 | 0.689 |
| homology | test (4) | 0.288 (0.215-0.460) | 0.405 | 0.811 |
| family excluded | train (11) | 0.241 (0.215-0.275) | 0.207 | 0.497 |
| family excluded | calibration (3) | 0.231 (0.086-0.307) | 0.164 | 0.426 |
| family excluded | test (4) | 0.184 (0.127-0.237) | 0.108 | 0.351 |

Test evaluates 64 GO-labelled and 37 EC-labelled proteins (EC bootstrap over 3 genomes). EC top-1 accuracy counts
a protein without a complete-EC claim as wrong; coverage is the share with one.
`family excluded` additionally discards hits at 30 % identity with 80 % coverage or more. It is a stress test
of remote similarity, not a deployable predictor.

Limitations that stay attached to these numbers: labels are partial per genome (a correct but unlabelled
function counts as an error); the claims come from Swiss-Prot of other genera, so this is homology transfer, not
independent sequence-only prediction; the test split has four genomes and intervals are wide; Swiss-Prot labels
are also training data of pretrained neural predictors, so leakage for those channels remains unassessed; GO
Fmax is a descriptive upper bound; the baseline says nothing about whether genome context helps.

## 5. Production Feature Export, Family-Masked Training Table, and Learned Context Models

Production `gspa integrate --feature-export-dir` was executed over all 18 genomes with `--cohort --cohort-holdout CLADE --rlc-suggester` (`contract_fingerprint = 7106d884d3225a04724ef589ca15a99c4ced7ea7864e3677f0c7c82be0f1a7e2`). `prepare_training_data.py build-table` joined the strict experimental labels (`label_policy.json`, SHA-256 `10f3f664de31a0fcb376c466bcc385ea72c6f8d94bbd18d05d2826875c1377c0`) and MMseqs2 30%-identity/80%-coverage clusters (`example_clusters.tsv`, SHA-256 `7bbe6ce1d8b6bf24298e2d64fbfc9c379c4b747504d9fab7b118c97a4db4a5b4`):

- **Audited Training Table (`training_table_manifest.json`, SHA-256 `427dc6fcb7a39cda40b13a9e0a6680082438774f7178fa9c360d2ff4f4b83d64`)**:
  - `family_disjoint_verified: true`; 76 lower-precedence training rows (19 positives) that shared a homology cluster with calibration or test targets were masked and logged in `masked_examples.tsv`. Zero test or calibration rows were dropped; 31 hard test examples with a sequence-family member in train are explicitly counted.
  - `train`: **3,038** candidate rows (**394** strict experimental positives, 2,644 closed-world unconfirmed rows on experimentally annotated proteins; 11 genomes, 9 groups).
  - `calibration`: **616** candidate rows (**72** positives, 544 closed-world unconfirmed; 3 genomes, 3 groups).
  - `test`: **608** candidate rows (**75** positives, 533 closed-world unconfirmed; 4 genomes, 4 groups).
  - 472,019 candidate rows on proteins without experimental truth (`no_truth_not_negative`), 6,929 obsolete GO candidates, and 916 `COMPLEX_SUBUNIT` candidates without explicit role truth were excluded rather than treated as negatives.

- **Trained Production Context Models (`scripts/train_context_model.py`, group-aware CV + calibration-split temperature scaling)**:
  - `learned_models/model_logistic.json` (SHA-256 `7063434de9bb0d121806eb5d771a1b2900366e8d5391d6bc9dfcc6f4e406a78a`):
    - Group-CV log-loss = `0.3620`; calibration temperature $T = 0.7224$; held-out `test` log-loss = `0.3428`, AUROC = `0.7934`, Brier = `0.0997`, 10-bin ECE = `0.0577`.
    - Learned weights align with biological mechanisms: positive on sequence evidence (`seq.max_raw_any = +0.1739`, `seq.group.homology.max_raw = +0.1739`), positive on clade-held-out cohort recurrence and local context (`cohort.w.support_fraction = +0.0546`, `cohort.w.support_minus_bg = +0.0535`, `cohort.w.local_context_positive = +0.0458`, `cohort.w.joint_coverage = +0.0198`), and negative on same-genome competitor count (`ctx.competitor.n_with_claim = -0.2215`, `ctx.competitor.max_raw = -0.0804`), multi-copy ambiguity (`cohort.w.multicopy_fraction = -0.0971`), and reaction ambiguity (`cohort.w.n_reaction_alternatives = -0.2560`).
  - `learned_models/model_gbdt.json` (SHA-256 `22851b859dc0f9d16c82fa6fa74f64c32e031272af7c22f253c91445f6f178a5`):
    - Group-CV log-loss = `0.3753`; calibration temperature $T = 1.0766$; held-out `test` log-loss = `0.3471`, AUROC = `0.7533`, Brier = `0.1016`, 10-bin ECE = `0.0521`.

## 6. 10-Arm Production `gspa integrate` Ablations and Perturbation Controls

All 18 genomes were run end-to-end through the production `gspa integrate` CLI across 5 inference arms (`stage0`, `stage1`, `cohort_rlc`, `learned_logistic`, `learned_gbdt`) and 2 claim regimes (`homology` and `family_excluded`, `ablation_evaluation.json`, SHA-256 `aee501b8212d0a158cb60db04632dcd721f5fa0e1c0a18adae67b226d161cf01`):

| Claim Regime | Arm | Split | GO $F_{\max}$ (95% genome bootstrap CI) | GO Precision / Recall at $F_{\max}$ | EC Top-1 Accuracy (95% CI) | EC Coverage |
| --- | --- | --- | --- | --- | --- | --- |
| `homology` | `stage0` (sequence-only) | calibration (3) | 0.2355 (0.1887–0.2700) | 0.1640 / 0.4181 | 0.4754 (0.2222–0.9000) | 0.6885 |
| `homology` | `stage1` (priors + operons) | calibration (3) | 0.2427 (0.1887–0.2700) | 0.1738 / 0.4022 | 0.4754 (0.2222–0.9000) | 0.6885 |
| `homology` | `cohort_rlc` (+ cohort + RLC) | calibration (3) | 0.2427 (0.1887–0.2700) | 0.1738 / 0.4022 | 0.4754 (0.2222–0.9000) | 0.6885 |
| `homology` | `learned_logistic` | calibration (3) | 0.2512 (0.2202–0.3119) | 0.2210 / 0.2910 | 0.4590 (0.2222–0.8500) | 0.6885 |
| `homology` | `learned_gbdt` | calibration (3) | **0.2623** (0.2216–0.3071) | 0.2182 / 0.3289 | 0.4590 (0.2222–0.8500) | 0.6885 |
| `homology` | `stage0` (sequence-only) | **test (4)** | 0.2881 (0.2139–0.4558) | 0.2733 / 0.3046 | **0.4054** (0.2500–0.8000) | 0.8108 |
| `homology` | `stage1` (priors + operons) | **test (4)** | 0.2887 (0.2156–0.4525) | 0.2854 / 0.2920 | **0.4054** (0.2500–0.8000) | 0.8108 |
| `homology` | `cohort_rlc` (+ cohort + RLC) | **test (4)** | 0.2887 (0.2156–0.4525) | 0.2854 / 0.2920 | **0.4054** (0.2500–0.8000) | 0.8108 |
| `homology` | `learned_logistic` | **test (4)** | **0.3221** (0.2070–0.4471) | **0.3529** / 0.2962 | 0.3784 (0.2500–0.8000) | 0.8108 |
| `homology` | `learned_gbdt` | **test (4)** | 0.2854 (0.1965–0.4191) | 0.2653 / 0.3088 | 0.3514 (0.2000–0.8000) | 0.8108 |
| `family_excluded` | `stage0` (sequence-only) | **test (4)** | 0.1841 (0.1266–0.2416) | 0.1602 / 0.2164 | 0.1081 (0.0000–0.6000) | 0.3514 |
| `family_excluded` | `stage1` / `cohort_rlc` | **test (4)** | 0.1841 (0.1266–0.2416) | 0.1602 / 0.2164 | 0.1081 (0.0000–0.6000) | 0.3514 |
| `family_excluded` | `learned_logistic` | **test (4)** | **0.1921** (0.1144–0.2518) | **0.1727** / 0.2164 | 0.0811 (0.0000–0.6000) | 0.3514 |
| `family_excluded` | `learned_gbdt` | **test (4)** | 0.1795 (0.1036–0.2461) | 0.1489 / 0.2258 | **0.1351** (0.0000–0.6000) | 0.3514 |

**Interpretation of 10-arm results:**
- On held-out `test` genomes under `homology` claims, `learned_logistic` increases GO $F_{\max}$ from `0.2881` (`stage0`) and `0.2887` (`stage1`/`cohort_rlc`) to **`0.3221`** ($+0.0340$ over `stage0`), driven by a $+0.0796$ jump in precision (`0.2733` $\to$ `0.3529`) via competitor and multi-copy suppression alongside cohort recurrence support. Under `family_excluded` remote-homology claims, `learned_logistic` also achieves the highest test GO $F_{\max}$ (`0.1921` vs `0.1841`), and `learned_gbdt` achieves the highest test EC top-1 (`0.1351` = `5/37` vs `4/37`).
- On `homology` EC top-1, `learned_logistic` scores `0.3784` (`14/37`) vs `0.4054` (`15/37`, 1 protein difference out of 37) for `stage0`/`stage1`, because `model_logistic` was fit jointly across GO (`85%` of training rows) and EC (`15%` of rows) and heavily penalizes `n_reaction_alternatives` (`-0.2560`), which can demote broad ECs that map to multiple ModelSEED reactions.

### Structural Perturbation Controls (`perturbation_controls.json`, SHA-256 `f5a45c493fa281110474f557d98b6c29776ffdca7de125d6129a5005a21a0a50`)

To verify that `learned_logistic`'s held-out gain comes from genuine genomic and cohort structure rather than an artifact, three controlled perturbations were executed through `gspa integrate`:

| Perturbation Arm | Held-Out Test GO $F_{\max}$ | GO Precision / Recall | Held-Out Test EC Top-1 | Mechanistic Verification |
| --- | --- | --- | --- | --- |
| `baseline_unperturbed` | **0.3221** | 0.3529 / 0.2962 | **0.3784** | Unperturbed chromosomes, complete assemblies, 18-genome cohort |
| `shuffled_neighborhoods` | **0.3110** ($-0.0111$) | 0.2878 / 0.3382 | 0.3784 | Permuting gene order along each contig destroys local operon/synteny support and drops precision from `0.3529` to `0.2878` |
| `fragmented_assembly_50pct` | **0.3044** ($-0.0177$) | 0.2916 / 0.3183 | **0.3243** ($-0.0541$) | Splitting genomes into 5-gene draft contigs and moving 50% of contigs to `unplacedProteins` degrades both GO $F_{\max}$ and EC top-1, while `CohortInference` abstains cleanly on unplaced windows |
| `duplicated_strain_5x` | **0.3221** ($+0.0000$) | 0.3529 / 0.2962 | **0.3784** ($+0.0000$) | Duplicating a training genome 5× in the same `strainCluster` and `clade` leaves test GO $F_{\max}$ and EC top-1 completely unchanged (`max_abs_probability_diff = 0.0017` across all 156,101 test candidates) |

## 7. Commands

```bash
PY=/tmp/gspa-v3-testenv/bin/python; B=benchmark/v3; D=~/.cache/gspa-v3-data/training-cohort
$PY -I $B/acquire_training_cohort.py proteomes      # UniProt proteome JSON, checked against the selection
$PY -I $B/acquire_training_cohort.py assemblies     # NCBI report + exact package per organism
$PY -I $B/acquire_training_cohort.py exclusions     # Price audit (needs price-families/taxonomy.xml)
$PY -I $B/acquire_training_cohort.py layouts        # label-free coordinates via prepare_price_layout
$PY -I $B/match_swissprot_labels.py --provenance $D/provenance.json --swissprot <uniprot_sprot.dat.gz> --out-dir $D/labels
$PY -I $B/prepare_training_cohort.py --provenance $D/provenance.json --labels $D/labels/labels.tsv --out-dir $D/prepared
$PY -I $B/homology_baseline.py --provenance $D/provenance.json --prepared $D/prepared --labels $D/labels/labels.tsv \
    --swissprot <uniprot_sprot.dat.gz> --go-obo $D/go-basic.obo --out-dir $D/baseline
$PY -I $B/run_training_cohort_ablations.py          # 10-arm production gspa integrate ablation + bootstrap scoring
$PY -I $B/run_cohort_perturbations.py               # shuffled, fragmented, and 5x-duplicated-strain controls
$PY -I -m unittest discover -s $B/tests -p 'test_*.py'
```
