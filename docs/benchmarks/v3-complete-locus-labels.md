# Complete-locus labels for the learned stage (row 4), 2026-11-02

Frozen plan: `tasks/v3-row4-complete-locus-labels-2026-11-02.md` (sha256 `9d23a106…`, recorded in
`benchmark/v3/results/complete-locus-labels-v1/PLAN_SHA256` at 13:21 before the first count; dated
deviations appended). The registry sha256 values were recorded in `REGISTRY_SHA256` at 13:36, before
any score was read. Results: `benchmark/v3/results/complete-locus-labels-v1/` (`SHA256SUMS`). Raw
work: `/mnt/data2/leechuck/gspa-v3-data/complete-locus-labels-v1/` on ws. Code:
`benchmark/v3/complete_locus/`. Compute: ws only; nothing ran on Shaheen.

A competition is every protein that one candidate source proposes for one function in one genome,
at any score. It is complete when every member has an experimental label with a PMID: positive
(Swiss-Prot 2025_03 experimental GO or ECO:0000269 EC; Price-148 for the Price genomes) or an
experimentally supported negative that the curator accepted after reading it. Unconfirmed members
are never negatives; one unlabelled member makes the competition incomplete. Membership came from
the claim and export files only; no score was read before the registry was frozen.

## Verdict against the frozen criteria

| id | criterion | result |
|---|---|---|
| C0 | plan sha before counts, registry sha before scores | PASS (13:21 and 13:36; the refit and the first score read came afterwards) |
| C1 | every label has accession, release, evidence code, PMID | PASS: 1,156 labels in 941 complete competitions; 3 positives without a PMID were excluded first, which made 3 competitions incomplete |
| C2 | no negative without an accepted review | PASS: 0 validator violations; 7 negatives accepted, 3 rejected |
| C3 | ≥ 20 complete competitions usable by the null layer (0 or 1 positive) | PASS: 812 |
| C4 | ≥ 5 of them in calibration or test genomes | PASS: 23 (cohort 11 calibration + 7 test; Price 2 calibration + 3 test), 16 distinct positive proteins |
| C5 | strict reading: contested (≥ 2 members), ≥ 20 and ≥ 5 in calibration/test | FAIL: 6 contested, 0 in calibration/test |
| C6 | builder certifies every declared cohort group, apart from masking losses | PASS: 66 of 68 certified; 2 training groups (EC 2.5.1.41) lost a member to family masking |
| C7 | null layer trained by the unchanged trainer at its default gates | FAIL: 0 all-negative training groups (gate 20); 48 single-positive (gate 20) |
| C8 | test NLL below the uniform baseline | NOT APPLICABLE (no null layer); the descriptive evaluation is below |

The plan's targets for the count are met (C3, C4), but only because single-candidate competitions
are allowed. These do not test ranking, and none of them is an all-negative case. The null layer
cannot be fitted from public data: across all 43 genomes there is no complete competition in which
no candidate is correct.

## Complete competitions

| set | split | single positive, one member | single positive, contested | all negative | several positives | distinct positive proteins (usable) |
|---|---|---:|---:|---:|---:|---:|
| cohort (18) | train | 49 | 1 | 0 | 2 | 36 |
| cohort | calibration | 11 | 0 | 0 | 0 | 8 |
| cohort | test | 7 | 0 | 0 | 2 | 3 |
| model genomes (3) | held out, development-examined | 734 | 5 | 0 | 125 | 498 |
| Price (22) | calibration | 2 | 0 | 0 | 0 | 2 |
| Price | test | 3 | 0 | 0 | 0 | 3 |
| total | | 806 | 6 | 0 | 129 | 550 |

Competitions overlap: one protein can complete several GO terms (for example AAM33000.1 completes
four nested aminoacylation terms). The model genomes have no split, so they are not counted as
calibration or test. Their universe is the role-fenced DIAMOND claim set, which the current learned
models do not support. The Price universes come from the TopEC, MAPred and EnzBERT claims. All five
Price competitions are TopEC singletons; MAPred and EnzBERT propose each EC on too many proteins to
complete any. Only the 66 certified cohort competitions can be scored by the learned stage.

The six contested competitions with exactly one positive:

| genome | function | positive | accepted negative(s) | negative source |
|---|---|---|---|---|
| S. solfataricus (cohort train) | GO:0004829 threonine-tRNA ligase | ThrRS-cat AAK42622.1 | ThrRS-ed AAK40712.1 | PMID 15240874 |
| B. subtilis | EC 6.3.1.20 and GO:0016979 lipoate-protein ligase | LplJ O07608 | LipM P54511 | PMID 20882995 |
| B. subtilis | GO:0004619 phosphoglycerate mutase | GpmI P39773 | YhfR O07617 | PMID 10869096 |
| E. coli | GO:0051989 coproporphyrinogen dehydrogenase | HemN P32131 | HemW P52062 | PMID 29282292 |
| M. tuberculosis | GO:0003938 IMP dehydrogenase | GuaB2 P9WKI7 | GuaB1 P9WKI3, GuaB3 P9WKI5 | PMID 21081761 |

Negatives reviewed (`registry/negative_reviews.tsv`, PubMed records hashed in
`inputs/pubmed_review_SHA256SUMS`):

| protein | function | disposition | reason |
|---|---|---|---|
| ThrRS-ed AAK40712.1 | GO:0004829 | accepted | no threonylation in the assays |
| ThrRS-ed AAK40712.1 | GO:0006435 (process) | rejected | an editing subunit can take part in the process; not denied by an activity assay |
| LipM P54511 | EC 6.3.1.20, GO:0016979 | accepted (2) | octanoyltransferase; "does not display lipoate protein ligase activity" |
| YhfR O07617 | GO:0004619 | accepted | "little if any PGM activity" |
| HemW P52062 | GO:0051989 | accepted | "does not have coproporphyrinogen III dehydrogenase activity in vitro" |
| GuaB1 P9WKI3, GuaB3 P9WKI5 | GO:0003938 | accepted (2) | "inactive in vitro" |
| GlpF P0AER0 | GO:0006833 water transport | rejected | Swiss-Prot reports low but non-zero water permeability |
| Rv2607 P9WIJ1 | GO:0004733 | rejected | denies only the PMP substrate; experimental EC 1.4.3.5 (PNP oxidation) is in the same term |

Rejecting Rv2607 removed the only all-negative competition found anywhere (a singleton).

## Why no all-negative competitions exist, and what is missing

| quantity | count |
|---|---:|
| competitions examined (cohort / model / Price) | 84,843 / 23,167 / 205,384 |
| ... with at least one experimental positive | 231 / 2,391 / 257 |
| experimental negatives joined to a protein of the 43 genomes (protein, function pairs) | 178 (cohort 40, model 137, Price 1) |
| ... whose function the source proposes for that protein (a competition contains it) | 35 competitions |
| ... complete before review | 8 (7 contested, 1 singleton); 7 after review |
| one unlabelled member short of complete (all-negative if that member were negative) | 5 cohort training, 8 model |
| cohort negatives in calibration / test genomes | 1 / 0 (no competition contains them) |

The trainer needs at least 20 all-negative and 20 single-positive training competitions, plus 5 of
each in calibration, before it fits the null layer. Even if every competition that is one member
short had that member experimentally confirmed negative, there would be at most 5 cohort training
and 8 model-genome all-negative competitions, and none in calibration. The limit is the public
negative data itself: Swiss-Prot and GOA record about 750 curated experimental activity negatives
across all prokaryotes, and 17 cohort proteins carry one. The calibration and test genomes
(Campylobacter, Haemophilus, Methanothermobacter; Methanosarcina, Dickeya, Clostridium,
Paramagnetospirillum) have essentially none. Reaching the gate would need one of the following:

- experimental negative assays (purified-protein activity tests or knockout complementation) that
  cover every proposed candidate of at least 20 functions in training genomes and 5 in calibration
  genomes, where no candidate has the function;
- a split that moves heavily characterized genomes such as E. coli, B. subtilis and M. tuberculosis
  into training and calibration, with a learned-stage export under the training predictor contract.
  This is a new split and a new export, not the existing ones, and it would still give at most about
  8 all-negative competitions.

## Descriptive evaluation (frozen family-excluded models, `models/default/`)

The refit with the new truth metadata reproduces the frozen default models exactly: the scoring,
calibration and threshold parameters are identical. The test metrics also repeat the published
ones (logistic: Brier 0.0908, AUROC 0.743, precision 0.25 at threshold 0.298; tree: Brier 0.0884,
AUROC 0.752, precision 0.286). This is expected, because the cohort labels are unchanged; only the
locus certificates are new. All 18 complete calibration/test cohort competitions are scored by both
models. Each has one experimentally confirmed positive and no other candidate, so in all 18 the true
answer is "this candidate", never "none".

| model | threshold | calibration: positive called | test: positive called | probability range of the confirmed positive (test) |
|---|---:|---:|---:|---|
| logistic | 0.298 | 4 of 11 | 0 of 7 | 0.021–0.279 |
| tree | 0.216 | 6 of 11 | 1 of 7 | 0.039–0.343 |

On the test genomes, the marginal models leave nearly every experimentally confirmed single
candidate below their display threshold. The evaluation has no contested or all-negative group, so
it measures recall on certain loci, not ranking or abstention.

## Files

- `benchmark/v3/results/complete-locus-labels-v1/`: `plan_frozen.md`, `PLAN_SHA256`,
  `REGISTRY_SHA256`; `registry/` (`complete_competitions.jsonl` with per-label provenance,
  `excluded_competitions.jsonl`, `negative_reviews.tsv`, cohort builder inputs `labels.tsv`,
  `locus_groups.tsv`, `label_policy.json`, step-1 summary and review queue); `fit/` (builder
  manifest with the locus audit, both frozen refits, `training.json`, `evaluation.json`, table hash);
  `audit/` (all 35 competitions containing a negative, step-1 output hashes and the identical rerun);
  `inputs/` (Swiss-Prot scan receipt, GOA QuickGO receipt, export restore receipt, negative-source
  inventories and hashes).
- Code: `benchmark/v3/complete_locus/{restore_export,scan_swissprot,build_candidates,fetch_goa_pmids,finalize_registry,run_training}.py`.

Reproduce (new output directories; the PMID download and PubMed retrieval are recorded and should
not be repeated against a changed GOA):

```bash
export TMPDIR="$HOME/tmp"; PY="$HOME/tmp/gspa-v3-testenv/bin/python"; W=<new work dir>
$PY benchmark/v3/complete_locus/restore_export.py --archive benchmark/v3/results/family-excluded-training-v1 --out $W/features
$PY -I benchmark/v3/complete_locus/scan_swissprot.py --swissprot <uniprot_sprot.dat.gz a6907c4d…> --accessions <accessions.txt> --faa <price all_proteins.faa> --out $W/entries.jsonl
$PY benchmark/v3/complete_locus/build_candidates.py ... --out $W/step1          # see the plan for inputs
$PY benchmark/v3/complete_locus/finalize_registry.py --step1 $W/step1 --reviews negative_reviews.tsv --goa <QuickGO dir> ... --out $W/registry
$PY benchmark/v3/complete_locus/run_training.py {build,train,evaluate} --features $W/features --registry $W/registry --go-obo <go-basic.obo> --default-models benchmark/v3/results/family-excluded-training-v1/models/default --out $W/fit
```
