# v3 training cohort: acquisition, label join and homology baseline

Status (2026-10-26): **the 19 selected proteomes were acquired, 18 resolved to exact assemblies, labels were
joined by verified identity, and the split and family audit passes. No feature export exists and no context
model was trained.** The inputs are not adequate to fit a stacker (531 labelled proteins, 97 in the test
split), so no trained-model claim is made anywhere. Everything below is reproducible with the commands in
section 5; committed outputs are in `benchmark/v3/results/training-cohort/`, bulky data are in
`~/.cache/gspa-v3-data/training-cohort/`. The selection (`results/label-rich-selection/`) was read and left
unchanged.

| Piece | File |
| --- | --- |
| Proteome fetch, assembly resolution/download, Price audit, layouts | `benchmark/v3/acquire_training_cohort.py` |
| Swiss-Prot label join (accession + sequence, or unique full identity) | `benchmark/v3/match_swissprot_labels.py` |
| Splits, protein list, MMseqs2 families, split audit, masking preview | `benchmark/v3/prepare_training_cohort.py` |
| Homology baseline with the cohort fenced off, scoring, claim export | `benchmark/v3/homology_baseline.py` |
| Tests (48) | `benchmark/v3/tests/test_acquire_training_cohort.py`, `test_swissprot_label_matching.py`, `test_homology_baseline.py` |

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

## 5. Commands

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
$PY -I -m unittest discover -s $B/tests -p 'test_*.py'      # whole v3 suite; 4 Price tests import via the repo root
```

The registry entry `swissprot_label_rich_v1` in `registry/training_cohort_registry.json` still states
`genomes_acquired: 0`; it was not edited here because it is shared. It should read 18 genomes, 499 strict-label
proteins matched.

## 6. What is next

1. Build the shadow jar (`./gradlew :gspa-cli:shadowJar`) and run `gspa integrate` with the per-genome claim files
   in `baseline/claims/homology/` and `--feature-export-dir`, then `prepare_training_data.py build-table` with
   `labels.tsv`, `splits.tsv`, `example_clusters.tsv` and a label policy that declares `experimental` for both GO
   and EC and `negative_policy: none`. Not run here: the manifest, ontology and essential-profile inputs of
   `integrate` have not been prepared for this cohort.
2. Replace the unresolved *M. gryphiswaldense* only by a documented, hash-checked assembly decision (for example the
   strain whose proteins the Swiss-Prot entries actually match), or accept four test genomes.
3. Declare a closed-world or candidate-ranking evaluation separately from exact-EC accuracy: with 97 test proteins
   and sparse experimental positives, only the latter two are defensible, and an unconfirmed alternative is never a
   proven negative.
4. Training a stacker needs far more labelled genomes than exist locally; this cohort is an evaluation and
   calibration reservoir. The existing registry/data docs list the other cohorts and why they do not supply context.
