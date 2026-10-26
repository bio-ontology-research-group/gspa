# v3 training-data audit and preparation

Status (2026-10-26): **audit and preparation tooling only. No trained model exists and no
local cohort is sufficient to train and independently test a production learned context
model.** Everything below was read from files on disk; counts are reproducible with
`benchmark/v3/audit_training_cohorts.py` (report committed at
`benchmark/v3/results/training-data-audit/cohort_profile.json`). Price148 is a
**development** set and appears here only as such.

| Piece | Location |
| --- | --- |
| Split audit + training-table builder | `benchmark/v3/prepare_training_data.py` (`audit-splits`, `build-table`) |
| Cohort registry (role of every dataset) | `benchmark/v3/registry/training_cohort_registry.json` |
| Reproducible cohort profile | `benchmark/v3/audit_training_cohorts.py` |
| Entry-date split (fixed, relabelled) | `benchmark/neural/make_temporal_split.py` |
| Tests | `benchmark/v3/tests/test_prepare_training_data.py`, `test_audit_training_cohorts.py`, `test_make_temporal_split.py` |

## 1. What data exists (exact counts)

| Cohort | Content | Usable for | Not usable for |
| --- | --- | --- | --- |
| **CAFA6 reconstruction** (`benchmark/neural/cafa6_recon`, Kaggle `train_terms.tsv`) | 82,404 training proteins, 537,027 pre-t0 GO pairs; 25,000-protein "validation" fold; 224,309-protein test superset; later-snapshot experimental annotations 562,499 pairs on 88,894 proteins (1,204,013 evidence lines; IPI 396,441, IDA 330,393, TAS 134,746, IMP 181,103 ...) | protein-level, sequence-evidence-only GO candidate training | context features: no genome layout, operon, assembly or cohort. A model fitted here cannot learn what context does |
| **suite1 CAFA temporal manifest** | 24,505 rows: LK-Partial 16,233, LK-Limited 7,512, Aspect-NK 282, Protein-NK 251, EC holdout227 161, EC holdout66-strict 66 | protein-level regime labels | independent test: 23,721 of 24,254 distinct proteins (97.8%) are in the CAFA6 training set (all LK rows but 56; only the new *terms* are held out) |
| **SwissProt holdout227** | 227 entries created 2025-02-05 .. 2026-01-28, 121 organisms, all with EC; raw CLEAN / EnzBERT / DeepEC matrices (227 x 5,242 / 4,743 / 2,802) | protein-level EC hold-out for predictors trained before 2025 | genome context; time basis is entry creation, not annotation (see 2) |
| **Price148** (22 genomes, 148 targets, six predictors) | RB-TnSeq fitness-validated exact EC | development / training; documentation of correction behaviour | untouched test; negatives (one validated EC per target) |
| **P03 panels** (13 and 21 genomes) | GOA experimental GO; cluster-resident | examined test or training | untouched test (v2 tuning saw them); most genomes carry tiny experimental truth (e.g. vcholerae 46 bytes, rprowazekii 25) |
| **`cafa_temporal_holdout_real.json`** | aggregate F-max per aspect (4,712 MF / 5,969 CC / 5,444 BP proteins) | nothing | no per-protein rows, undefined t0 |
| **`full_genome_eval_gbff_*.pkl`** (12 files) | 73 genomes, 52,814 proteins, per-protein F-max and 54,642 score changes from the old runner | nothing | not raw candidates; truth source unrecorded; METEOR arm is the Python re-implementation |
| **suite2 multikingdom 3,100** | taxonomy + counts for 3,100 taxa | taxonomy-aware split design | 2,268 rows are eggNOG *predictions*; the 832 GOA rows are partial archaeal proteomes (annotated proteins only, `n_proteins == n_annotated`) |
| **Swiss-Prot flat file** | 573,661 entries, 14,803 organisms | label-rich organism selection | whole-genome context (reviewed proteins only) |

Organisms with enough Swiss-Prot experimental labels to serve as a labelled genome (reviewed proteins only, so partial proteomes):

| Threshold (proteins per organism) | experimental GO | EC with ECO:0000269 | EC, any evidence |
| --- | --- | --- | --- |
| >= 25 | 88 (19 bacteria, 4 archaea, 56 eukaryotes) | 58 (15 bacteria) | 1,220 |
| >= 50 | 49 (8 bacteria) | 29 (9 bacteria) | 1,006 |
| >= 100 | 31 (6 bacteria, 25 eukaryotes) | 15 (4 bacteria) | 868 |

Experimentally labelled whole genomes in distinct families are therefore on the order of tens, not thousands. Swiss-Prot EC numbers carry no experimental guarantee: 868 organisms reach 100 EC-bearing entries but only 15 reach 100 with an experimental ECO tag.

## 2. Findings that change the protocol

1. **GOA record dates are not first-assertion dates.** Of 514,240 (protein, term) pairs that are already in the pre-t0 `train_terms.tsv` snapshot, 52,406 carry a record date on or after 2024-01-01 and 24,628 on or after 2025-06-01 (median 2014-05, but 95th percentile 2025-04 and maximum 2026-06-13: records are re-dated on re-curation/re-release). Conversely 14,035 of the 38,216 pairs absent from the snapshot are dated before 2024. A date cut therefore mislabels known annotations as "new" and old ones as "prospective". Temporal claims must come from the difference of two dated snapshots, and even then only the pairs absent from the earlier snapshot, with t0 and release recorded.
2. **The SwissProt entry-date split was not an annotation-time split.** `make_temporal_split.py` used the `DT ... integrated into UniProtKB/Swiss-Prot` date of the *entry*; the file carries no per-annotation date. Its docstring claimed a per-annotation first-seen date and called the side files "held-out". The "train" side contained annotations added to old entries after the cutoff. Also: the EC extractor required `DE   EC=` but real files indent `DE            EC=...`, so the script emitted **zero EC rows** (0 of 13,259 EC lines in the first 3M lines matched the guard); undated entries silently went to the training side; `3.5.1.n2` was truncated. Fixed: outputs are `*.integrated_le_cutoff.tsv` / `*.integrated_gt_cutoff.tsv` with an evidence and an integration-date column and a `.provenance.json` stating `is_annotation_time_split: false`; undated entries are excluded and counted; only complete four-level ECs are kept; `--evidence experimental` is available. No other file in the repository references the script or its old output names.
3. **The "validation" fold is a subset of training** (25,000 of 25,000 in the 82,404). It is a CV fold, not an independent split.
4. **81,556 of 88,894 annotated test-superset proteins are training proteins.** Only 7,338 proteins (10,043 pairs) are not in the training set; 38,216 pairs on training proteins are absent from the snapshot. Evaluating on the superset without removing known pairs is leakage.
5. **Price148 stays development.** The registry forbids it in `test`; `audit-splits` and `build-table` refuse such a table.

## 3. Label policy (implemented, tested)

* Positive = an allowed-evidence assertion naming the same protein and function. Default evidence = the CAFA experimental set; IEA or non-experimental evidence is dropped and counted (`label_records_dropped`) unless the dataset policy lists it explicitly.
* GO is exact-term unless `--go-obo` is given, then propagated over `is_a` and `part_of` (roots dropped, alt ids mapped); the choice is written to the manifest.
* **Missing labels are not negatives.** A row with no positive is excluded (`no_truth_not_negative` for proteins with no truth of that type, `unconfirmed_not_negative` for labelled proteins). A 0 label exists only where the dataset policy declares `negative_policy: closed_world_declared` with a written justification, only on proteins with at least one assertion of that type, and every such row carries `label_basis: closed_world_unconfirmed`.
* **Roles are not guessed.** EC truth labels only `CATALYTIC` candidates, GO truth only `UNSPECIFIED` / `COMPLEX_SUBUNIT`; accessory, transporter and regulator candidates are excluded (`role_not_labelled`), never negatives.
* **No EC exclusivity, no role/null negatives.** `locus_group_id` is set only from an explicit locus-group file whose groups declare complete truth *and* a policy carrying `locus_truth_complete`; otherwise the trainer correctly reports the null layer as untrained.
* `truth_recall_ceiling` reports, per genome and type, how many truth pairs have any candidate row at all.

## 4. Executable preparation path

```bash
PY=/tmp/gspa-v3-testenv/bin/python
# 0. profile the local cohorts (read-only, ~1 min)
$PY -I benchmark/v3/audit_training_cohorts.py --recon-dir <cafa6_recon> --kaggle-dir <cafa6/kaggle-official> \
    --suite1 benchmark/manifests/suite1_cafa_temporal_manifest.tsv --holdout-sample <swissprot_holdout_sample.tsv> \
    --swissprot-dat <uniprot_sprot.dat.gz> --out benchmark/v3/results/training-data-audit/cohort_profile.json
# 1. raw candidates: the production engine writes one export per genome (no labels)
gspa integrate --batch-manifest M.tsv ... --feature-export-dir features/      # features/<genome>/candidate_features.jsonl
# 2. author three small inputs: splits.tsv (genome_id split group source [species]),
#    labels.tsv (genome_id protein_id function_id evidence), label_policy.json
# 3. audit the split on its own (exit 2 on any violation)
$PY benchmark/v3/prepare_training_data.py audit-splits --splits splits.tsv --proteins proteins.tsv --protein-clusters clusters.tsv
# 4. build the trainer input + manifest of exact counts (refuses if the audit fails)
$PY benchmark/v3/prepare_training_data.py build-table --feature-dir features/ --contract features/<g>/feature_contract.json \
    --labels labels.tsv --splits splits.tsv --label-policy label_policy.json --go-obo go.obo --out-dir table/
# 5. train (frozen before test; see docs/concepts/learned-context.md)
$PY scripts/train_context_model.py train --table table/training_table.jsonl --contract <contract> --out-dir model/
```

`label_policy.json` shape: `{"datasets": {"<source>": {"EC": {"allowed_evidence": "experimental", "negative_policy": "none"}, "GO": {...}}}, "locus_truth_complete": "<justification, optional>"}`. Every genome's `source` must be a registry id. `build-table` verifies that each export's rows hash matches its `export_manifest.json`, `labels_included` is false, and the contract fingerprint matches; it refuses exporter rows that already carry a label/split/group.

## 5. What an independent split claim can and cannot say

Supportable once a real table passes the audit:

* train, calibration and test are disjoint by genome and by declared family/species group (enforced by `audit-splits` and again by the trainer);
* no registered development source is in the test split; protein ids do not overlap across splits; if a homology map was supplied it covers every protein and no cluster spans splits (the report records whether this was checked);
* labels are traceable per row to `(dataset, genome, protein, function, evidence)` and to an input hash; every 0 is stamped with its basis.

Not supportable from what exists today:

* "untouched test": Price148 and the P03 panels were examined during v2; `examined` sources in test give a disjointness claim only (`untouched_test_claim_supportable` is false);
* "prospective", "annotation-time" or "future" anything from GOA record dates or Swiss-Prot entry dates;
* "no pretraining leakage": predictor channels (ESM/ProstT5/CLEAN/DeepEC) may have seen test proteins; the report always states `pretrained_model_leakage: not_assessed`;
* genus-held-out claims without a taxonomy-keyed `group` and `--disjoint genome_id,group,species`;
* any context-effect claim from CAFA6-type data (no context).

## 6. Remaining tasks

1. **Choose the genome-level cohort.** The only genome-level labelled data are Price148 (development), the P03 panels (examined, mostly tiny truth) and label-rich Swiss-Prot organisms (tens). Recommended: train/calibrate on Price148 + P03 genomes grouped by genus, and build a **new** held-out set of 15-30 label-rich organisms (>= 50 experimental GO or ECO:0000269 EC proteins, distinct genera/families from every training genome) whose reference proteomes are fetched with exact accessions and run through `gspa integrate --feature-export-dir`. These are partial-truth, not whole-genome-complete, so metrics stay on labelled proteins.
2. Generate raw claims for those proteomes (root owns remote runs; this task did not start any).
3. Homology/family clusters: none ships with any cohort. Build with a pinned tool (e.g. MMseqs2 at a declared identity) over all proteins across splits, then pass `--protein-clusters`.
4. Taxonomy `group` and `species` for each genome (suite2 metadata carries phylum..species for 3,100 taxa).
5. Declare t0 and GOA/UniProt releases for any temporal claim and use snapshot differences (two dated releases) rather than dates.
6. Decide per dataset whether a closed-world policy is defensible (Price has one validated EC per target but promiscuous enzymes exist; keep `none` until the root decides).
7. Locus groups: only from masked-reaction experiments with the full candidate list, not the truncated top-5 in `lro50_real_unimatrix.json`.
8. CAFA6 cohort: if used for the sequence block, remove known pairs, drop the 25,000-protein "validation" fold as a split, and report it as sequence-only.
