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

Organisms with enough Swiss-Prot labels to serve as a labelled genome (reviewed proteins only, so partial proteomes). Evidence tiers are never merged: *strict* = GO EXP/IDA/IPI/IMP/IGI/IEP/H* or EC tagged ECO:0000269; *CAFA-eligible* additionally admits the curator-asserted, non-experimental TAS/IC and ECO:0000303/305 (`benchmark/v3/evidence_codes.py`). `labelled_strict` counts a protein once if it has strict GO or strict EC. Totals include eukaryotes; bacteria/archaea in parentheses.

| Threshold (proteins per organism) | strict GO | strict GO or EC (`labelled_strict`) | strict EC (ECO:0000269) | CAFA-eligible GO | EC, any evidence |
| --- | --- | --- | --- | --- | --- |
| >= 25 | 83 (19 bacteria, 4 archaea) | 123 (33 bacteria, 13 archaea) | 58 (15 bacteria, 7 archaea) | 88 | 1,220 |
| >= 50 | 45 (8 bacteria, 1 archaea) | 59 (14 bacteria, 3 archaea) | 29 (9 bacteria, 1 archaea) | 49 | 1,006 |
| >= 100 | 29 (6 bacteria) | 31 (7 bacteria) | 15 (4 bacteria) | 31 | 868 |

The earlier version of this table called the CAFA-eligible GO column "experimental"; it overstated the strict count (88 vs 83 at >= 25). Experimentally labelled bacterial/archaeal genomes are on the order of tens, not thousands, and each carries tens of labelled proteins. Swiss-Prot EC numbers carry no experimental guarantee: 868 organisms reach 100 EC-bearing entries but only 15 reach 100 with an experimental ECO tag.

## 2. Findings that change the protocol

1. **GOA record dates are not first-assertion dates.** Of 514,240 (protein, term) pairs that are already in the pre-t0 `train_terms.tsv` snapshot, 52,406 carry a record date on or after 2024-01-01 and 24,628 on or after 2025-06-01 (median 2014-05, but 95th percentile 2025-04 and maximum 2026-06-13: records are re-dated on re-curation/re-release). Conversely 14,035 of the 38,216 pairs absent from the snapshot are dated before 2024. A date cut therefore mislabels known annotations as "new" and old ones as "prospective". Temporal claims must come from the difference of two dated snapshots, and even then only the pairs absent from the earlier snapshot, with t0 and release recorded.
2. **The SwissProt entry-date split was not an annotation-time split.** `make_temporal_split.py` used the `DT ... integrated into UniProtKB/Swiss-Prot` date of the *entry*; the file carries no per-annotation date. Its docstring claimed a per-annotation first-seen date and called the side files "held-out". The "train" side contained annotations added to old entries after the cutoff. Also: the EC extractor required `DE   EC=` but real files indent `DE            EC=...`, so the script emitted **zero EC rows** (0 of 13,259 EC lines in the first 3M lines matched the guard); undated entries silently went to the training side; `3.5.1.n2` was truncated. Fixed: outputs are `*.integrated_le_cutoff.tsv` / `*.integrated_gt_cutoff.tsv` with an evidence and an integration-date column and a `.provenance.json` stating `is_annotation_time_split: false`; undated entries are excluded and counted; only complete four-level ECs are kept; `--evidence experimental` is available. No other file in the repository references the script or its old output names.
3. **The "validation" fold is a subset of training** (25,000 of 25,000 in the 82,404). It is a CV fold, not an independent split.
4. **81,556 of 88,894 annotated test-superset proteins are training proteins.** Only 7,338 proteins (10,043 pairs) are not in the training set; 38,216 pairs on training proteins are absent from the snapshot. Evaluating on the superset without removing known pairs is leakage.
5. **Price148 stays development.** The registry forbids it in `test`; `audit-splits` and `build-table` refuse such a table.

## 3. Label policy (implemented, tested)

* Positive = an allowed-evidence assertion naming the same protein and function. **Evidence categories are separate** (`benchmark/v3/evidence_codes.py`): `experimental` (the default, strict: EXP/IDA/IPI/IMP/IGI/IEP/H*, ECO:0000269), `cafa_eligible` (adds the curator-asserted TAS/IC and ECO:0000303/305, which are *not* experimental), `all`, or an explicit list. TAS/IC dropped by the strict default are counted (`evidence_curated_inferential_not_experimental`), and every positive row records `label_evidence` and `label_evidence_categories`. The same constants drive `make_temporal_split.py` (`--evidence experimental|cafa_eligible|all`; EC rows now carry their ECO tag, so the filter applies to EC too) and the cohort profile.
* GO is exact-term unless `--go-obo` is given, then propagated over `is_a` and `part_of` (roots dropped, alt ids mapped to their primary id, obsolete terms dropped and counted as `obsolete_go_term` / `obsolete_go_candidate`). The OBO is read with the existing `v3_metrics.load_obo`: the former private parser never reset its stanza state, so `[Typedef]` stanzas and `is_obsolete` flags could attach to the wrong term. The OBO hash and propagation mode are written to the manifest.
* **Missing labels are not negatives.** A row with no positive is excluded (`no_truth_not_negative` for proteins with no truth of that type, `unconfirmed_not_negative` for labelled proteins). A 0 label exists only where the dataset policy declares `negative_policy: closed_world_declared` with a written justification, only on proteins with at least one assertion of that type, and every such row carries `label_basis: closed_world_unconfirmed`.
* **Roles are not guessed.** EC truth labels only `CATALYTIC` candidates, GO truth only `UNSPECIFIED`. A `COMPLEX_SUBUNIT` (or any other) role is labelled only from role truth that says so: a labels-file `role` column whose role the dataset policy declares under `role_truth` *with its source*, or - when the policy opts in with `complex_membership_from_cc` (and a source) and a GO OBO is given - an allowed-evidence cellular-component assertion to a protein-containing-complex term (GO:0032991 descendant) for a candidate naming that same complex term. An unrelated MF/BP label never labels a subunit. Every positive records `label_role_basis`; roles without truth are excluded (`role_not_labelled`, `complex_subunit_without_role_truth`), never negatives.
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
# 3. audit the split on its own (exit 2 on any violation); proteins.tsv: genome_id protein_id [kind = target|context]
$PY benchmark/v3/prepare_training_data.py audit-splits --splits splits.tsv --proteins proteins.tsv --protein-clusters clusters.tsv \
    --disjoint genome_id,group,genus     # taxonomy keys (species/genus/family) need that column for EVERY genome
# 4. build the trainer input + manifest of exact counts (refuses if the audit fails)
$PY benchmark/v3/prepare_training_data.py build-table --feature-dir features/ --contract features/<g>/feature_contract.json \
    --labels labels.tsv --splits splits.tsv --label-policy label_policy.json --go-obo go.obo \
    --example-clusters example_clusters.tsv --disjoint genome_id,group,genus --out-dir table/
# 5. train (frozen before test; see docs/concepts/learned-context.md)
$PY scripts/train_context_model.py train --table table/training_table.jsonl --contract <contract> --out-dir model/
```

`label_policy.json` shape: `{"datasets": {"<source>": {"EC": {"allowed_evidence": "experimental", "negative_policy": "none"}, "GO": {...}}}, "locus_truth_complete": "<justification, optional>"}`. Every genome's `source` must be a registry id. `build-table` validates every `<genome>/` export against the production writer (`FeatureExporter`): manifest schema `gspa.context-feature-export/1` with its `rows`/`rows_file`/`contract_file` entries, `rows_sha256`, `labels_included: false`, the shipped `feature_contract.json` equal to the supplied contract (id, fingerprint, every name and group), row schema `gspa.context-feature-row/1`, per-row contract id and fingerprint, the exact `features` column order, finite values, the `available_contexts` keys, `candidate_key = function_key|ROLE`, `genome_id` equal to its directory, and no label/split/group/example_id. Any drift (a stale schema, a missing key) stops the build. `FeatureExportBuilderSpec` (Groovy) runs the real `gspa integrate --feature-export-dir` route and feeds the output through this builder; the Python tests use a fixture with the identical layout, and `GSPA_REAL_EXPORT_DIR` points the Python suite at any real export directory.

**Context proteins vs labelled examples.** Whole-genome context (all proteins of a genome, used only for features) cannot be made cluster-disjoint across genomes: conserved core families occur in every genome. The audit therefore treats `kind=context` proteins as informational (`context_sharing`) and enforces disjointness only for labelled `kind=target` proteins. With `--example-clusters` (`protein_id[, genome_id], cluster`, covering every labelled protein) `build-table` masks labelled examples of a lower-precedence split whose family also occurs among higher-precedence labelled examples (train against calibration/test, calibration against test; test examples are never dropped). Every masked example is written to `masked_examples.tsv` (example, split, label, cluster, clashing split, reason) and counted per split, with positives, in `family_masking`; hard test cases whose family occurs in train are kept and counted; an unmapped labelled protein, or a split left without positives, stops the build; the kept table is re-verified family-disjoint. A genus-disjoint main split uses `--disjoint genome_id,group,genus`; a family-disjoint evaluation adds `family`. Without `--example-clusters` the manifest states that masking was not applied.

## 5. What an independent split claim can and cannot say

`audit-splits` reports `checked_axes` and `unchecked_axes`. An independence claim is supportable only when genome id, group, source registry, protein-id overlap and homology (labelled targets) were all checked and clean (`independent_test_claim_supportable`, with a generated `independent_test_claim_text` naming what was and was not checked). A split table alone is never enough: without a protein list and a cluster map those axes are listed as unchecked.

Supportable once a real table passes the audit with all axes checked:

* train, calibration and test are disjoint by genome, by declared group, and (when requested and present for every genome) by species/genus/family. Duplicate genome rows with conflicting group/source/taxonomy are refused even inside one split; requested taxonomy that is missing is refused, not assumed;
* no registered development source is in the test split; labelled target protein ids and families do not cross splits (context proteins may share homologs and are only reported);
* labels are traceable per row to `(dataset, genome, protein, function, evidence, role basis)` and to an input hash; every 0 is stamped with its basis.

Not supportable from what exists today:

* "untouched test": the registry records only *known* past use. The absence of a recorded use is not evidence that data are untouched, so the report always carries `untouched_test_claim: not_assessed`; Price148 and the P03 panels are known to have been examined;
* "prospective", "annotation-time" or "future" anything from GOA record dates or Swiss-Prot entry dates;
* "no pretraining leakage": predictor channels (ESM/ProstT5/CLEAN/DeepEC) may have seen test proteins; the report always states `pretrained_model_leakage: not_assessed`;
* any context-effect claim from CAFA6-type data (no context).

## 6. Selected label-rich organisms (selection manifest; no genomes acquired)

`benchmark/v3/select_label_rich_genomes.py` reads the per-organism table written by `audit_training_cohorts.py --swissprot-dat <snapshot> --taxa-tsv` and writes `benchmark/v3/results/label-rich-selection/{selection_manifest.json, selected_organisms.tsv}`. The input is the local Swiss-Prot flat file `uniprot_sprot.dat.gz` (573,661 entries, 687,123,790 bytes, sha256 `a6907c4dc83ff62253fcd94c5aed342b254295fde57b7fbe76f2c21d0231cbdc`, recorded in `swissprot_profile.json` next to the table). No predictor output is an input.

Criteria (all recorded in the manifest): Bacteria or Archaea; at least 20 proteins with strict experimental GO or ECO:0000269 EC (`labelled_strict`; TAS/IC never count); one organism per species and per genus; one UniProt proteome that covers at least half of the organism's Swiss-Prot entries (otherwise the labels cannot be tied to one exact genome); and the **whole genus** of every Price148, P03 panel13/panel21 and benchmark-panel organism excluded (`registry/evaluation_exclusions.json`, built from `benchmark/panel/panel_manifest.tsv`, `suite2_proteomes25_manifest.tsv` and the panel13 tags, with synonyms such as Escherichia/Shigella and Mycobacterium/Mycolicibacterium). Every exclusion is listed with its reason.

Result: **19 organisms (9 bacteria, 10 archaea) in 19 genera and 17 families, 577 strictly labelled proteins in total (20-58 per organism).** 1,528 bacterial/archaeal organisms fall below the label threshold; of the candidates above it, 48 were excluded (23 genus overlaps with Price/P03/panels, 19 exact panel taxids, 4 without a dominant proteome, 2 genus cap). The proposed split assigns whole families (a hash-ordered proposal, not frozen; re-run `audit-splits` after acquisition): 11 train, 3 calibration, 5 test organisms.

| Split | Organisms (NCBI taxid; strict Swiss-Prot labels) |
| --- | --- |
| test | *Methanosarcina mazei* 192952 (27); *Clostridium acetobutylicum* 272562 (25); *Paramagnetospirillum magneticum* 342108 (24); *Dickeya dadantii* 198628 (22); *Magnetospirillum gryphiswaldense* 431944 (22) |
| calibration | *Haemophilus influenzae* 71421 (36); *Campylobacter jejuni* 192222 (31); *Methanothermobacter marburgensis* 79929 (29) |
| train | *Saccharolobus solfataricus* 273057 (58); *Thermococcus kodakarensis* 69014 (49); *Haloferax volcanii* 309800 (44); *Sulfolobus acidocaldarius* 330779 (39); *Synechococcus elongatus* 1140 (35); *Methanococcus maripaludis* 267377 (25); *Sulfurisphaera tokodaii* 273063 (25); *Acinetobacter baylyi* 62977 (24); *Aeropyrum pernix* 272557 (21); *Thermoplasma acidophilum* 273075 (21); *Legionella pneumophila* 272624 (20) |

Caveats that stay attached to this list:

* It is a *selection*, not a dataset: no genome, feature export or labels table exists yet (registry entry `swissprot_label_rich_v1` records `genomes_acquired: 0`).
* Each organism is a partial proteome with tens of labels, so metrics stay on labelled proteins, and 577 proteins are far too few to fit a stacker; this is an evaluation/calibration reservoir. GOA experimental annotations of the same UniProt proteomes can add labels at acquisition time and must be counted separately.
* Swiss-Prot labels are inputs of the pretrained predictors, so leakage is `not_assessed`.
* **Assembly mapping is `unresolved` for every organism.** The manifest carries the dominant UniProt proteome id (for example UP000001974 for *S. solfataricus*), but the exact NCBI assembly and its protein FASTA and coordinates have not been fetched or matched, and no stand-in genome is used.
* The 22 Price host assemblies carry no organism names in the repository, so `price_genus_check` is `pending` for every selected genus until they are resolved against NCBI.
* *Thermococcus kodakarensis* shares a family with the excluded *Pyrococcus* panel genomes (`family_shared_with_excluded_panel_organism`).

## 7. Remaining tasks

1. **Acquire the 19 genomes** (root owns remote jobs): resolve each dominant UniProt proteome to its exact NCBI assembly, fetch the full protein FASTA and coordinates, record unresolved mappings instead of guessing, and resolve the 22 Price host assemblies to genera to clear `price_genus_check`.
2. Generate raw claims and `--feature-export-dir` exports for those proteomes (nothing was started here).
3. Homology/family clusters: none ships with any cohort. Build them with a pinned tool (for example MMseqs2 at a declared identity) over the labelled proteins (`--example-clusters`) and, for the audit, over all proteins with `kind=context|target`.
4. Taxonomy columns (`genus`, `family`) in `splits.tsv` for every genome; the selection manifest provides them for the 19 organisms, Price and P03 genomes need theirs.
5. Declare t0 and GOA/UniProt releases for any temporal claim and use snapshot differences (two dated releases) rather than dates.
6. Decide per dataset whether a closed-world policy is defensible (Price has one validated EC per target but promiscuous enzymes exist; keep `none` until the root decides).
7. Locus groups: only from masked-reaction experiments with the full candidate list, not the truncated top-5 in `lro50_real_unimatrix.json`.
8. CAFA6 cohort: if used for the sequence block, remove known pairs, drop the 25,000-protein "validation" fold as a split, and report it as sequence-only.

Training is not complete and no calibrated model exists.
