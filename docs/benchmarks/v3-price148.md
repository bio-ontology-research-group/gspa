# Price-148 production inputs (`benchmark/v3/prepare_price.py`)

Exports the inputs of the Price RB-TnSeq enzyme-identity benchmark for the v3 runner. It does not run `gspa` and computes no metric; it only converts a frozen predictor artifact into `gspa integrate` claims, labels and runner manifests, with a report that discloses what the data cannot support. Tested with `/tmp/gspa-v3-testenv/bin/python` (numpy only is needed; pandas is needed only to read original DataFrame pickles).

```bash
/tmp/gspa-v3-testenv/bin/python benchmark/v3/prepare_price.py export \
    --artifact benchmark/manifests/reference_data/ibex_extracted/price149_all6_compact.pkl.gz \
    --manifest benchmark/manifests/suite6_price149_manifest.tsv --out /tmp/gspa-v3-price-inputs
```

Generated output is not committed. Options: `--variants` (default: the six single predictors and `ensemble_pool4`; `ensemble_all6` available), `--namespace auto|always|never`, `--include-ec-top`, `--raw-matrix-map JSON --trust-pickle` (full matrices).

## What the inputs are

| Item | Source and handling |
|---|---|
| Artifact | `price149_all6_compact.pkl.gz`, extraction `v2-2026-10-19` (ibex job 53987410). Six predictors (CLEAN, DeepProZyme, EnzBERT, MAPred, GraphEC, TopEC) × 22 assemblies, 148 targets. Hashed in the report |
| Targets | `target_scores`: every EC with score ≥ 0.01 plus the 100 highest, kept **unmodified and unthresholded** (weak candidates retained). CLEAN scores for some targets peak at ~6e-5; they are exported as raw model scores, not rescaled, and carry no `calibrated_prob`. The engine applies its default sigmoid(4r−2) to all sources, so CLEAN posteriors stay low. A learned calibration is future work and was not fitted on Price |
| Context | `neighbor_ecs` (top-5 ECs ≥ 0.15 of the ±10 neighbouring proteins) exported as claims for those neighbours, flagged `origin: neighbor_ecs`; `--include-ec-top` adds the 8 best proteins per EC (`ec_top_scores`). Precedence per (protein, EC): target_scores > ec_top_scores > neighbor_ecs, so no pair is claimed twice by one source |
| Evidence type | The engine knows only `clean`; the other five sources would be dropped silently, so `evidence_type` is written explicitly: sequence deep learning for CLEAN, DeepProZyme, EnzBERT, MAPred; structure deep learning for GraphEC, TopEC. Fixed before any Price result was inspected; it selects the combiner's reliability class |
| Labels | Manifest `ground_truth_ecs` → `labels/<assembly>.ec.tsv` (`protein_id`, `ec`) and `labels/price_label_audit.tsv`. Complete four-level ECs only; multi-EC rows keep every complete EC |

## Identifiers

Claims and labels use the **canonical protein id**: `raw_target_stats[entry].matched_id`, the proteome id the predictors scored (not the Price/RefSeq entry id, which differs for 30 rows). DeepProZyme drops version suffixes for 19 `WP_` targets (`WP_063460136` vs `WP_063460136.1`); within an assembly a bare accession is mapped to its versioned form, in targets and in neighbour/EC-top ids, and the rewrites are counted in the report. An entry whose predictors disagree on the protein after that is *not* exported or labelled under a guess (`alias_conflicts`). The canonical id was checked against the manifest's `proteome_protein_id` (agreement for all 148). If a canonical id occurs in two assemblies, `--namespace auto` prefixes every id with `<assembly>__`; claims and labels are namespaced together. All outputs are sorted and serialised deterministically (re-export is byte-identical, independent of input order).

## Exclusions (reported with reasons)

| Case | Rule |
|---|---|
| `NP_384885` (GCA_000006965.1) | `proteome_protein_id` is the replicon `NC_003047.1`, not a protein, and no predictor scored it: excluded, 148 of 149 rows remain |
| Partial EC in truth (`1.1.1.-`, fewer than four levels) | dropped from the label set and listed under `partial_ec_drops`; a row with no complete EC is excluded |
| Predictors disagree on the matched protein | excluded from labels and claims, listed |
| Assembly or entry absent from the artifact | excluded, listed |

## Truth versus inference

`build_claims` takes the artifact only; `build_labels` takes the manifest and the alias map only. Truth is used for writing labels and for the accession-mapping audit (`alias_agrees_with_manifest`), never for claim content, retention, priors, pathway conditions or manifests' engine arguments. Metabolic strata, seed reaction ids and universal-model flags (all truth-derived) are not read into the claims. A test exports with entirely different truth columns and requires byte-identical claim files.

Limitation that cannot be removed here: the artifact was extracted for the Price targets, so the claim files contain the target proteins and their neighbourhoods only. Every target is an enzyme with a known EC and its presence is a property of the design. Accuracy is measured among these known enzymes (top-1 among their own candidate ECs), not as enzyme detection.

## Truncation, coordinates and full context

The compact artifact is **not the whole proteome** (6,233 proteins in the first genome; 75 carry CLEAN claims there) and holds **no gene coordinates**. None are invented from gene order, so no `genome_layout` is exported and layout-dependent arms (RLC, genomic-context priors beyond what claims alone carry) are not in the manifests. `export_report.json` states this under `truncation` and lists per predictor and genome the proteome size, claim proteins, context proteins, claims by origin, duplicate merges, score mismatches between sources, non-complete EC claims and the number of targets with tied top scores.

The original full-matrix schema is supported: `--raw-matrix-map` takes `{predictor: {assembly: path}}` with `.npz` (`index`, `columns`, `values`) or the original pandas pickle (`--trust-pickle`; unpickling executes code). Rows are proteins, columns ECs with optional `EC:` prefix. The same retention rule (≥ 0.01 or top 100 per protein) yields whole-proteome claims, `proteome_complete: true`, and the manifests drop the truncation notice. Exact raw files, with sha256 and size, are in `raw_sources_needed_for_full_context` (ibex, read-only, owner `kexin`); the layout/coordinates must come from the genome annotation (GFF), which is not in the artifact. The METEOR archive paths are listed too, but METEOR is not used.

## Arms and matched evidence

One manifest per predictor and per ensemble (`manifests/v3_price_<variant>.json`); inside a manifest every arm reads the identical claims file, so evidence is matched across arms by construction: `stage0` (`--enable-priors none`, the baseline; with one claim per pair its posterior is a monotone function of the predictor's raw score) and `stage1` (default priors with the KEGG reference files). `ensemble_pool4` is exactly the union of the four single-predictor claim sets (CLEAN, DeepProZyme, MAPred, EnzBERT; each (protein, EC, source) once, enforced, so the baseline is never added twice). Single predictors are compared by running their manifests separately; per-arm claims inside one runner manifest are deliberately not supported. Arms that require a genome layout or whole-proteome context are added only after full-context export. Thresholds are not tuned on Price: the fixed 0.5 threshold is a grid constant, and with raw CLEAN/GraphEC scores and the default sigmoid calibration, top-1 accuracy is the meaningful metric.

## Status of the generated export (local run, `/tmp/gspa-v3-price-inputs`)

148 labelled targets in 22 genomes, 151 label pairs, 1 exclusion. Claims: CLEAN 15,676, DeepProZyme 20,451, EnzBERT 16,260, MAPred 8,322, GraphEC 358,315 (75 MB), TopEC 2,352, `ensemble_pool4` 60,709. `gspa integrate` (Gradle launch, `--lite`) accepted the exported GCF_001633105.1 claims for the pooled and TopEC manifests: every claimed pair reached the output, with up to four supporting sources per pair. No benchmark has been run on these inputs. The predictors' training sets may contain Price enzymes; this is not assessed.
