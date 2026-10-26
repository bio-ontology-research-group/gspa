# Price-148 production inputs (`benchmark/v3/prepare_price.py`)

Exports the inputs of the Price RB-TnSeq enzyme-identity benchmark for the v3 runner. It does not run `gspa` and computes no accuracy metric; it converts a frozen predictor artifact into `gspa integrate` claims, labels and runner manifests, with a report that discloses what the data cannot support. Tested with `/tmp/gspa-v3-testenv/bin/python` (numpy only; pandas is needed only to read original DataFrame pickles).

```bash
/tmp/gspa-v3-testenv/bin/python benchmark/v3/prepare_price.py export \
    --artifact benchmark/manifests/reference_data/ibex_extracted/price149_all6_compact.pkl.gz \
    --manifest benchmark/manifests/suite6_price149_manifest.tsv --out /tmp/gspa-v3-price-full-inputs
/tmp/gspa-v3-testenv/bin/python benchmark/v3/prepare_price.py check-roundtrip --claims C.jsonl --integrated OUT.tsv
```

Generated output is not committed. Options: `--variants` (default: six single predictors and `ensemble_pool4`; `ensemble_all6` available), `--namespace auto|always|never`, `--include-ec-top`, `--no-target-full` (truncated candidates only; not recommended), `--evidence-types JSON`, and for full matrices `--raw-matrix-map JSON --trust-pickle --retention all_nonzero|threshold_topk --min-score --top-k`.

## What the inputs are

| Item | Source and handling |
|---|---|
| Artifact | `price149_all6_compact.pkl.gz`, extraction `v2-2026-10-19` (ibex job 53987410), sha256 `03eb5c51…ebe48`. Six predictors (CLEAN, DeepProZyme, EnzBERT, MAPred, GraphEC, TopEC) × 22 assemblies, 148 targets |
| Target candidates | **Every nonzero entry of `target_full`** (float32, targets × `_meta.ec_vocab`), not the artifact's truncated `target_scores` (≥ 0.01 plus the 100 highest). Verified before use: matrix shape equals (targets, vocabulary), row order equals `target_ids`, vocabulary is unique, every value is finite and in [0, 1], and every `target_scores` entry equals the matrix value within float32 rounding (rel 1e-6). Any disagreement (altered value, row order, vocabulary, a zeroed candidate, an EC outside the vocabulary) **rejects the export** instead of falling back. If `target_full` is absent or structurally unusable the exporter falls back to `target_scores` and records `truncated: true` with the reason; coverage is then reported as not full |
| Why it matters | The truncated lists dropped most candidates for the dense predictors: CLEAN 760,154 of 774,954 nonzero target candidates (98 %), DeepProZyme 386,267 of 401,820, GraphEC 396,951 of 750,508 (as little as 13 % of the score mass kept in some genomes); EnzBERT, MAPred and TopEC lose none. Coverage is computed from the artifact alone (counts and score mass); no retention choice looks at labels |
| Vocabulary | 7,433 entries: 7,021 complete ECs, 41 preliminary ECs (`1.1.1.n11`) and 371 hierarchical classes (`1`, `1.1`, `1.1.1`). Preliminary ECs are real predictor outputs and are kept (they can never be a complete-EC top-1). Hierarchical classes are never candidates; no nonzero score was found on them, and any would be counted in `dropped_hierarchical_partial_nonzero`. `ec_vocab` is never printed or copied into reports |
| Calibration | Raw model scores exported unchanged; no `calibrated_prob`. The engine's default calibrator is **provisional**: `sigmoid(4·raw − 2)` for every source. It is neither the identity nor learned: raw 0 maps to 0.12, small scores are compressed (CLEAN peaks at ~6e-5 for some targets, so their posteriors sit near 0.12). A learned calibration is future work and was not fitted on Price |
| Context | `neighbor_ecs` (top-5 ECs ≥ 0.15 of the ±10 neighbouring proteins) exported as claims for those neighbours, `origin: neighbor_ecs`; `--include-ec-top` adds the 8 best proteins per EC. Precedence per (protein, EC): target > ec_top_scores > neighbor_ecs, so no pair is claimed twice by one source; shared neighbours merge |
| Evidence type | The engine knows only `clean` as a source name; the other five would be dropped silently, so `evidence_type` is written explicitly. **Provisional architecture assumption, not verified against the models' sources** (the repository has no model cards; no network was available): CLEAN, DeepProZyme, EnzBERT, MAPred are treated as sequence-based, GraphEC and TopEC as structure-based. The assignment was chosen without using the truth columns; it was **not** fixed before Price was looked at (earlier audits examined Price results), and it selects the combiner's reliability class. Override with `--evidence-types` and re-export. Whether any endpoint's training data overlaps the Price enzymes is unknown |
| Labels | Manifest `ground_truth_ecs` → `labels/<assembly>.ec.tsv` (`protein_id`, `ec`) and `labels/price_label_audit.tsv`. Complete four-level ECs only; multi-EC rows keep every complete EC |

This is a development reuse of a dataset that earlier audits already examined, not a blind test. No threshold, weight or calibration was tuned on it.

## Identifiers

Claims and labels use the **canonical protein id**: `raw_target_stats[entry].matched_id`, the proteome id the predictors scored (the Price entry id differs for 30 rows). DeepProZyme drops version suffixes for 19 `WP_` targets; within an assembly a bare accession is mapped to its versioned form (targets, neighbour and EC-top ids), and the rewrites are counted. An entry whose predictors still disagree is not exported or labelled under a guess. The canonical id agrees with the manifest's `proteome_protein_id` for all 148. If a canonical id occurs in two assemblies, `--namespace auto` prefixes every id with `<assembly>__` (claims and labels together). Outputs are sorted and serialised deterministically; re-export is byte-identical, independent of input order.

## Exclusions (reported with reasons)

| Case | Rule |
|---|---|
| `NP_384885` (GCA_000006965.1) | `proteome_protein_id` is the replicon `NC_003047.1`, not a protein, and no predictor scored it: excluded, 148 of 149 rows remain |
| Partial EC in truth (`1.1.1.-`, fewer than four levels) | dropped from the label set, listed under `partial_ec_drops`; a row with no complete EC is excluded |
| Predictors disagree on the matched protein | excluded from labels and claims, listed |
| Assembly or entry absent from the artifact | excluded, listed |

## Truth versus inference

`build_claims` takes the artifact only; `build_labels` takes the manifest and the alias map only. Truth is used for writing labels and for the accession-mapping audit (`alias_agrees_with_manifest`), never for claim content, retention, priors, pathway conditions or engine arguments. Metabolic strata, seed reaction ids and universal-model flags are not read into the claims. A test exports with entirely different truth columns and requires byte-identical claim files.

Limitation that cannot be removed here: the artifact was extracted for the Price targets, so the claim files contain the target proteins and their neighbourhoods only. Every target is an enzyme with a known EC; accuracy is measured among these known enzymes, not as enzyme detection.

## Truncation, coordinates and full context

The compact artifact is **not the whole proteome** and holds **no gene coordinates**. None are invented from gene order, so no `genome_layout` is exported and layout-dependent arms are not in the manifests. `export_report.json` states this under `truncation` and lists per predictor and genome the target-candidate source, nonzero counts, candidates lost by the artifact's own truncation, score mass kept, claim and context proteins, duplicate merges, source-score mismatches, non-complete EC claims and tied top scores. `candidate_coverage` summarises per predictor; `full_candidate_coverage` is true only if every genome exported all nonzero candidates, and the runner manifests say so in their description (otherwise "PRUNED or truncated").

The original full-matrix schema is supported through `--raw-matrix-map` (`.npz` with `index`, `columns`, `values`, or the original pandas pickle with `--trust-pickle`; unpickling executes code). Rows are proteins, columns ECs with optional `EC:` prefix; shape, NaN/inf/negative/>1 values are validated. **Retention policy is explicit**: `all_nonzero` (default; nothing pruned) or `threshold_topk` (per protein every score ≥ `--min-score` plus the `--top-k` highest; lossy). A pruned export reports `pruned: true`, `full_candidate_coverage: false` and "PRUNED" in the manifest description. The whole matrix is processed one protein row at a time but held in memory as loaded.

**Raw file retrieval manifest:** `raw_files_needed.tsv` in the export directory (and a short `raw_retrieval_manifest` summary in the report) lists the exact ibex path, sha256, size and mtime of each of the 132 raw prediction files (6 predictors × 22 assemblies, 23.6 GB), read from the artifact metadata, under `/ibex/scratch/projects/c2014/kexin/funcarve/…` (owner `kexin`, read-only). Gene coordinates must come from the genome annotation (GFF), which is not among them. METEOR files are not needed.

## Arms, matched evidence and raw-rank preservation

One manifest per predictor and per ensemble (`manifests/v3_price_<variant>.json`); inside a manifest every arm reads the identical claims file, so evidence is matched across arms by construction: `stage0` (`--enable-priors none`) and `stage1` (default priors with the KEGG reference files). `ensemble_pool4` is exactly the union of the four single-predictor claim sets (CLEAN, DeepProZyme, MAPred, EnzBERT; each (protein, EC, source) once, enforced). Single predictors are compared through separate manifests. Arms that need a genome layout or whole-proteome context are added only after full-context export.

**Do not equate `stage0` with the unmodified raw ranking without checking.** The engine's posterior is a sigmoid of the score, and the batch writer in this worktree serialises posteriors with `%.6f`, which collapses very small CLEAN differences into ties (and a tie can flip top-1). Measured on GCF_001633105.1 with the full CLEAN claims: 15,227 of 15,244 adjacent raw-distinct candidate pairs became posterior ties. The root's `Double.toString` serialisation fix should remove this; `check-roundtrip` verifies it for any engine build: it compares the raw order of a single-source claims file with the posterior order and reports ties created, inversions and changed top-1, and exits non-zero unless the raw order is fully preserved. A raw-score baseline may be claimed only after it passes. The exporter itself writes raw scores bit-exactly (a test with distinct float32 values one ULP apart checks the JSON round trip). Thresholds are not tuned on Price: the fixed 0.5 is a grid constant, and with the provisional calibrator top-1 accuracy is the meaningful metric.

## Status of the generated export (`/tmp/gspa-v3-price-full-inputs`, not committed)

148 labelled targets in 22 genomes, 151 label pairs, 1 exclusion; all six predictors with `full_candidate_coverage: true` and no fallback. Claims: CLEAN 775,830, DeepProZyme 406,718, EnzBERT 16,260, MAPred 8,322, GraphEC 755,266 (156 MB), TopEC 2,352, `ensemble_pool4` 1,207,130 (246 MB). No benchmark has been run on these inputs.
