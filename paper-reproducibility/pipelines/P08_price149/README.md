# P08 — Price-149 RB-TnSeq EC benchmark (six enzyme predictors)

Manuscript elements: Table 6, Fig. 6b (Price-149 / 22 host genomes).
Compute: ibex (extraction of the raw predictor outputs, Slurm `debug`, account
`pi-hohndor`) + local CPU (evaluation, ~8 min).

## Contents

| File | Role |
| --- | --- |
| `run.sh` | `submit` (ibex: 6-task extraction array + dependent merge job), `extract <pred>`, `merge`, `fetch` (local: rsync pickle + sha256 TSV), `eval` (local) |
| `extract_price149_all6.py` | Reproducible extractor (v2). Writes `price149_all6_compact.pkl.gz` with full target score vectors and per-predictor raw-format statistics in `_meta` |
| `input_sha256.tsv` | sha256, size and mtime of the 6 × 22 raw prediction pickles and of the manifest, computed on ibex inside extraction array job 53987408_[0-5] (identical in the earlier runs 53986797 and 53987356) |
| `manifest_rows.tsv` | The same inputs in the `paper-reproducibility/manifest/data_manifest.tsv` schema, plus the compact pickle (to be merged into the central manifest) |
| `OUTPUTS.md` | Result files written and their columns |

Evaluation entry point (in the repo, not copied): `scripts/eval_price149_real_all6.py`.

## Reproduce

```bash
# on ibex (login node only submits; all loading of the pickles happens in sbatch jobs)
mkdir -p /ibex/scratch/projects/c2014/hohndor/price149_repro_2026-10-19/P08_price149
scp run.sh extract_price149_all6.py ibex:/ibex/scratch/projects/c2014/hohndor/price149_repro_2026-10-19/P08_price149/
ssh ibex 'bash /ibex/scratch/projects/c2014/hohndor/price149_repro_2026-10-19/P08_price149/run.sh submit'
# locally, after the merge job has finished
bash run.sh fetch     # keeps the previous pickle as price149_all6_compact.pkl.gz.old
bash run.sh eval      # writes results/pillar2_price149_{all6_real,by_stratum,per_protein_real}.tsv
```

The extraction takes about 1.5 min per predictor (6 parallel array tasks, ≤7.3 GB RSS),
the merge 30 s. The local evaluation solves one minimal-glucose LP per genome and
predictor (132 solves) and takes ~8 min on one core.

## Remote inputs (read-only; colleague's tree)

Root `F=/ibex/scratch/projects/c2014/kexin/funcarve` (owned by the c2014 project; do not write).

| Predictor | Raw files (22 host assemblies) | File dates (mtime) |
| --- | --- | --- |
| CLEAN (2023) | `$F/paperA_2026/baseline_preds/price22_CLEAN_clean/<GCA>/<GCA>_CLEAN_confidence.pkl` | 2026-06-26/27 |
| DeepProZyme (2023) | `$F/dpec2_result/result_price/<GCx_number>_DeepECv2_t5.pkl` (76 files in the folder; matched by accession number) | 2025-10-26 – 2026-07-04 |
| EnzBERT (2023) | `$F/tfpc/resultprice_newg/<GCA>_enzbert_predictions.pkl` | 2026-01-11 |
| GraphEC (2024) | `$F/graphec_price_new/<GCA>_GraphEC.pkl` | 2026-06-03 |
| MAPred (2025) | `$F/mapred_price_new/<GCA>_MAPred.pkl` | 2026-05-18 |
| TopEC (2025) | `$F/topec_price_new/<GCA>_TopEC.pkl` | 2026-06-04 |

Also read: `$F/meteor_v8/data/all_ancestors.txt` (7,433-EC vocabulary),
`$F/meteor_v8/data/enzymeobsolete.pkl` (transferred-EC map),
`$F/meteor_v8_evw_p2mu3_run/meteor_out_price/<pred>/meteor_preds_<GCA>.pkl` (METEOR v8
active/muted EC sets, written 2026-07-26), and the manifest
`/ibex/scratch/projects/c2014/hohndor/suite6_price149_manifest.tsv` (149 rows; sha256
`dbbbb6f9…c2dcc6`). It differs from the local `benchmark/manifests/suite6_price149_manifest.tsv`
(sha256 `9c4ecf68…0881c1`) only in `cognate_exchange_rxn` / `cognate_substrate_name`: the
ibex copy keeps the assay substrates, the local copy is blinded to glucose minimal medium
(`EX_cpd00027_e` for all rows). The extractor reads only `price_entry`, `host_assembly`
and `proteome_protein_id` (identical in both); the evaluation reads the local copy.
The tool versions themselves (model checkpoints) are not recorded next to the pickles;
the dates above are the only provenance available on ibex.

Writes go only to `/ibex/scratch/projects/c2014/hohndor/price149_repro_2026-10-19/`
(`P08_price149/`, `logs/`, `out/parts/`, `out/price149_all6_compact.pkl.gz`,
`out/input_sha256.tsv`).

## Raw formats (inspection job 53986007)

All six raw files are pandas DataFrames, rows = proteins in proteome order, columns =
4-level ECs (no partial ECs; GraphEC/MAPred/TopEC share one 7,062/7,068-label set with an
`EC:` prefix). No column falls outside the 7,433-EC vocabulary for any genome.

| Predictor | Shape (GCA_000006965.1) | dtype | Non-zero ECs per protein | Scores |
| --- | --- | --- | --- | --- |
| CLEAN | 6,233 × 5,242 | float64 | all 5,242 (dense) | 1e-58 – 0.999; median 0 ECs ≥ 0.01, max 8 |
| DeepProZyme (`_t5`) | 6,233 × 2,802 | float32 | all 2,802 (dense) | 2e-25 – 1; median 72 ECs ≥ 0.01 |
| EnzBERT | 6,233 × 4,742 | float64 | exactly 100 (top-100 softmax) | row sums 0.16–1 |
| GraphEC | 6,233 × 7,062 | float32 | median 5,071 (dense sigmoid) | max 0.999; median 228 ECs ≥ 0.01; hub labels such as 4.2.99.18 and 3.6.1.31 score highly on many proteins |
| MAPred | 6,233 × 7,062 | float32 | median 13, max 183 (scores < 0.01 zeroed) | 0.01 – 1 |
| TopEC | 6,233 × 7,068 | float64 | **exactly 1** | the argmax class with its probability (0.06–1) |

Findings on the suspected defects:

* **TopEC is one EC per protein in the raw files themselves**, not an extraction artefact.
  Every copy on ibex (`topec_price/` from 2026-04-05, `topec_price_new/` and the identical
  tax-named copies in `paperA_2026/baseline_preds/price22_TopEC_clean/` from 2026-06-04)
  stores only the argmax EC and its softmax probability; the full distribution was not
  kept. TopEC therefore stays in the benchmark as a single-call predictor (Top-1 = Top-5).
* **GraphEC is a dense multi-label output** and the v1 extraction read it correctly; its
  low accuracy (18/148 Top-1) is a property of the raw scores (hub labels outrank the
  true EC). v1 truncated each target to the 100 highest-scoring ECs; v2 keeps all ECs
  ≥ 0.01 (median 2,141 per target), which only affects threshold-based metrics.
* **GraphEC vs MAPred for GCA_000006965.1**: identical byte size (176,253,393) because both
  are float32 DataFrames of identical shape (6,233 × 7,062, same index and columns), but the
  content differs (md5 `b75a9a7e…` vs `c216a9af…`; `DataFrame.equals` False; MAPred has
  13 non-zero ECs per protein, GraphEC 5,071). `graphec_price_new/GCA_000006965.1` is
  byte-identical to `graphec_price/GCA_000006965`; `mapred_price_new` is a float32 copy of
  `mapred_price` (same values).
* Top-1 EC of every one of the 148 targets is identical between v1 and v2 for all six
  predictors.
* Manifest row `NP_384885` (GCA_000006965.1) has `proteome_protein_id = NC_003047.1`, a
  replicon accession, so it matches no protein in any predictor; 148 of the 149 targets are
  evaluated (as in v1).

## Extraction schema (v2, `extraction_version = v2-2026-10-19`)

Per predictor and host assembly: `target_full` (float32, 148-row matrix over the 7,433-EC
vocabulary: the lossless full vector of each target), `target_scores` (all ECs with
score ≥ 0.01 plus the 100 highest-scoring ECs, unrounded; lossless for Top-k, MRR within
rank 100 and threshold F-max at τ ≥ 0.01), `ec_top_scores` (8 best proteins per EC with
score ≥ 0.01, genome-wide), `neighbor_ecs` (top-5 ECs ≥ 0.15 of the ±10 neighbouring
proteins), METEOR active/muted ECs, `raw_target_stats`. `_meta` holds the vocabulary,
manifest sha256, per-file sha256/size/mtime and per-genome raw statistics (orientation,
column counts, score range, non-zeros per protein). Obsolete ECs are mapped to their
current EC (max-merged) before reindexing.

v1 differences: v1 capped `target_scores` at 100 ECs and rounded to 5 decimals (CLEAN's
small scores became 0.0 ties); otherwise identical genome-level fields.

## Evaluation arms (`scripts/eval_price149_real_all6.py`)

* **Baseline**: the predictor's own target vector.
* **METEOR v8 (Iter 1 / Iter 2)**: nudge from the archived METEOR v8 sets for that predictor
  and genome: `s ← s + 0.4·s(1−s)` if the EC is in `meteor_active_ecs`, `s ← s − 0.4·s(1−s)`
  if in `meteor_muted_ecs`, unchanged otherwise; Iter 2 applies the nudge twice. METEOR
  itself is not re-run.
* **GSPA SLR (Minimal Glucose)**, **SLR + RLGC**, **Single-Pred (Minimal Glucose)**,
  **Full GSPA (Iter 1–3)**: single-predictor arms (only this predictor's genome matrix,
  ±10 neighbours, the universal S-matrix and one minimal-glucose LP per genome).
* **Stage 1 Noisy-OR Ensemble (Sequence Only)**, **Full GSPA Multi-Pred (Minimal Glucose /
  Iter 3)**: `p0 = 0.42·s_baseline + 0.58·noisyOR(CLEAN 0.38, DeepProZyme 0.34, MAPred 0.34,
  EnzBERT 0.28)` — the composition is written to the `pooled_predictors` column. For CLEAN,
  DeepProZyme, MAPred and EnzBERT rows the baseline is also inside the pool.
* **pool_without_this_baseline**: the same Stage 0 formula with the row's baseline removed
  from the pool and its direct score set to 0 (`0.75·noisyOR(others)`), so the difference to
  the Stage 1 row is the baseline's own contribution. For GraphEC and TopEC (not pooled)
  it is the four-predictor pool alone.
* McNemar: exact two-sided binomial test (`scipy.stats.binomtest`) on the discordant Top-1
  pairs against the row's Baseline (`n_improved_vs_base` = b, `n_degraded_vs_base` = c).

## v1 → v2 result changes (`results/pillar2_price149_all6_real.tsv`)

Top-1 is unchanged for every predictor's Baseline. Changes come only from v1's 5-decimal
rounding and top-100 cap:

* CLEAN Baseline Top-5 90 → 95 and per-protein F-max 0.529 → 0.555 (v1 rounded CLEAN's
  small scores to 0.0, creating ties that were broken alphabetically); CLEAN Multi-Pred
  Iter 3 Top-1 68 → 70 (b/c 14/6 → 14/4, p 0.115 → 0.031).
* GraphEC Multi-Pred arms lose 2–4 Top-1 hits (Stage 0 59 → 56, Multi-Pred Iter 3
  67 → 65) because GraphEC's full ≥ 0.01 candidate list now enters the pool.
* DeepProZyme, MAPred, TopEC and EnzBERT rows are unchanged up to ±0.002 in F-max
  (EnzBERT METEOR Iter 1 Top-5 80 → 78).
* `mcnemar_p` is now printed with 4 significant digits (v1 rounded to 6 decimals, so the
  GraphEC/TopEC Multi-Pred p-values showed as 0).

`pool_without_this_baseline` shows that for GraphEC and TopEC the four-predictor pool on
its own (68/148 Top-1, macro F-max 0.474) is better than every Multi-Pred arm that adds
GraphEC (≤ 65) or TopEC (≤ 59); their Multi-Pred gains over Baseline come from the pool.

The evaluation script's per-protein and dataset-level F-max loops were rewritten as
sorted sweeps (same thresholds, same results: re-running on the v1 pickle reproduces the
committed v1 table exactly except for the p-value print precision); the quadratic version
did not finish within an hour on the v2 GraphEC vectors.
