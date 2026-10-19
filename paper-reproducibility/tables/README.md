# Generated manuscript tables

`make tables` (from `paper-reproducibility/`) runs

    python3 tables/make_tables.py --repo ..

and writes one file per manuscript table:

| Label | Output | Results files read |
| --- | --- | --- |
| `tab:cafa_splits` | `paper/tables/cafa_splits.tex` | `results/cafa_whole_genome_time_and_homology_splits.tsv` |
| `tab:ablation` | `paper/tables/ablation.tex` | `results/pillar1_2_ablation_summary.tsv`, `results/pillar1_2_perturbation_summary.tsv`, `results/pillar1_2_ablation_per_genome.tsv` (consistency column) |
| `tab:complex_lso` | `paper/tables/complex_lso.tex` | `results/pillar3_complex_lso.tsv`: rows by `decoy_design` x `genome` (`ALL (Weighted Mean)`, `Out-of-Sample 5 (Weighted Mean)`, `synechocystis`) x `method_id` (input_attenuated, real_diamond_channel, real_interproscan_channel, operon_only, coherence, coherence_operon, full_stage3; `noisy_or` and `coherence_operon_stage2` omitted); bold = column max within each design x cohort sub-block |
| `tab:supp_lso_mcnemar` (supplement) | `paper/tables/supp_lso_mcnemar.tex` | `results/pillar3_complex_lso_mcnemar.tsv` (`n_cases`, `b_a_only`, `c_b_only`, `p_exact_two_sided`) |
| `tab:anatomy` | `paper/tables/anatomy.tex` | `results/pillar3_eukaryote_anatomy.tsv`, `results/pillar3_relax_vs_remove.tsv` (`scope` = `injected_cohort_all_proteins` / `injected_cohort_injected_proteins`), `results/pillar3_eukaryote_endogenous_vs_injected.tsv`; only rows with `input_source` = `real` or `jar:*` (never `legacy:*` / `synthetic`). Unperturbed row: endogenous E1-E5 violations and measured CAFA protein-centric F-max (`endog_e1_fmax`, `endog_e4_remove_fmax`, `endog_e5_relax_fmax`) |
| `tab:metabolism` | `paper/tables/metabolism.tex` | `results/pillar2_lro15_summary.tsv`, `results/pillar4_lro_50_cases.jsonl` (single-feature rows via `benchmark/unified/pillar4/evaluate_pillar4.py:evaluate_lro_single_feature_baselines` / `rank_lro_by_feature`), `results/pillar4_v2_multiseed_summary.tsv` (Stage 3 row: benchmark v2, de-identified prompts, 3-seed mean) |
| `tab:price149` | `paper/tables/price149.tex` | `results/pillar2_price149_all6_real.tsv` (`top1_count`, `fmax_dataset_macro`; newest `extraction_version` only; "Pool without this predictor" = method `pool_without_this_baseline`) |
| `tab:bacdive` | `paper/tables/bacdive.tex` | `results/pillar3_bacdive_track1_summary.tsv`, `results/pillar3_bacdive_track2_summary.tsv` |

Supplement: `tab:supp_perturbation_designs` -> `paper/tables/supp_perturbation_designs.tex` from
`results/pillar1_2_perturbation_designs.tsv` (one row per design x jar; per config A1/A6/A7/A8/A9 the
survivors summed over the nine genomes, out of the summed `injected_pairs`, and the mean
`fmax_micro_exp`). `supplement.tex` is expected to `\input{tables/supp_perturbation_designs.tex}`.

Supplement, further generated tables (each `paper/tables/<label>.tex`, `\input` by `supplement.tex`):

| Label | Suppl. table | Results files read |
| --- | --- | --- |
| `tab:supp_panels` | S2 | `results/pillar1_suite2_panel13_21_sota.tsv` (rows with `source=computed` carry numbers; other rows are printed as "not evaluated") |
| `tab:supp_cafa_wg_splits` | S4 (longtable, 130 rows) | `results/cafa_whole_genome_time_and_homology_splits.tsv` (incl. `violations_050`, `hom_scale_selected`) |
| `tab:supp_1000bacteria` | S8 | `results/pillar1_suite2_1000bacteria_sota.tsv` (`source=computed` rows; other rows printed as pending) |
| `tab:supp_mags` | S11 | `results/pillar3_dark_matter_mags_v2.tsv` (per-phylum rows and `ALL`) |
| `tab:supp_multikingdom_domains` | S9 | `results/multikingdom_summary_by_domain.tsv`, `results/multikingdom_bootstrap_ci.tsv` (written by `pipelines/P15_multikingdom/bootstrap_ci.py`), `results/pillar4_v2_trigger_rates.tsv` (`SUMMARY_*` pooled capped rate) |
| `tab:supp_invariant_audit` | S12 | `results/pillar2_invariant_audit_summary.tsv`, `results/pillar2_curated6_gems_real.tsv` (`MEAN_6GEMS`), `results/pillar2_price149_by_stratum.tsv` (newest `extraction_version`, CLEAN) |
| `tab:supp_bacdive_folds` | S13 | `results/pillar3_bacdive_per_fold_cv.tsv`, `results/pillar3_bacdive_track1_summary.tsv`; strain counts as in `tab:bacdive` |
| `tab:supp_atb_survey` | S16 | `results/allthebacteria/insight3b_phylum_summary.tsv`, `insight5b_phylum_nogd_vs_divergent.tsv` (alignment columns only; phyla with < 15 genomes summed) |
| `tab:supp_atb_integrator` | S17 | `results/allthebacteria/atb_java_summary_stageD_snapshot.json`, `atb_java_summary.json`, `insight3_per_genome_java.tsv` (median runtime) |
| `tab:supp_stage3_v2` | S14 | `results/pillar4_v2_multiseed_summary.tsv`, `results/pillar4_v2_per_category.tsv` |

Supplementary table numbers follow LaTeX order in `supplement.tex`; the numbers above are those of the
2026-10-19 revision.

Each output holds the complete `\begin{table*} ... \end{table*}` block.

## What is generated and what is static

* **Body (every numeric cell)**: read from the results files above, selected by
  column values (split id, method, config, genome, cohort, ...), never by row
  position. The mapping (table, block, row, column) -> (file, row selector,
  column, derivation) for every cell is in `table_sources.tsv`.
* **Derived cells** (documented in the `derivation` column of
  `table_sources.tsv`): ablation consistency (mean over genomes of
  `max(0, 1 - 0.01 V_g)`), perturbation survivors and suppression (means over
  the nine genomes), out-of-sample 5-species LSO rows (n-weighted means over
  the non-calibration genomes), anatomy totals / percentages, LRO Hit@1
  counts (`round(hit1_all * n_total)`), Price-149 gains, BacDive McNemar p
  values (recomputed from b, c by an exact binomial test where the TSV stores
  an underflowed 0), and BacDive strain counts (parsed from the dataset id).
* **Static text**: caption, column header and footer live in
  `templates/<label>.head.tex` / `templates/<label>.foot.tex`. They were
  extracted once from `paper/main.tex` (`make_tables.py --extract-templates`).
  After the switch described below, edit captions there, not in main.tex.
  Numbers inside captions (e.g. `$20{,}697$ ORFs`) are static and are not
  checked by either tool.
* **Unsourced cells** (no results file holds them) are omitted and listed on
  stderr and in `UNSOURCED.md`, which is regenerated on every run together
  with the cell-by-cell diff against the block currently in `main.tex`.

## Bold

* `tab:cafa_splits`, `tab:ablation`, `tab:complex_lso`, `tab:bacdive`:
  computed column-wise maximum within a block (minimum for violation and
  surviving-claim counts), on the rounded values, ties all bold. Comparison
  sets: cafa_splits excludes the InterPro-LR row (as the current manuscript
  does); bacdive compares {gapseq, GSPA Track 1} and {RF, GBDT, L1-SVM,
  GSPA + GBDT} separately in the Track 1 blocks, all rows in the other blocks,
  and bolds the McNemar cell with the largest b - c.
* `tab:price149`: row-wise, the arm with the highest Top-1 count among METEOR,
  SLR + RLGC, single-predictor, Stage 0 pooling, pool without this predictor
  and multi-predictor (Standalone excluded); ties all bold; F-max is not
  bolded separately.
* `tab:metabolism` (Full-GSPA / Stage 3 rows) keeps the current row
  highlighting; a `BOLD-SEMANTICS` warning is printed where a highlighted cell
  is not the block maximum.
* `tab:anatomy` has no bold data cells.

## Schema changes in results files

If a results file gains a column that splits one logical row into several
(e.g. `design`, `decoy_design`, `input_source`), the script warns
`AMBIGUOUS ...` and uses the first match. Choose the intended rows with

    python3 tables/make_tables.py --repo .. --filter design=<value> [--filter ...]

(the filter applies to every file that has that column).

## Adopting the generated tables in main.tex (not done yet)

Once the diff in `UNSOURCED.md` has been reviewed, replace each table block in
`paper/main.tex`, from its `\begin{table*}[t!]` line to its `\end{table*}`
line, by a single line, e.g.

    \input{tables/cafa_splits.tex}

(`main.tex` lives in `paper/`, so the path is relative to `paper/`; the same
for `ablation`, `complex_lso`, `anatomy`, `metabolism`, `price149`,
`bacdive`). Leave the `\label{tab:...}` inside the generated file; references
keep working. After the switch, `make tables` must be run before building the
paper; validation then compares against the block in `git HEAD:paper/main.tex` (or
`supplement.tex`), which is only meaningful while HEAD still holds the hand-written table.
Before switching `tab:bacdive`, drop the urease sentence from
`templates/bacdive.head.tex` or source the urease rows.
