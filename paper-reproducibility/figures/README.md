# Manuscript figures

`generate_figures.py` regenerates Figures 1-6 of `paper/main.tex` from `results/*.tsv`
and writes `paper/figures/<stem>.pdf` and `.png` (300 dpi). Run it from the repository
root or via `make figures` in `paper-reproducibility/`:

    python3 paper-reproducibility/figures/generate_figures.py        # Figures 1-6
    python3 paper-reproducibility/figures/generate_figures.py 3 6    # selected figures

Every plotted value is selected with an explicit `column=value` row selector (the same
selectors as `../tables/table_sources.tsv`); a selector that matches no row or more than
one row stops the run. Each run also writes `figure_sources.tsv`, one line per plotted
value: figure, panel, source file, row selector, source column, derivation, value (a run
of selected figures replaces only those figures' rows). No
number is typed into the scripts; the only constants are benchmark design parameters
(chance level 1/6 in Fig. 4a = one true subunit among six candidates) and the display
threshold of four MAGs per phylum in Fig. 3d.

`paper/generate_figures.py`, `paper/generate_figure6.py` and `paper/generate_figure7.py`
are thin wrappers that call this folder.

| File | Content |
| --- | --- |
| `common.py` | palette, typography, sub-figure grid, row selection with provenance |
| `fig1_architecture.py` | Fig. 1 (text-box diagram, text corrections only) |
| `fig2_ablation.py`, `fig3_metabolism.py`, `fig4_complex_anatomy.py`, `fig5_agent.py`, `fig6_fba_bacdive.py` | Figs. 2-6 |
| `generate_figure7.py` | Fig. 7 (AllTheBacteria), **not run**; awaits the P11 rerun. Only the repository path was changed; the current `fig7_*` files are from the earlier results |
| `figure_sources.tsv` | provenance of every plotted value (written by each full run) |

## Visual grammar

All figures use one palette (the dataviz reference palette, light mode, validated with
its `validate_palette.js`): the GSPA refinement stages form a one-hue blue ramp from light
(first prior) to dark (full GSPA); unconstrained pooling, no-prior and single-feature
references are neutral grey; external tools and rule baselines take the non-blue
categorical slots (orange, aqua, yellow); paired conditions (decoy design, capped versus
uncapped) are blue and orange. Typeface Liberation Sans (Arial metrics), 12 pt on an
11.8 x 8.3 in canvas, which prints near 7 pt at `0.96\textwidth`. Panel letters a-d sit
with a sentence-case title at the top left of each panel.

## Panels and source rows

### Figure 1 - `fig1_architecture.{pdf,png}`
Text-box diagram kept from the earlier version (redesign is a separate task). Text fixes:
`Qwen3.8-27B` and the Stage 3 trigger line "Inspects about N% of loci under the C-anomaly
cap", with N = round of `pct` `pooled=` in `pillar4_v2_trigger_rates.tsv` `genome=SUMMARY_ALL`
and C = max `n_anomalies` over per-genome rows with `cap_reached=True`.

### Figure 2 - `fig2_pillar1_2_ablation`
| Panel | Source | Rows | Columns |
| --- | --- | --- | --- |
| a | `pillar1_2_ablation_summary.tsv` | `config=A1_noisy_or ... A9_stage2_relax` (nine-genome means; `*_non_ecoli8` rows not used) | `composite_score`, `completeness`, `complex_coherence`, `pathway_coherence` |
| b | `pillar1_2_ablation_per_genome.tsv` | each `genome`; `config=A1_noisy_or` and `config=A9_stage2_relax` | `composite_score` |
| c | `pillar1_2_perturbation_designs.tsv` | `design=v1` and `design=v2`; `jar=gspa-unified.jar; claims_ids=normalised`; `config` in A1, A6, A7, A8, A9; each genome | `surviving_violations` (total row = sum over genomes; denominator = sum of `injected_pairs`) |
| d | `pillar1_suite2_panel13_21_sota.tsv` | `panel=Panel13_Prok13`; the six methods of `tab:ablation` | experimental truth: `go_fmax_exp_prop` (micro), `go_fmax_cafa_exp_prop` (CAFA); all-evidence truth: `go_fmax_micro_all`, `go_fmax_cafa_all` |

### Figure 3 - `fig3_pillar3_metabolism`
| Panel | Source | Rows | Columns |
| --- | --- | --- | --- |
| a | `pillar2_lro15_summary.tsv` | `benchmark_cohort=Ecoli_MG1655_Strict_256; stratum=All_Strata; protocol_H=Unfiltered`; methods `Sequence_Only_Candidate_Rank`, `GSPA_Strict_Homology_Masked_TrackA`, `RLGC_Kernel_Density_Only`, `GSPA_Stage1_M2_Coherence_RLGC`, `ESM2_PLM_Centroids_Only`, `GSPA_Full_M4_RLGC_Coherence_PLM` (prefix match) | `hit1_all`, `hit1_cov`, `hit3_cov`, `n_cov`, `n_total` |
| b | `pillar2_lro15_summary.tsv` | `benchmark_cohort=Ecoli_MG1655_Full_414; protocol_H=Unfiltered; stratum` in A/B/C; methods `gapseq_Unassigned_Operon_Density`, `RLGC_Kernel_Density_Only`, `GSPA_Full_M4_RLGC_Coherence_PLM` | `hit1_all`, `n_total` |
| c | `pillar2_lro15_summary.tsv`; `pillar4_v2_ml_baselines.tsv`; `pillar4_v2_multiseed_summary.tsv` | `benchmark_cohort=Bsubtilis_168_LRO_144` Track A / Track B; `method=Single Feature: ...` (three rows); `seed_label=MEAN_SD` | `hit1_all`, `hit1_cov`; `lro50_hit1_pct`/100; `lro_hit1_pct` (mean +/- s.d.)/100 |
| d | `pillar3_dark_matter_mags_v2.tsv`; `pillar3_dark_matter_mags_v2_per_mag.tsv` (P14) | `phylum_or_group` with `n_mags` >= 4 (all 17 phyla at present) plus `phylum_or_group=ALL`; per-MAG points: rows with `phylum=<group>` (all rows for ALL) | bars `mean_fills_per_mag_q050`, points `new_fills_q050` (fixed jitter, seed 20261019); paired markers `pathway_completeness_before` / `pathway_completeness_after`; labels `n_mags` |

### Figure 4 - `fig4_pillar3_complex_anatomy`
| Panel | Source | Rows | Columns |
| --- | --- | --- | --- |
| a | `pillar3_complex_lso.tsv` | `genome=ALL (Weighted Mean)`; `decoy_design` in `unrelated_loci`, `within_operon`; the seven `method_id` rows of `tab:complex_lso` (`noisy_or` and `coherence_operon_stage2` omitted, as in the table) | `hit_at_1`, `n_cases` |
| b | `pillar3_complex_lso.tsv` | `decoy_design=unrelated_loci`; six genomes; `method_id` in `input_attenuated`, `coherence`, `coherence_operon`, `full_stage3` | `f1_at_05`, `n_cases` |
| c | `pillar3_eukaryote_anatomy.tsv` | `input_source=real`; six species | `n_violations_e1_noisy_or`, `_e2_go_taxon`, `_e3_stage1_anat`, `_e4_stage2_remove`, `_e5_stage2_relax` (E4 and E5 drawn as one bar when identical), `n_cross_clade_injected` |
| d | `pillar3_relax_vs_remove.tsv` | `input_source=real; scope=injected_cohort_all_proteins`; `mode` in `E1_noisy_or`, `E4_stage2_remove`, `E5_stage2_relax` | `cafa_fmax` (the pooled micro F-max of `tab:anatomy`); mean row = mean over six species |

### Figure 5 - `fig5_pillar4_agent` (Stage 3, benchmark v2)
| Panel | Source | Rows | Columns |
| --- | --- | --- | --- |
| a | `pillar4_v2_evaluation.tsv`; `pillar4_v2_multiseed_summary.tsv` | the configurations of `tab:agent` (best supervised per feature set = highest `f1`, then `accuracy`, among `... + Guardrail` rows); guarded agent = `seed_label=MEAN_SD` | `precision`, `recall`, `f1`, `accuracy` |
| b | `pillar4_v2_evaluation.tsv` | same rows; guarded agent = max over the three seed rows | `violations_accepted` |
| c | `pillar4_v2_pg_tp_hit1.tsv`; `pillar4_v2_ml_baselines.tsv`; `pillar4_v2_multiseed_summary.tsv` | heuristics, zero-shot, fast, unguarded, three guarded seeds (mean, population s.d.); LRO-50 single features and agent rows | `hit1_pct`/100; `lro50_hit1_pct`/100; `lro_hit1_pct`/100 |
| d | `pillar4_v2_trigger_rates.tsv` | `genome=SUMMARY_bacteria/archaea/eukaryote/virus/ALL` | `pct` and `pct_uncapped` (`pooled=` value), `taxon` (`n_genomes=`) |

### Figure 6 - `fig6_pillar2_3_fba_bacdive`
| Panel | Source | Rows | Columns |
| --- | --- | --- | --- |
| a | `pillar2_invariant_audit_summary.tsv` | `method` = METEOR v8, Grow-Support Threshold Repair, GSPA StoichiometricLocusRefiner (the raw-threshold row is not plotted) | `mean_added_rxns`, `mean_flux_certified_added`, `mean_uncertified_padding`, `mean_imbalanced_added`; feasibility strip `lp_biomass_feasible_pct_own_medium` (empty for GSPA, whose own medium is the audit medium) and `lp_biomass_feasible_pct` |
| b | `pillar2_price149_all6_real.tsv` | `extraction_version` = newest (`v2-2026-10-19` at the time of writing); six `baseline`s; `method` in `Baseline`, `pool_without_this_baseline`, `+ Stage 1 Noisy-OR Ensemble (Sequence Only)`, `+ Full GSPA Multi-Pred (Iter 3)` | `top1_count`, `fmax_dataset_macro`, `n_proteins` |
| c | `pillar2_lro15_summary.tsv` | `benchmark_cohort=Ecoli_MG1655_Full_414; stratum=All_Strata; method=GSPA_Full_M4_RLGC_Coherence_PLM*`; `protocol_H` in Unfiltered, le70, le50, le30_Dark | `hit1_cov`, `cov_frac`, `n_cov`, `n_total` |
| d | `pillar3_bacdive_track1_summary.tsv`, `pillar3_bacdive_track2_summary.tsv` | the Track 1 and Track 2 rows of `tab:bacdive` (same `dataset`/`method` selectors); 10-EC diagnostic panel rows not plotted; no urease rows exist in these files | `f1`, `roc_auc` |
