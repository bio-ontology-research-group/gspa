# Provenance map: manuscript items to results files, pipelines and generators

One row per table, figure or text section of `paper/main.tex` and `paper/supplement.tex`
(numbering of the 2026-10-19 revision; supplementary numbers follow LaTeX order).

* **Pipeline**: folder under `paper-reproducibility/pipelines/` that (re)computes the results
  file. P01, P06 and P09 are listed in `paper-reproducibility/README.md` but have no folder yet;
  their scripts are named instead.
* **Generator**: the script that turns the results file into the manuscript item.
  `make_tables.py` = `paper-reproducibility/tables/make_tables.py` (writes `paper/tables/<label>.tex`);
  figure scripts are under `paper-reproducibility/figures/`.
* **Status**: `computed` = the results file is written by a committed script from raw outputs
  listed in the pipeline's `OUTPUTS.md` / `manifest_rows.tsv`; `archived` = the results file is
  produced by a committed script from archived inputs that are not regenerated here, or is a
  transcription; `pending` = regeneration not finished, numbers in the manuscript are provisional;
  `hand` = the manuscript item is typed in the `.tex` source (descriptive tables only).

## Main text

| Item | Results file(s) | Pipeline | Generator | Status |
| --- | --- | --- | --- | --- |
| Table 1 (`tab:cafa_splits`) | `cafa_whole_genome_time_and_homology_splits.tsv`, `cafa_whole_genome_per_species.tsv` | P01 (`scripts/run_cafa_whole_genome_splits.py`; split-retrained caches built on orix by `benchmark/whole_genome_cafa/scripts/run_orix_wg_cafa.sbatch` and `scripts/train_eval_external_sota_wg_cafa.py`) | `make_tables.py` (`cafa_splits`) | computed (being re-run after the validation-alignment cache fix) |
| Table 2 (`tab:ablation`) | `pillar1_2_ablation_summary.tsv`, `pillar1_2_ablation_per_genome.tsv`, `pillar1_2_perturbation_summary.tsv` | P02 (`ablation/01_run_ablation_grid.py` … `05_compile_results.py`, `run_perturbation.py`, `collect_results.py`) | `make_tables.py` (`ablation`) | computed |
| Table 3 (`tab:complex_lso`) | `pillar3_complex_lso.tsv` | P04 (`benchmark/unified/pillar3/run_pillar3_complex_lso.py`) | `make_tables.py` (`complex_lso`) | computed |
| Table 4 (`tab:anatomy`) | `pillar3_eukaryote_anatomy.tsv`, `pillar3_relax_vs_remove.tsv`, `pillar3_eukaryote_endogenous_vs_injected.tsv` (rows `input_source=real`) | P05 (`benchmark/unified/pillar3/run_pillar3_eukaryotic_anatomy.py`) | `make_tables.py` (`anatomy`) | computed |
| Table 5 (`tab:metabolism`) | `pillar2_lro15_summary.tsv`, `pillar4_lro_50_cases.jsonl`, `pillar4_v2_multiseed_summary.tsv` (Stage 3 row) | P06 (`benchmark/leave_reaction_out/run_ablation_rlgc*.py`; summary by `benchmark/unified/pillar2_fba_lro_price149/run_pillar2_all.py`); P10 (Stage 3 row) | `make_tables.py` (`metabolism`; single-feature rows via `benchmark/unified/pillar4/evaluate_pillar4.py`) | archived (per-case LRO rankings; METEOR row asserted, see `provenance` column); Stage 3 row computed |
| Table 6 (`tab:price149`) | `pillar2_price149_all6_real.tsv` | P08 (`extract_price149_all6.py` on ibex, `scripts/eval_price149_real_all6.py`) | `make_tables.py` (`price149`) | computed (METEOR arm from archived repairs) |
| Table 7 (`tab:bacdive`) | `pillar3_bacdive_track1_summary.tsv`, `pillar3_bacdive_track2_summary.tsv` | P09 (`benchmark/unified/pillar3_bacdive/run_bacdive_tracks1_2.py`) | `make_tables.py` (`bacdive`) | computed (inputs: gapseq comparison tables) |
| Table 8 (`tab:agent`) | `pillar4_v2_evaluation.tsv`, `pillar4_v2_multiseed_summary.tsv`, `pillar4_v2_ml_baselines.tsv`, `pillar4_v2_pg_tp_hit1.tsv` | P10 (`benchmark/unified/pillar4/build_pillar4_dataset_v2.py`, `run_pillar4_v2_agent.py`, `evaluate_pillar4_v2.py`) | `make_tables.py` (`agent`) | computed |
| Figure 1 (`fig:architecture`) | `pillar4_v2_trigger_rates.tsv` | P10 | `figures/fig1_architecture.py` | computed |
| Figure 2 (`fig:pillar12`) | `pillar1_2_ablation_summary.tsv`, `pillar1_2_ablation_per_genome.tsv`, `pillar1_2_perturbation_designs.tsv`, `pillar1_suite2_panel13_21_sota.tsv` | P02, P03 | `figures/fig2_ablation.py` | computed |
| Figure 3 (`fig:complex_anatomy`) | `pillar3_complex_lso.tsv`, `pillar3_eukaryote_anatomy.tsv`, `pillar3_relax_vs_remove.tsv` | P04, P05 | `figures/fig4_complex_anatomy.py` | computed |
| Figure 4 (`fig:metabolism`) | `pillar2_lro15_summary.tsv`, `pillar4_v2_ml_baselines.tsv`, `pillar4_v2_multiseed_summary.tsv` | P06, P10 | `figures/fig3_metabolism.py` | archived (LRO) / computed (Stage 3) |
| Figure 5 (`fig:fba_bacdive`) | `pillar2_invariant_audit_summary.tsv`, `pillar2_lro15_summary.tsv`, `pillar2_price149_all6_real.tsv`, `pillar3_bacdive_track1_summary.tsv`, `pillar3_bacdive_track2_summary.tsv` | P07, P06, P08, P09 | `figures/fig6_fba_bacdive.py` | computed (LRO panel archived) |
| Figure 6 (`fig:agent`) | `pillar4_v2_evaluation.tsv`, `pillar4_v2_ml_baselines.tsv`, `pillar4_v2_multiseed_summary.tsv`, `pillar4_v2_pg_tp_hit1.tsv`, `pillar4_v2_trigger_rates.tsv` | P10 | `figures/fig5_agent.py` | computed |
| Section "Sequence dark matter across 72,533 bacterial genomes" | `allthebacteria/insight3b_phylum_summary.tsv` (alignment columns), `allthebacteria/atb_java_summary*.json` | P11 | prose | computed (alignment survey, Stage D snapshot, run (a)); full Stage D pending |
| Prose numbers | all of the above | – | `paper-reproducibility/checks/check_numbers.py` with `checks/number_ledger.tsv` | – |

## Supplement

| Item | Results file(s) | Pipeline | Generator | Status |
| --- | --- | --- | --- | --- |
| Section S1, Table S1 (`tab:supp_suites`) | (descriptive; names the files below) | all | `supplement.tex` | hand |
| Table S2 (`tab:supp_panels`) | `pillar1_suite2_panel13_21_sota.tsv` (`source=computed` rows; others shown as not evaluated) | P03 (`p03_eval.py`, `p03_panel21.py`, `p03_aggregate.py`) | `make_tables.py` (`supp_panels`) | computed |
| Section S1.1 (Stage 0 calibration) | – (code constants: `gspa-core/.../integration/CalibrationTable.groovy`, `EvidenceCombiner.groovy`, `EvidenceType.groovy`) | – | `supplement.tex` | hand |
| Section S2, Table S3 (`tab:supp_hyperparams`), S2.1–S2.3 (rescue floors, coherence weights, options) | – (code constants: `gspa-core/.../integration/prior/*.groovy`, `IntegrationState.groovy`, `IterativeRefiner.groovy`, `gspa-cli/.../IntegrateCommand.groovy`, `scripts/stoichiometric_locus_refiner.py`) | – | `supplement.tex` | hand |
| Table S4 (`tab:supp_cafa_wg_splits`) | `cafa_whole_genome_time_and_homology_splits.tsv` (130 rows, incl. `violations_050`, `hom_scale_selected`) | P01 | `make_tables.py` (`supp_cafa_wg_splits`) | computed (being re-run) |
| Table S5 (`tab:supp_holdout125`) | `scenario1_swissprot_ec_holdout125_real.tsv` | none | `scripts/eval_scenario1_holdouts_real.py` (TSV); table typed in `supplement.tex` | archived |
| Table S6 (`tab:supp_95genomes`) | `scenario1_95genomes_ec_homology_real.tsv` | none | `scripts/eval_scenario1_holdouts_real.py` (TSV); table typed in `supplement.tex` | archived |
| Table S7 (`tab:supp_bridge`) | `gspa-core/src/main/resources/taxon-constraints/*.tsv` (row counts) | – | `supplement.tex` | hand |
| Table S8 (`tab:supp_1000bacteria`) | `pillar1_suite2_1000bacteria_sota.tsv` (six tool rows `source=computed`; GSPA rows `source=unavailable` until the array finishes) | P12 (`p12_structural_suite.py`, `03_eval_chunk.sh`, `05_aggregate.py`) | `make_tables.py` (`supp_1000bacteria`) | computed (tool rows); pending (GSPA rows) |
| Table S9 (`tab:supp_multikingdom_domains`) | `multikingdom_summary_by_domain.tsv`, `multikingdom_per_genome_eval.tsv`, `multikingdom_bootstrap_ci.tsv`, `pillar4_v2_trigger_rates.tsv` | P15 (`bootstrap_ci.py`; runs by `scripts/run_multikingdom_batch_unimatrix.py`); P10 (trigger rates, `benchmark/unified/pillar4/trigger_rates_v2.py`) | `make_tables.py` (`supp_multikingdom_domains`) | computed (bootstrap, trigger rates); archived (per-proteome F-max, cluster run) |
| Table S10 (`tab:supp_bgc`) | `pillar2_3_antismash_106bgc_specificity.tsv` (to be replaced by `pillar2_3_antismash_bgc_specificity_v2.tsv`) | P13 (`jobs/00_ref` … `03_evaluate.sbatch`) | table typed in `supplement.tex` | pending |
| Section S4 and Table S11 (`tab:supp_mags`) | `pillar3_dark_matter_mags_v2.tsv`, `pillar3_dark_matter_mags_v2_per_mag.tsv` | P14 (`jobs/00_select` … `02_aggregate.sbatch`, `benchmark/unified/pillar3/mag_pathway_gaps.py`) | `make_tables.py` (`supp_mags`) | computed |
| Section S4 and Table S12 (`tab:supp_invariant_audit`) | `pillar2_invariant_audit_summary.tsv`, `pillar2_invariant_audit_108.tsv`, `pillar2_curated6_gems_real.tsv`, `pillar2_price149_by_stratum.tsv` | P07 (`run_pillar2_all.py:run_108_genome_invariant_audit`), P08 | `make_tables.py` (`supp_invariant_audit`) | computed (METEOR solutions archived; curated-GEM comparison archived) |
| Section S5 and Table S13 (`tab:supp_bacdive_folds`) | `pillar3_bacdive_per_fold_cv.tsv`, `pillar3_bacdive_track1_summary.tsv` | P09 (`benchmark/unified/pillar3_bacdive/run_bacdive_tracks1_2.py`) | `make_tables.py` (`supp_bacdive_folds`) | computed |
| Section S6 and Table S14 (`tab:supp_stage3_v2`) | `pillar4_v2_multiseed_summary.tsv`, `pillar4_v2_per_category.tsv`; text: `pillar4_v2_id_mapping.tsv`, `pillar4_v2_leak_checks.tsv`, `pillar4_v2_transcripts_*.jsonl` | P10 | `make_tables.py` (`supp_stage3_v2`) | computed |
| Table S15 (`tab:supp_cafa6`) | `suite1_cafa6_prospective_summary.tsv` (from `cafa_retrained_time_and_homology_splits.tsv`, split `Time_Split_CleanA`) | none | `scripts/eval_scenario1_holdouts_real.py` (TSV); table typed in `supplement.tex` | archived |
| Section S7, Table S16 (`tab:supp_atb_survey`) | `allthebacteria/insight3b_phylum_summary.tsv`, `allthebacteria/insight5b_phylum_nogd_vs_divergent.tsv` (alignment columns only: genomes, species, ORFs, dark, zero-alignment, twilight; derived from `rep_max_pident.tsv`) | P11 (Stages A–B) | `make_tables.py` (`supp_atb_survey`) | computed (alignment survey) |
| Section S7, Table S17 (`tab:supp_atb_integrator`) | `allthebacteria/atb_java_summary_stageD_snapshot.json`, `allthebacteria/atb_java_summary.json`, `allthebacteria/insight3_per_genome_java.tsv` | P11 (`benchmark/allthebacteria/stageD_java_integrate_batch.py`, `aggregate_atb_java.py`) | `make_tables.py` (`supp_atb_integrator`) | computed (snapshot; Stage D over all batches running) |
| Section S8, Table S18 (`tab:supp_perturbation_designs`) | `pillar1_2_perturbation_designs.tsv` | P02 | `make_tables.py` (`supp_perturbation_designs`) | computed |
| Section S9, Table S19 (`tab:supp_lso_mcnemar`) | `pillar3_complex_lso_mcnemar.tsv` | P04 | `make_tables.py` (`supp_lso_mcnemar`) | computed |
| Section S10, Table S20 (`tab:supp_agent_mcnemar`) | `pillar4_v2_mcnemar.tsv` | P10 | `make_tables.py` (`supp_agent_mcnemar`) | computed |
| Table S21 (`tab:supp_trigger_rates`) | `pillar4_v2_trigger_rates.tsv` | P10 | `make_tables.py` (`supp_trigger_rates`) | computed |

## Results files no longer used by the manuscript

`pillar4_agent_benchmark_120.jsonl`, `pillar4_agent_evaluation.tsv`, `pillar4_multiseed_summary.tsv`,
`pillar4_ml_baselines.tsv` (Stage 3 benchmark v1, superseded by P10);
`pillar3_dark_matter_97mags.tsv` (untraceable transcription, replaced by the P14 v2 tables); `pillar6_empty_quarter_vignette.tsv` (no inputs or generating script in the repository; its table was
removed from the supplement); `scenario1_cafa_temporal_holdout_real.tsv` (runtimes hard-coded in
`scripts/eval_scenario1_holdouts_real.py`; its panel was removed from Table S4);
`runtime_memory_empirical_scaling.tsv`, `all_8_scenarios_sota_iterative_runtime_memory.tsv`
(transcribed timing archives); `allthebacteria/insight*.tsv` without `_java` (Stage C heuristic; only the alignment columns of `insight3b` and `insight5b` are used).
