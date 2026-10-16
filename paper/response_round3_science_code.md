# Point-by-Point Response to Round 1 Technical & Scientific Reviews (Fable 5.1 & GPT Astra Validation Audits)

We thank both reviewers (**Fable 5.1** across `paper/review_round3_fable_science_code.txt`, `paper/review_round3_fable_validation.txt`, `paper/review_round3_fable_validation_v2.txt`, and `paper/review_round3_fable_validation_v3.txt`, and **GPT Astra** across `paper/review_round3_astra_science_code.txt`, `paper/review_round3_astra_validation.txt`, and `paper/review_round3_astra_validation_v2.txt`) for their exhaustive line-by-line code, equation, benchmark, and statistical audits. Every code, mathematical, dataset, provenance, and manuscript issue identified across all validation audits has been resolved, verified by unit tests (`./gradlew test`, 306 tests passing), re-executed across the repository scripts, and synchronized across `results/*.tsv`, `results/PROVENANCE_MAP.md`, `CHANGELOG.md`, `paper/main.tex` (`paper/main.pdf`), and `paper/supplement.tex` (`paper/supplement.pdf`).

---

## Part I. Explicit Item-by-Item Disposition of `review_round3_fable_validation_v3.txt` (Items 1–14)

- **Fable v3 Item 1 (BacDive Track 1 & Track 2 formulas, models, features, and hash splits in `paper/main.tex` Section 2.4 & Section 4.8 vs. `benchmark/unified/pillar3_bacdive/run_bacdive_tracks1_2.py`):**
  - **Resolved:** Synchronized Section 2.4 (`line 291`) and Section 4.8 (`lines 564–571`, `Eq. 8`) in `paper/main.tex` with `benchmark/unified/pillar3_bacdive/run_bacdive_tracks1_2.py`:
    - Documented `build_cooccurrence_coupling`: blended call matrix $B(s, y) = 0.55\,g(s, y) + 0.30\,m(s, y) + 0.15\,c(s, y)$ on carbon and $0.55\,g + 0.28\,m + 0.17\,c$ on enzymes, column-centered nonnegative cosine similarity $\Psi_0(y, y')$ with zero diagonal, row-normalized $\Psi(y, y')$, and unsupervised cross-assay coupling score $\psi(s, y) = \sum_{y'} \Psi(y, y')(B(s, y') - 0.5)$.
    - Documented **GSPA Track 1 (100% Blind)** (`Eq. 8`): calibrated probabilities $p_{\text{gap}} \in \{0.12, 0.88\}$, $p_{\text{ms}} \in \{0.22, 0.80\}$, $p_{\text{cm}} \in \{0.30, 0.72\}$ on carbon ($p_{\text{gap}} \in \{0.14, 0.86\}$, $p_{\text{ms}} \in \{0.24, 0.78\}$, $p_{\text{cm}} \in \{0.28, 0.74\}$ on enzymes) pooled in log-odds space with $(w_{\text{gap}}, w_{\text{ms}}, w_{\text{cm}}, w_\psi) = (0.56, 0.28, 0.16, 0.65)$ on carbon and $(0.55, 0.27, 0.18, 0.55)$ on enzymes.
    - Documented the exact 5-fold genus-blocked CV hash (`int(md5(genus)[:8], 16) mod 5`), baseline $(3+M)$-dimensional feature vector $\mathbf{x}_{\text{base}}(s, y) = [g, m, c, \mathbf{e}_y]$, classifier hyperparameters (`Traitar` `LinearSVC` + `CalibratedClassifierCV` with $C=0.3$ on carbon and $C=0.2$ on enzymes; `BacDive RF` $200/150$ trees at max depth $8$; `DiASPora` `GradientBoostingClassifier(150, 4, 0.05)` on carbon and `HistGradientBoostingClassifier(100, 4, 0.08)` on enzymes), fold-local pooled score $p_{\text{pool}}(s, y)$, $(10+4M)$-dimensional $\mathbf{x}_{\text{GSPA}}(s, y)$ vector, and hybrid blend $P_{\text{Hyb1}} = 0.88\,P_{\text{ML}} + 0.12\,p_{\text{pool}}$.
    - Documented the exact Track 2 30% revealed split (`int(md5(f"{id}_{assay}_{idx}")[:8], 16) mod 100 < 30`), `GSPA Track 2 Bayesian Anchor` ($P_{\text{Track2}} = \sigma_{\text{log}}(0.72\operatorname{logit}(p_{\text{pool}}) + 0.28\operatorname{logit}(s_{b,\text{emp}}) + 0.85\psi_{\text{rev}})$ on carbon and $\sigma_{\text{log}}(0.70\operatorname{logit}(p_{\text{pool}}) + 0.30\operatorname{logit}(g_{e,\text{emp}}))$ on enzymes, thresholds tuned on the revealed 30% split), and `GSPA Track 2 + DiASPora Hybrid` ($P_{\text{Hyb2}} = 0.65\,P_{\text{ML2}} + 0.35\,P_{\text{Track2}}$).

- **Fable v3 Item 2 (`Price-149` responsibility, outer loop Eq. 7, and multi-predictor probability blending in `paper/main.tex` Section 4.6 vs. `scripts/eval_price149_real_all6.py`):**
  - **Resolved:** Synchronized Section 4.6 (`lines 547–553`, `Eq. 7`) in `paper/main.tex` with `scripts/eval_price149_real_all6.py` (`lines 364–465`):
    - Documented $r_j(p) = w_h / (0.05 + w_h)$ where $w_h = \text{odds}_{\text{seq}}^{0.5}\,\exp(0.85\,S_{\text{1hop}} + 0.40\,S_{\text{2hop}})\,\exp(1.2\,C_{\text{syn}})\,g_{\text{need}}$ and $g_{\text{need}} = \max(0.10, 1 - 0.85\,s_{\text{other},\max})$.
    - Documented the recursive outer loop (`Eq. 7`) $L_{\text{outer}}^{(1)} = (1 - 0.70)L^{(0)} + 0.70 L_{\text{target}}$, $L_{\text{outer}}^{(2)} = (1 - 0.42)L_{\text{outer}}^{(1)} + 0.42(L_{\text{target}} + 0.25\Delta_{\text{tot}})$, and $L_{\text{outer}}^{(3)} = (1 - 0.252)L_{\text{outer}}^{(2)} + 0.252(L_{\text{target}} + 0.35\Delta_{\text{tot}})$, and the probability-space multi-predictor blend $P_{\text{multi}}^{(t)}(p, u) = 0.55\,P_{\text{single}}^{(t)}(p, u) + 0.45\,P_{\text{ens}}(p, u)$.

- **Fable v3 Item 3 (`Price-149` dataset macro/micro $F_{\max}$ and TopEC single-predictor narrative in `paper/main.tex` Section 2.3 vs. `results/pillar2_price149_all6_real.tsv`):**
  - **Resolved:** Updated Section 2.3 (`line 224`) in `paper/main.tex` so every dataset macro and micro $F_{\max}$ matches `results/pillar2_price149_all6_real.tsv` exactly (`CLEAN` `0.4150 / 0.5067` $\to$ `0.4216 / 0.4932` single-pred $\to$ `0.4706 / 0.5831` multi-pred; `MAPred` `0.4337 / 0.5284` $\to$ `0.4494 / 0.5439` $\to$ `0.4847 / 0.5877`; `DeepProZyme` `0.4426 / 0.5300` $\to$ `0.4648 / 0.5503` $\to$ `0.4832 / 0.5867`; `EnzBERT` `0.3886 / 0.4561` $\to$ `0.3834 / 0.4362` $\to$ `0.4241 / 0.5112`; `GraphEC` `0.2143 / 0.2466` $\to$ `0.2188 / 0.2432` $\to$ `0.3682 / 0.4351`; `TopEC` `0.0238 / 0.0328` across single-predictor and multi-predictor rows since TopEC emits top-1 calls on only $5$ of $148$ proteins).

- **Fable v3 Item 4 (`Holdout-125` and 95-genome EC homology-masking narrative in `paper/main.tex` Section 2.1 vs. `results/scenario1_swissprot_ec_holdout125_real.tsv` and `results/scenario1_95genomes_ec_homology_real.tsv`):**
  - **Resolved:** Updated Section 2.1 (`lines 132, 136`) in `paper/main.tex` to state explicitly that on `Holdout-125`, `DeepProZyme + GSPA Single-Pred` improves Top-5 accuracy (`32.8%` vs. `32.0%`) and MRR (`0.2779` vs. `0.2750`) while Top-1 accuracy remains `23.2%` (`29/125`) and micro $F_{\max}$ is `0.2703` (vs. `0.2715`), and that on the 95-genome bacterial panel, single-predictor refinement improves CLEAN across all four homology tiers (`+77` to `+139` net Top-1 hits, $p < 10^{-6}$) while multi-predictor refinement improves all three predictors across all four tiers (`+36` to `+155` net Top-1 hits).

- **Fable v3 Item 5 (Hyperparameter sensitivity analysis: 17 configurable parameters in `scripts/run_hyperparameter_sensitivity.py` and `results/hyperparameter_sensitivity_non_ecoli.tsv`):**
  - **Resolved:** Created and executed `scripts/run_hyperparameter_sensitivity.py` over the **17** genuine configurable pipeline parameters (`12` CLI flags in `IntegrateCommand.groovy` and `5` parameters in `scripts/stoichiometric_locus_refiner.py`), regenerated `results/hyperparameter_sensitivity_non_ecoli.tsv`, and updated `paper/main.tex` (`line 506`), `paper/supplement.tex` (Section S1 and Table S1), and `results/PROVENANCE_MAP.md`.

- **Fable v3 Item 6 (LSO complex completion numbers in `paper/main.tex` vs. `CHANGELOG.md`):**
  - **Resolved:** Updated `paper/main.tex` (`lines 67, 81, 184, 411`) to report both the 6-genome unweighted mean (`86.5%` Hit@1, $F_1 = 0.678$) and the 175-case weighted mean (`92.0%` Hit@1, `100.0%` Recall@0.50, $F_1 = 0.692$ in `M3` and `93.1%` Hit@1, $F_1 = 0.728$ in `M4`), and updated `CHANGELOG.md` (`line 36`) to match `results/pillar3_complex_lso.tsv` and `CoherencePrior.groovy` (`+0.85` / `+0.5525`).

- **Fable v3 Item 7 (`Table 1` caption and Stage 3 anomaly breakdown in `paper/main.tex` & `CHANGELOG.md`):**
  - **Resolved:** Updated `Table 1` (`lines 112–113`) in `paper/main.tex` to state $\alpha_2 = 1.0, \alpha_5 = 1.5$ and the out-of-sample 8-genome Stage 1 row (`0.7749` [+0.0029], CAFA $F_{\max}$ `0.7852` [+0.0014], GAEF `0.8599` [+0.1656]), and updated Section 2.5 (`lines 340, 361, 363`) so all promotion/relaxation/demotion/abstention counts match `results/pillar4_resolutions_c3_guarded.jsonl`.

- **Fable v3 Item 8 (108-genome invariant audit in `benchmark/unified/pillar2_fba_lro_price149/run_pillar2_all.py`, `results/pillar2_invariant_audit_108.tsv`, and `results/pillar2_invariant_audit_summary.tsv`):**
  - **Resolved:** Re-executed `run_108_genome_invariant_audit(refiner)` in `benchmark/unified/pillar2_fba_lro_price149/run_pillar2_all.py` with zero magic constants, regenerating `results/pillar2_invariant_audit_108.tsv` and `results/pillar2_invariant_audit_summary.tsv`: METEOR v8 adds `313.35` non-core reactions per genome (`51.01` flux-carrying additions + `262.34` [`83.72%`] zero-flux padding; `51.01 + 262.34 = 313.35`), whereas GSPA SLR adds `5.70` (`5.70` flux-certified, `0.00` uncertified padding, `300.6` flux-carrying active out of `2,815.2` total active, `2.29` FVA-essential, `1.96` KO-essential). Synchronized `paper/main.tex` (`lines 67, 81, 226, 411`), `paper/supplement.tex` (`Table S5`), `paper/generate_figure6.py`, and `CHANGELOG.md`.

- **Fable v3 Item 9 (`results/pillar2_curated6_gems_real.tsv` and `Table S5`):**
  - **Resolved:** Regenerated `results/pillar2_curated6_gems_real.tsv` from `run_pillar2_all.py` containing only the 21 real rows from `panelB_curated6.json`, `panelB_ec.json`, and `recon_dpz_ec.json` (`CarveMe`, `Reconstructor`, `METEOR v8`), and updated `paper/main.tex` (`line 226`) and `paper/supplement.tex` (`Table S5`) accordingly.

- **Fable v3 Item 10 (Section 4 `Methods` code-to-text alignment across items 10.1–10.7):**
  - **10.1 (`GenomicDensityField.groovy:60`):** Updated `Eq. 3` (`line 491`) in `paper/main.tex` to $K_{\text{bp}}(\Delta\text{bp}) = \exp\bigl(-\Delta\text{bp}^2 / (2 h_{\text{bp}}^2)\bigr)$.
  - **10.2 (`StoichiometricLocusRefiner.groovy:83–97`):** Updated Section 4.6 (`line 547`) to state $C_{\text{syn}} = \text{clip}(B_{\text{ctx}}(p, u) + B_{\text{coh}}(p, u), 0, 1)$ and $S_{\text{ESM2}}$ equal to the maximum raw claim score `sc.rawScore` over ESM/ProtT5/PLM/CLEAN supporting claims.
  - **10.3 (`EvidenceCombiner.groovy:51–70` and `EvidenceType.groovy:42–97`):** Updated Section 4.1 (`line 455`) to list the exact `EvidenceCombiner.defaultReliability()` weights (`0.60` for `GENOMIC_REGION_ML`, `SEQUENCE_REGION_ML`, `DARK_MATTER`, `REACTION_LOCAL_CONTEXT`, `CROSS_GENOME_TRANSFER`, and `ML_RANKER`) and exact `EvidenceType.correlationGroup()` strings (`region_features` for `SEQUENCE_REGION_ML` and `viral` for `GENOMIC_REGION_ML`).
  - **10.4 (`GenomicContextPrior.groovy:135–146`):** Updated $B_6$ (`lines 500–502`) to require $P^{(k-1)}(p', f) > 0.50$ for the co-operonic heteromeric complex partner check.
  - **10.5 (`ConsistencyPrior.groovy:98–100`):** Updated $B_3$ (`line 482`) to reference `checker.check(checkSet).minimalUnsatCores`.
  - **10.6 (`AnomalyExtractor.groovy:180–191`):** Updated Section 4.7 (`line 557`) to state that `AnomalyExtractor.groovy` triggers `INCOMPLETE_COMPLEX` when an obligate heteromeric complex (`'n'`) has exactly one active anchor protein ($\ge 0.50$) and $\ge 1$ sub-threshold candidate in $[\theta_{\text{floor}}, 0.50)$.
  - **10.7 (`go-basic.obo`):** Updated Section 4.7 (`line 559`) to state `48,101 GO term stanzas (40,122 active and 7,979 obsolete) in go-basic.obo`.

- **Fable v3 Item 11 (CAFA temporal holdout counts in `paper/main.tex` Section 2.1):**
  - **Resolved:** Clarified in Section 2.1 (`line 126`) and `Table 8` (`line 438`) of `paper/main.tex` that `results/scenario1_cafa_temporal_holdout_real.tsv` evaluates the `1,920`-protein bacterial temporal holdout (`800` proteins `<30%` identity), while `results/suite1_cafa6_prospective_summary.tsv` evaluates the `769`-protein multi-species prospective holdout (`412` NK and `357` LK).

- **Fable v3 Item 12 (Provenance documentation in `results/PROVENANCE_MAP.md`, `paper/supplement.tex` Section S2, `paper/main.tex` Section 4.4, and `results/pillar4_resolutions_c3_guarded.jsonl`):**
  - **Resolved:**
    1. Updated `results/PROVENANCE_MAP.md` to explicitly document which tables are generated by local repository scripts and which tables (`pillar3_eukaryote_endogenous_vs_injected.tsv`, `pillar2_antismash_bgc_106genomes.tsv`, `pillar4_extremophile_benchmark.tsv`, `pillar3_dark_matter_97mags.tsv`, `pillar1_suite2_panel13_21_sota.tsv`, `suite1_cafa6_prospective_summary.tsv`, `all_8_scenarios_sota_iterative_runtime_memory.tsv`, `runtime_memory_empirical_scaling.tsv`) are static transcribed archives of the `unimatrix01` cluster runs.
    2. Updated `paper/supplement.tex` Section S2 and Table S2 (`lines 112, 116`) to cite `scripts/extract_anatomical_and_taxon_constraints.py`, `TaxonConstraints.groovy`, `AnatomicalConstraints.groovy`, and `SatConsistencyChecker.groovy` (`gspa-core/src/main/resources/taxon-constraints/` and `anatomical-constraints/`).
    3. Removed the unarchived 100-pair manual audit sentence from `paper/main.tex` Section 4.4 (`line 511`).
    4. Populated the 3 missing `anomaly_type` fields (`ANOM_0018`, `ANOM_0063`, `ANOM_0076`) in `results/pillar4_resolutions_c3_guarded.jsonl` from `results/pillar4_agent_benchmark_120.jsonl`.

- **Fable v3 Item 13 (Removal of byte-identical duplicate TSVs in `results/`):**
  - **Resolved:** Removed all 8 duplicate TSV files (`pillar1_suite1_cafa_temporal.tsv`, `pillar1_4_suite2_proteomes25_protocol_H_D.tsv`, `pillar3_metabolism_lro.tsv`, `pillar5_agent_vs_rerankers.tsv`, `pillar2_price149_all6_summary.tsv`, `pillar2_price149_all6_baselines.tsv`, `pillar3_eukaryotic_anatomy.tsv`, `pillar2_price149_stratum_real.tsv`) via `git rm -f` and updated all benchmark/figure scripts, `Table 8` in `paper/main.tex`, and `results/PROVENANCE_MAP.md` so only canonical TSV filenames are written and referenced (`0` duplicate MD5 hashes across `results/*.tsv`).

- **Fable v3 Item 14 (Explicit item-by-item disposition in `paper/response_round3_science_code.md`):**
  - **Resolved:** Provided in Parts I, II, and III of this document.

---

## Part II. Explicit Item-by-Item Disposition of `review_round3_fable_validation_v2.txt` (Items 1–24)

1. **Fable v2 #1 (`CalibrationTable.groovy` Platt scaling vs. isotonic regression):** Resolved in `paper/main.tex` Section 2.1 (`line 122`) and Section 4.1 (`line 455`).
2. **Fable v2 #2 (`EvidenceType.groovy` default reliability weights):** Resolved in `paper/main.tex` Section 4.1 (`line 455`), matching `EvidenceCombiner.defaultReliability()`.
3. **Fable v2 #3 (`EvidenceCombiner.groovy` correlation-group Noisy-OR formula Eq. 1):** Resolved in `paper/main.tex` Section 4.1 (`lines 455–460`, `Eq. 1`).
4. **Fable v2 #4 (`IterativeRefiner.groovy` log-odds clamping and damping Eq. 2):** Resolved in `paper/main.tex` Section 4.2 (`lines 464–472`, `Eq. 2`).
5. **Fable v2 #5 (`EssentialityPrior.groovy` and `CoherencePrior.groovy` formulas $B_1, B_2$):** Resolved in `paper/main.tex` Section 4.2 (`lines 474–481`).
6. **Fable v2 #6 (`GapFillingPrior.groovy`, `ReactionLocalContextSuggester.groovy`, `GenomicDensityField.groovy` Eq. 3):** Resolved in `paper/main.tex` Section 4.2 (`lines 491–499`, `Eq. 3`).
7. **Fable v2 #7 (`GenomicContextPrior.groovy` $B_6$ and `ConsistencyEnforcer.groovy` Stage 2A Eq. 4):** Resolved in `paper/main.tex` Section 4.2 (`lines 500–503`) and Section 4.4 (`lines 515–521`, `Eq. 4`).
8. **Fable v2 #8 (`StoichiometricLocusRefiner` 108-genome certificate audit):** Resolved in `benchmark/unified/pillar2_fba_lro_price149/run_pillar2_all.py`, `results/pillar2_invariant_audit_108.tsv`, `results/pillar2_invariant_audit_summary.tsv`, and `paper/main.tex` Section 2.3 & Section 4.5.
9. **Fable v2 #9 (`attribute_locus_responsibility` Eq. 6 in `scripts/stoichiometric_locus_refiner.py`):** Resolved in `paper/main.tex` Section 4.5 (`lines 540–547`, `Eq. 6`).
10. **Fable v2 #10 (`eval_price149_real_all6.py` Eq. 7 outer loop):** Resolved in `paper/main.tex` Section 4.5 (`lines 547–553`, `Eq. 7`).
11. **Fable v2 #11 (`Qwen/Qwen3.8-27B-FP8` model identifier in `run_gspa_agent.py`):** Resolved in `paper/main.tex` Section 4.6 (`line 559`), `paper/supplement.tex` (`line 257`), and `CHANGELOG.md` (`line 30`).
12. **Fable v2 #12 (`evaluate_pillar4.py` 3-feature and 8-feature tabular ML baselines):** Resolved in `paper/main.tex` Section 2.5 (`line 340`).
13. **Fable v2 #13 (`run_bacdive_tracks1_2.py` Track 1 & Track 2 equations):** Resolved in `paper/main.tex` Section 4.7 (`lines 564–571`, `Eq. 8`).
14. **Fable v2 #14 (`pillar1_2_ablation_per_genome.tsv` per-genome essentiality and complex coherence deltas):** Resolved in `paper/main.tex` Section 2.1 (`line 122`).
15. **Fable v2 #15 (`multikingdom_per_genome_eval.tsv` 88-genome 95% CIs by domain):** Resolved in `paper/main.tex` Section 2.1 (`line 126`).
16. **Fable v2 #16 (`pillar1_suite2_1000bacteria_sota.tsv` Iter 3 macro $F_{\max}$ `0.6937`):** Resolved in `results/pillar1_suite2_1000bacteria_sota.tsv` and `paper/main.tex` Section 2.1 (`line 126`).
17. **Fable v2 #17 (`multikingdom_summary_by_domain.tsv` runtime and heap memory synchronization):** Resolved in `results/multikingdom_summary_by_domain.tsv` and `paper/main.tex` Section 2.1 (`line 126`).
18. **Fable v2 #18 (`suite1_cafa6_prospective_summary.tsv` aspect counts and NK/LK split):** Resolved in `paper/main.tex` Section 2.1 (`line 126`).
19. **Fable v2 #19 (`scenario1_95genomes_ec_homology_real.tsv` single-pred vs. multi-pred separation):** Resolved in `paper/main.tex` Section 2.1 (`line 136`) and `paper/supplement.tex` (`Table S4`).
20. **Fable v2 #20 (`pillar3_eukaryote_endogenous_vs_injected.tsv` Iter 1 vs. Iter 3 per-species deltas):** Resolved in `paper/main.tex` Section 2.4 (`line 277`).
21. **Fable v2 #21 (`AnatomicalPrior.groovy` $-\alpha_{\text{anat}} = -3.5$ and $\tilde{q}(c) \ge 0.95$):** Resolved in `paper/main.tex` Section 2.4 (`line 277`) and Section 4.2 (`line 483`).
22. **Fable v2 #22 (`pillar4_resolutions_c3_guarded.jsonl` seed-42 anomaly outcome counts):** Resolved in `paper/main.tex` Section 2.5 (`lines 340, 361, 363`).
23. **Fable v2 #23 (`paper/supplement.tex` Table S1, Table S2, Table S5 synchronization):** Resolved in `paper/supplement.tex` (`Tables S1, S2, S5`).
24. **Fable v2 #24 (`results/PROVENANCE_MAP.md` and `CHANGELOG.md` synchronization):** Resolved in `results/PROVENANCE_MAP.md` and `CHANGELOG.md`.

---

## Part III. Explicit Item-by-Item Disposition of `review_round3_astra_validation_v2.txt` (Items 1–18)

1. **Astra v2 #1 (Source-code package paths in `paper/main.tex` Section 4):** Resolved across Section 4 (`gspa.integration.*`, `gspa.integration.prior.*`, `gspa.integration.fba.*`, `gspa.config.EssentialFunctions`).
2. **Astra v2 #2 (`run_bacdive_tracks1_2.py` Track 1 & Track 2 mathematical specification):** Resolved in `paper/main.tex` Section 4.7 (`lines 564–571`, `Eq. 8`).
3. **Astra v2 #3 (`CalibrationTable.groovy` Platt scaling parameters):** Resolved in `paper/main.tex` Section 4.1 (`line 455`).
4. **Astra v2 #4 (`EvidenceCombiner.groovy`, `IterativeRefiner.groovy`, and `ConsistencyEnforcer.groovy` formulas):** Resolved in `paper/main.tex` Sections 4.1, 4.2, and 4.4 (`Eqs. 1, 2, 4`).
5. **Astra v2 #5 (`SparseRoleAggregator` and `attribute_locus_responsibility` in `stoichiometric_locus_refiner.py`):** Resolved in `paper/main.tex` Section 4.5 (`Eqs. 5–6`).
6. **Astra v2 #6 (108-genome FVA/KO certificate execution and active vs. flux-carrying reaction counts):** Resolved in `benchmark/unified/pillar2_fba_lro_price149/run_pillar2_all.py`, `results/pillar2_invariant_audit_108.tsv`, `results/pillar2_invariant_audit_summary.tsv`, and `paper/main.tex` Section 2.3 & Section 4.5.
7. **Astra v2 #7 (`Qwen/Qwen3.8-27B-FP8` checkpoint name and 3/8-feature ML baseline description):** Resolved in `paper/main.tex` Section 2.5 (`line 340`) and Section 4.6 (`line 559`).
8. **Astra v2 #8 (BacDive 5-fold genus-blocked CV hash and Track 2 30% split):** Resolved in `paper/main.tex` Section 4.7 (`lines 564–571`).
9. **Astra v2 #9 (Table 1 per-genome essentiality and complex coherence deltas):** Resolved in `paper/main.tex` Section 2.1 (`line 122`).
10. **Astra v2 #10 (88-genome multi-kingdom 95% CIs across domains):** Resolved in `paper/main.tex` Section 2.1 (`line 126`).
11. **Astra v2 #11 (`1000-bacteria` Iter 3 macro $F_{\max}$ `0.6937`):** Resolved in `results/pillar1_suite2_1000bacteria_sota.tsv` and `paper/main.tex` Section 2.1 (`line 126`).
12. **Astra v2 #12 (`multikingdom_summary_by_domain.tsv` runtime and memory synchronization):** Resolved in `results/multikingdom_summary_by_domain.tsv` and `paper/main.tex` Section 2.1 (`line 126`).
13. **Astra v2 #13 (`suite1_cafa6_prospective_summary.tsv` aspect counts):** Resolved in `paper/main.tex` Section 2.1 (`line 126`).
14. **Astra v2 #14 (`pillar3_eukaryote_endogenous_vs_injected.tsv` Iter 1 vs. Iter 3 deltas):** Resolved in `paper/main.tex` Section 2.4 (`line 277`).
15. **Astra v2 #15 (`AnatomicalPrior.groovy` $-\alpha_{\text{anat}} = -3.5$):** Resolved in `paper/main.tex` Section 2.4 (`line 277`) and Section 4.2 (`line 483`).
16. **Astra v2 #16 (`Price-149` single-predictor vs. multi-predictor table and text alignment):** Resolved in `paper/main.tex` Section 2.3 (`line 224`) and `Table 3`.
17. **Astra v2 #17 (`pillar4_resolutions_c3_guarded.jsonl` seed-42 outcome counts):** Resolved in `paper/main.tex` Section 2.5 (`lines 340, 361, 363`).
18. **Astra v2 #18 (`paper/supplement.tex`, `results/PROVENANCE_MAP.md`, and `CHANGELOG.md` synchronization):** Resolved across `paper/supplement.tex`, `results/PROVENANCE_MAP.md`, and `CHANGELOG.md`.
