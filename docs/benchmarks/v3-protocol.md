# GSPA v3 production-benchmark protocol

Status: foundation (runner, metrics, manifest schema, one real manifest). No large dataset has been executed with it yet. Nothing here modifies the v2 tables in `results/` or `paper/`.

## Principles

1. **The production engine makes every prediction.** `benchmark/v3/v3_runner.py` invokes `gspa integrate --batch-manifest` and contains no inference formula. A benchmark arm is a set of `gspa integrate` flags, nothing else. Python re-implementations of GSPA stages (the P08 Noisy-OR/METEOR-nudge arms, the AllTheBacteria heuristic) are not admissible as v3 arms.
2. **Labels never enter inference.** Labels are opened only in the evaluation phase, after every arm has finished. The runner rejects a manifest whose label file is also a claims file, appears in any `gspa` argument, or whose claims records carry label-like keys (`label`, `truth`, `is_true`, `y_true`, ...).
3. **Splits are validated before anything runs.** A leakage finding aborts the run (`status: aborted: leakage`) and is written to `report.json`.
4. **Failures are data.** Each (arm, genome) cell is `executed`, `skipped: <reason>` or `failed: <reason>`; engine errors never raise, and `--strict` turns any non-executed cell into a non-zero exit.
5. **Everything is hashed.** `report.json` records SHA-256 of the manifest, runner and metric code, git commit, gspa binary, every claims/label/aux input, reference files, batch manifests, and every integrated output.
6. **Descriptive vs deployable metrics are labelled.** Global Fmax chooses its threshold on the evaluated labels (`threshold_on_test: true`): it is a descriptive upper bound. Fixed-threshold precision/recall, calibration and correction/harm use a threshold fixed in the manifest before the run (default 0.5).

## Running

```bash
export GSPA_JAR=/data/hohndor/gspa/bin/gspa-unified.jar   # or a freshly built shadow jar
export GSPA_DATA=/data/hohndor/gspa/proteomes GSPA_REF=/data/hohndor/gspa/reference
python3 benchmark/v3/v3_runner.py validate benchmark/v3/manifests/v3_panel13_go.json
python3 benchmark/v3/v3_runner.py run      benchmark/v3/manifests/v3_panel13_go.json --out out/v3-panel13 --strict
python3 benchmark/v3/v3_runner.py evaluate benchmark/v3/manifests/v3_panel13_go.json --out out/v3-panel13   # re-score only
python3 -m unittest discover -s benchmark/v3/tests -v
```

`--dry-run` writes the batch manifests and exact commands without executing. Unset `${VAR}` placeholders stay literal, so their files are reported `missing`/`skipped` instead of silently resolving. For a Gradle-launched engine set `"args_style": "gradle"` and `"command": ["./gradlew", "-q", ":gspa-cli:run"]`. The opt-in smoke test (`GSPA_V3_ENGINE_SMOKE=1`) does this against the real CLI.

## Manifest (`schema_version: 1`)

| Key | Meaning |
|---|---|
| `experiment_id`, `description` | free text identifiers |
| `dataset_status` | `real` or `synthetic`. Required; synthetic results must never be pooled with real ones in a table |
| `gspa.command`, `binary`, `shared_args`, `args_style`, `cwd`, `timeout_sec` | how to launch the engine; `binary` is hashed |
| `arms[]` | `{name, args}`; `args` are extra `gspa integrate` flags |
| `baseline_arm` | the matched baseline all other arms are compared with |
| `genomes[]` | `genome_id`, `split` (`train`/`calibration`/`test`), `claims`, optional `labels: {go, ec}`, `taxon`, `essential_profile`, `operons`, `gaps`, `species`, `cluster`, `train_protein_ids` |
| `train_proteins`, `homology_clusters` | optional global training-id list and `protein<TAB>cluster` map |
| `reference` | hashed reference files; `go_obo` enables is_a/part_of propagation of GO truth and predictions |
| `settings` | `threshold` (on the 0.01..0.99 grid), `bootstrap` replicates, `seed` |

Only `test` genomes are run and scored. `train`/`calibration` genomes exist so that overlap can be checked against them.

## Leakage controls (implemented in `leakage_problems`)

* No genome id, species or declared cluster is shared between fit (`train`, `calibration`) and `test` genomes.
* No protein id in a test genome's claims occurs in the training set: union of `train_proteins`, each `train_protein_ids`, and the claims of fit genomes. This is deliberately conservative: because protein ids are FASTA seqids, a collision across genomes is reported even when it is coincidental, and must be resolved by prefixing ids.
* With `homology_clusters`, no cluster may contain both a training and a test protein.
* Label/claims separation as in principle 2.

Not enforced by code, and required of every manifest author: (a) predictor channels inside the claims must themselves be homology-masked or leave-target-out (e.g. DIAMOND against a reference with the benchmark proteins removed); (b) hyperparameters (`--theta`, calibration tables) must be frozen before looking at test labels and their provenance listed in `description`; (c) no training on test labels. The v2 calibration genome (*E. coli* K-12) must stay out of test sets that claim out-of-sample results.

## Metrics (`v3_metrics.py`)

All metrics are computed on the **labelled universe**: proteins with at least one label of the scored kind. Predictions for other proteins are ignored because a missing label is not a negative.

* **Complete-EC top-1.** For each labelled protein, the highest-scoring *complete* EC (four numeric levels; `1.1.1.-` and `3.5.1.n2` cannot be top-1) is correct if in the protein's true EC set. Ties go to the lexicographically smaller EC. Proteins with no complete-EC prediction count as not correct; coverage and precision-when-called are reported separately.
* **Fixed-threshold precision / recall / F1.** Micro over (protein, term) pairs, call = `score >= tau`, recall denominator = all true pairs.
* **Global Fmax.** Pairs pooled over all genomes, sweep over 0.01..0.99, ties to the lowest threshold. Descriptive (threshold selected on test labels).
* **Calibration.** ECE (10 equal-width bins), Brier score, reliability table over predicted pairs on the labelled universe. Under-reports calibration for terms outside the labelled universe, so read it with the ontology caveat below.
* **Correction / harm against the matched baseline.** Over the pairs called by either arm plus all true pairs: *corrected* = baseline wrong and arm right, *harmed* = baseline right and arm wrong, *net* = difference, exact McNemar p. For EC top-1 the same counts per protein.
* **Genome bootstrap.** Genomes are resampled with replacement (default 2000 replicates, seed 13, percentile 95 % CI). Paired deltas use identical resamples for arm and baseline. With fewer than 2 genomes the CI is NaN, not zero-width.
* **GO propagation.** If `reference.go_obo` is set, truth and predictions are propagated over `is_a` and `part_of`, ancestors take the maximum descendant score, alt ids are mapped, obsolete terms and the three roots are dropped. Without it, GO scoring is exact-term and must be labelled as such.

Known limits: GO truth from GOA is mostly IEA for all-evidence sets (about 85 % in the P03 panel). Use experimental-evidence truth as the primary label and report all-evidence only as secondary. Experimental truth is small (323 proteins across the 13-genome panel), so CIs are wide and many genomes contribute very little.

## Inventory of established datasets (checked against the repository, 2026-10-26)

"Real" = labels come from an external experimental or curated source and the repository holds or names the files. "Present" below means a path exists in this worktree; cluster files were not checked (no remote access in this task).

| Benchmark | Data / manifest | Status | What is missing for a v3 run |
|---|---|---|---|
| **GO reference panels (13 and 21 genomes)** | P03: claims `${GSPA_DATA}/bench10/claims/*_dpi_merged.jsonl`, `bench9/with_iprscan/*_claims_merged.jsonl`; experimental truth `bench10/truth/*_truth_exp.tsv`, `truth_dual/*_truth_exp.tsv`; hashes in `paper-reproducibility/manifest/data_manifest.tsv`. **Manifest: `benchmark/v3/manifests/v3_panel13_go.json`** | Real (UniProt/GOA). Cluster-resident | Run on unimatrix01. The 21-genome panel (`/data/hohndor/gspa-neural/panel`) needs a second manifest. Taxon for 5 genomes is a higher-rank ancestor (hierarchy lacks the species) |
| **Price-148/149 enzyme identity correction (RB-TnSeq fitness-validated ECs)** | `benchmark/manifests/suite6_price149_manifest.tsv`; six predictors' raw EC matrices in `benchmark/manifests/reference_data/ibex_extracted/price149_all6_compact.pkl.gz` (P08); `scripts/eval_price149_real_all6.py` | Real labels. 148 of 149 evaluable (`NP_384885` carries a replicon id) | The v2 arms are Python re-implementations, so they do not qualify. Needs an adapter that writes the six predictors' scores as `EC` claims JSONL (one claims file per host assembly, ±10-neighbour context from the layout) and then `gspa integrate`. TopEC is single-call (one EC per protein in the raw files); GraphEC/TopEC must not be credited with pool gains |
| **CAFA prospective / homology-tiered panels** | `benchmark/manifests/suite1_cafa_temporal_manifest.tsv`; `reference_data/ibex_extracted/cafa_temporal_holdout_real.json`; `scripts/run_cafa_whole_genome_splits.py` | Real temporal holdout. Validation report items 8, 9, 27: stale-alignment calibration bug, undisclosed `hom_scale`, undefined t0 / GOA release | Re-run the six splits with a corrected validation alignment (orix) and write claims per split. Define and record t0 and the GOA release. Remove or justify `hom_scale` |
| **Curated operons and complexes** | P04: `docs/benchmarks/complex-leave-subunit-out.md` cites Complex Portal and CORUM; inputs on unimatrix01 (`P04_inputs.sha256.tsv`) | Mixed. Complex membership is curated; the v2 inputs are synthetic attenuated claims (validation item 14) | No operon gold standard (e.g. experimentally validated operon set) was located in the repository. Need within-operon decoys and a real-claims variant; the v2 `real_diamond_channel` arm is the starting point |
| **Eukaryote constraints (anatomy/taxon)** | P05 (`--input-source real`), `docs/benchmarks/eukaryotic-anatomy.md` | Real Swiss-Prot labels; the stress cohort is injected (synthetic by construction) | Keep the injection and the real-prediction cohort as separately labelled `dataset_status` experiments; real predictors needed for the stress cohort |
| **BacDive phenotypes** | `benchmark/manifests/suites7_8_bacdive_manifest.tsv`; `reference_data/ibex_extracted/gapseq_carbon_source_2332_real.tsv`, `gapseq_bacdive_enzyme_17893_real.tsv` | Real labels (BacDive); inputs are the gapseq 2021 supplementary tables, not reruns of CarveMe / ModelSEED / gapseq | Not an EC/GO claims benchmark; needs a phenotype scorer on top of `gspa integrate` pathway output. Strain-filtering counts (529 / 3,017) unexplained (item 18). No separate tuning set exists, so no blind claim |
| **Growth / knockout (fitness)** | Price RB-TnSeq (above); `benchmark/unified/pillar2_fba_lro_price149/` | Real | Stage 2B is a Python sidecar (`scripts/stoichiometric_locus_refiner.py`) needing the ModelSEED matrices (gitignored); `--fba-data-dir` bundle must be deposited. No Keio-style knockout/growth gold set was located |
| **Community cross-feeding** | Code only: `CrossfeedingAnalyzer`, `CommunityPredictor` | No benchmark data located | A labelled community dataset (measured cross-feeding or co-culture growth) must be chosen and registered |
| **Secondary metabolism** | P13 (antiSMASH 8, 108 genomes), `docs/benchmarks/secondary-metabolism-mags.md` | MIBiG/antiSMASH real; the 106-BGC and 97-MAG tables were removed as unsupported (item 4) | Regenerate with committed scripts; antiSMASH output is a tool-derived reference, not curated truth, unless restricted to MIBiG |
| **AllTheBacteria** | | **Excluded by instruction** | n/a. P11 and Stage C heuristics are not used |

Synthetic vs real rule: any experiment whose claims or labels are generated by code in this repository (injected perturbations, attenuated claims, Swiss-Prot-derived "predictor" scores) is `dataset_status: synthetic`, is reported in its own table, and is never the sole support for a claim.

## Required baseline and ablation matrix

Every real experiment reports the following arms, each as `gspa integrate` flags, on identical genomes and claims, with the baseline matched:

| Arm | Flags (relative to baseline) | Question |
|---|---|---|
| `stage0` (baseline) | `--enable-priors none` | pooled evidence only |
| `stage1` | default priors | whole-genome priors |
| `stage1_minus_<prior>` ×6 | `--enable-priors` without one of essentiality, coherence, consistency, anatomical, gap_filling, genomic_context | contribution of each prior |
| `stage1_stage2a_relax` / `_remove` / `_minimal-flip` | `--enforce-stage2 --enforce-consistency-mode {relax,remove,minimal-flip}` | repair modes (`relax` vs deletion) |
| `stage1_rlgc` | `--rlc-suggester` (+ catalog for M2) | local-context rescue |
| `stage1_stage2b` | `--enable-fba-refiner --fba-data-dir ...` | flux-certified repair |
| `external_best_single` | the strongest single channel's claims alone through `--enable-priors none` | pooling gain |
| `perturbation_control` | the same claims with injected cross-clade pairs, `dataset_status: synthetic` | survival of injected errors (separate table) |

Also required: a **no-op control** (baseline against itself must give zero corrections/harms and delta CIs of exactly 0, asserted in tests) and a **seed/replicate check** (re-run the bootstrap with a second seed; CIs must agree to the reported precision).

### Adaptive FBA agent versus simple planners (equal budgets)

Stage 3 is not admissible as an improvement unless it beats simple planners under the **same budget**. The comparison uses the same anomaly list (`--anomalies-out`, produced by `gspa integrate` itself, never a hand-built manifest; see validation item 6) and the same per-genome budget of oracle calls: `B` FBA solves and `B` LLM tokens/steps.

| Planner | Rule |
|---|---|
| `agent_guarded` | `--enable-agent` (model, seeds, and temperature recorded; transcripts logged) |
| `agent_unguarded` | `--enable-agent --agent-no-guardrail` (ablation only) |
| `rank_by_posterior` | take the top-`B` anomalies by Stage 1 final log-odds |
| `rank_by_single_feature` | top-`B` by ESM-2 centroid cosine / best single channel (validation item 7: this alone reached Hit@1 96 % on LRO-50) |
| `random` | `B` anomalies at random, ≥ 100 seeds |
| `greedy_fba_gain` | repeatedly pick the candidate with the largest FBA objective gain until `B` solves are spent |
| `stage2b_only` | no agent, deterministic Stage 2B |

Report Hit@k/MRR or pair-level correction and harm per planner, with genome-bootstrap CIs and exact McNemar against `rank_by_single_feature`. The budget (solves, tokens, wall time) is logged per planner. An agent result is reported only if it exceeds the best non-agent planner at equal budget; otherwise it is reported as a negative result.

## Blockers for a full production execution

1. **Engine availability.** A jar built from this branch (the manifest points at `${GSPA_JAR}`; the paper's `gspa-unified.jar` is v2.1.0, which predates v2.2.x rescue changes). Using a different jar than the code under evaluation must be stated in `description`.
2. **Cluster data.** `${GSPA_DATA}` / `${GSPA_REF}` data live on unimatrix01 (and the CAFA splits on orix). Absent locally; only hashes in `data_manifest.tsv` are available, and they were not re-verified here.
3. **Adapters not written.** Price-149 claims export, CAFA split claims, Complex Portal/CORUM LSO with real claims, BacDive phenotype scorer, community cross-feeding scorer.
4. **Labels not located:** operon gold standard, knockout/growth gold set, cross-feeding gold set.
5. **Frozen hyperparameters.** `--theta` and calibration tables for the production arms must be fixed and hashed before the test labels are read; the P03/P02 values changed during the v2 audit (items 13, 25).
6. **Stage 2B bundle** (ModelSEED matrices) and **Stage 3 model endpoint** (Qwen3.8-27B-FP8 on 2× RTX 4090 per item 21) for the agent arms.
7. **The twelve audited recommendations** were not in the task brief. The matrix above covers the validation report's open "To do" items; the mapping to the twelve must be supplied.
8. **Not yet verified on real data:** the integrated TSV parsing and batch layout were checked against the real CLI on a 3-claim fixture only; `v3_panel13_go.json` has not been executed.
