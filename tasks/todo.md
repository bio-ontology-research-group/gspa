# Todo: FM operon understanding (gspa)

Companion checklist for `tasks/plan.md`. Strike tasks as they complete.

## Phase 1 — gLM operon caller proven on real data

- [x] **1.A** Sidecar scaffold + Groovy wrapper (landed in branch `phase11-crossgenome`)
- [ ] **1.B** Pixi env + gLM/ESM2 weights on `/mnt/data/u/hohndor/` (ORIX)
- [ ] **1.C** Implement `real_run()` against `y-hwang/gLM` API in `run_glm_operon.py`
- [ ] **1.D** Single-genome sanity check on B. subtilis (recover `trpEDCFBA` etc.)
- [ ] **1.E** `--operon-caller {heuristic,glm}` switch in `run_integrate_full_priors.sh`
- [ ] **1.F** 13-genome manifest + array sbatch, both modes integrated
- [ ] **1.G** F-max + claims-fired delta report → `benchmark/RESULTS.md` addendum

### CHECKPOINT 1 — phase-1 go/no-go

- [ ] `./gradlew test` green
- [ ] All 13 genomes integrated in both modes
- [ ] Mean micro F-max Δ ≥ +0.005
- [ ] No genome regresses by > 0.01
- [ ] User sign-off

→ GO: phase 2. NO-GO: write post-mortem, archive branch, stop.

## Phase 2 — GENOMIC_CONTEXT_FM evidence type (gated on Checkpoint 1)

- [ ] **2.A** `benchmark/glm/phase2/build_glm_centroids.py` — gLM contextualized centroid catalog from the 500-genome KAUST panel (leakage-checked vs. 13-genome benchmark + EQ MAGs)
- [ ] **2.B** `--predictor glm-centroid` mode in `run_neural_predictors.py`
- [ ] **2.C** Wire as `GENOMIC_LANGUAGE_MODEL` claims via `02b_parse_predictors_to_claims.py`
- [ ] **2.D** Default Platt calibration registered + B. subtilis end-to-end check
- [ ] **2.E** 13-genome ablation, all-protein and dark-matter slice
- [ ] **2.F** Phase-2 closeout report + go/no-go on phase 3

### CHECKPOINT 2 — phase-2 go/no-go

- [ ] `./gradlew test` green
- [ ] Dark-matter F-max Δ ≥ +0.01
- [ ] All-protein F-max Δ ≥ 0
- [ ] User sign-off

→ GO: phase 3. NO-GO: keep phase 2 as permanent, drop phase 3, stop.

## Phase 3 — `BF(O, P)` augmentation (gated on Checkpoint 2)

- [ ] **3.A** `benchmark/glm/phase3/build_pathway_operon_corpus.py` — known-pathway operon centroids (MetaCyc-curated)
- [ ] **3.B** Plumb corpus into `IntegrationState` and `DarkMatterSuggester`
- [ ] **3.C** Augment `computeRefinedBayesFactor` with embedding-distance term (opt-in flag)
- [ ] **3.D** Dark-matter ablation on 5 gapseq genomes (precision/recall delta)
- [ ] **3.E** Phase-3 closeout (sensitivity sweep, default flag recommendation)

### CHECKPOINT 3 — phase-3 closeout

- [ ] `./gradlew test` green
- [ ] Dark-matter precision Δ ≥ 0 on 5 gapseq genomes
- [ ] Default flag decision recorded in `RESULTS.md`
- [ ] User sign-off

## Resolved decisions (2026-05-05)

- [x] Gating bars accepted as written for now (revisit if numbers are ambiguous).
- [x] **Catalog reference for 2.A**: 500-genome KAUST panel.
- [x] **Branch + folder layout**: stay on `phase11-crossgenome`; phase-2 work under `benchmark/glm/phase2/`, phase-3 work under `benchmark/glm/phase3/`. ORIX outputs under `/mnt/data/u/hohndor/gspa-glm/phase{1,2,3}/`.

## Still open (do not block — decide when they arise)

- [ ] Stratify glm-centroid calibration by GO aspect (BP / MF / CC) — decide in 2.D from raw-score distribution.
- [ ] Track gLM2 (TattaBio) as an upgrade path — revisit before 3.E.

## Out of scope for this plan

- Fine-tuning gLM on KAUST data
- Migrating to gLM2 (mixed-modality)
- Empty-Quarter MAGs run (consumer of phase 2 / phase 3 output, not part of the gating)
- Methods paper write-up (separate workstream once phases land)
