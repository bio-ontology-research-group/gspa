# Implementation Plan: FM-based operon understanding for gspa

Companion to `SPEC.md`. Status: **Draft, awaiting approval.**

## Overview

Land foundation-model (gLM) operon understanding in three gated phases. Phase 1 swaps the intergenic-distance operon detector for a gLM caller and proves the FM signal on the 13-genome benchmark. Phase 2 — gated on phase-1 F-max — emits gLM-derived per-protein function predictions as a new evidence channel into Phase 7 Noisy-OR integration. Phase 3 — gated on phase-2 — augments the Phase 8 dark-matter `BF(O, P)` with embedding-space distance to known-pathway operons.

Phase 1's scaffolding (sidecar, wrapper, tests) has already landed; what remains is the real-inference body, ORIX environment setup, the actual benchmark run, and the go/no-go report. The plan below covers all three phases with explicit checkpoints between them.

## Architecture decisions

0. **One branch, per-phase subfolders.** Stay on `phase11-crossgenome`. Phase-1 artifacts use existing locations (`benchmark/neural/run_glm_operon.py`, the Groovy wrapper under `gspa-predictors/.../context/`). Phase 2 and phase 3 each get a dedicated `benchmark/glm/phase{2,3}/` directory for catalog builders, ablation drivers, and reports — keeps the three gated phases visually and logistically separable while still on one branch. ORIX outputs mirror the layout under `/mnt/data/u/hohndor/gspa-glm/phase{1,2,3}/`.
1. **`GLMOperonPredictor extends OperonPredictor`** (already implemented). Only `detectOperons` is overridden; `transferAnnotations`, `transferScore`, `minOperonSize`, and the metadata schema are inherited so the phase-1 ablation isolates one variable: which proteins are co-operonic. Confidence is a pre-filter on operons, never a per-claim weight.
2. **Sidecar over JNI / ONNX.** The Python sidecar pattern (`benchmark/neural/run_*.py`) is established and the JVM wrapper shells out via `ProcessBuilder`. No native Java for gLM.
3. **Storage on ORIX under `/mnt/data/u/hohndor/`.** Weights and outputs live there per `orix-workbench/README.md`. Local repo holds only scripts and small fixtures; no compute on login nodes; PI partition `pi-hohndor` for protected runs, `freecycle` (`--gres=gpu:h200:1`) for opportunistic.
4. **Reuse `EvidenceType.GENOMIC_LANGUAGE_MODEL`** (already defined, correlation group `ml_genomic`). No enum / schema migration in phases 2 or 3.
5. **For the GENOMIC_CONTEXT_FM channel (phase 2): start with centroid-kNN.** Mirror the existing `build_esm2_centroids.py` and `run_neural_predictors.py --predictor esm2-centroid` pattern, but using gLM contextualized embeddings. Defer trained classifier heads until centroid-kNN is shown to be worth it.
6. **For phase-3 `BF(O, P)` augmentation: distance to a known-pathway operon corpus.** Build the corpus once from MetaCyc-grade reference genomes; cache centroids; lookup is a small NPZ load.
7. **Vertical slicing.** Each task delivers one complete path end-to-end (data → model → integration → metric). No "build all the data, then all the model, then all the wiring" horizontal phases.
8. **Cache embeddings, do not regenerate.** Phase 1 emits both ESM2 protein embeddings and gLM contextualized embeddings; phase 2 consumes them without re-running gLM.

## Dependency graph

```
Phase 1  (gLM operon caller)
  ├─ 1.A Sidecar scaffold + wrapper                 [DONE — landed last turn]
  ├─ 1.B ORIX env (pixi recipe + gLM/ESM2 weights)         ──┐
  │                                                          ▼
  ├─ 1.C real_run() implementation against y-hwang/gLM API   ──┐
  │                                                          ▼
  ├─ 1.D Single-genome sanity check (B. subtilis on ORIX)    ──┐
  │                                                          ▼
  ├─ 1.E `--operon-caller {heuristic,glm}` switch in
  │      run_integrate_full_priors.sh                        ──┐
  │                                                          ▼
  ├─ 1.F 13-genome manifest + array sbatch                   ──┐
  │                                                          ▼
  └─ 1.G Report F-max + claims-fired delta + go/no-go    ◄──── CHECKPOINT 1

Phase 2  (GENOMIC_CONTEXT_FM evidence type)
  ├─ 2.A gLM centroid catalog builder (Swiss-Prot)           ──┐
  │                                                          ▼
  ├─ 2.B Per-protein gLM-centroid sidecar predictor          ──┐
  │       (reads protein_embeddings.npz from phase-1)        ▼
  ├─ 2.C Wire as GENOMIC_LANGUAGE_MODEL claims via
  │      02b_parse_predictors_to_claims.py                   ──┐
  │                                                          ▼
  ├─ 2.D Calibration + small test on B. subtilis             ──┐
  │                                                          ▼
  ├─ 2.E 13-genome ablation, dark-matter slice               ──┐
  │                                                          ▼
  └─ 2.F Report F-max delta + dark-matter recall delta   ◄──── CHECKPOINT 2

Phase 3  (BF(O, P) augmentation in DarkMatterSuggester)
  ├─ 3.A Known-pathway operon corpus (MetaCyc genomes)       ──┐
  │                                                          ▼
  ├─ 3.B Pass corpus into IntegrationState / DarkMatterSuggester ┐
  │                                                          ▼
  ├─ 3.C Augment computeRefinedBayesFactor with               ──┐
  │      centroid-distance term (opt-in flag)                ▼
  ├─ 3.D Ablation: dark-matter suggestions delta             ──┐
  │                                                          ▼
  └─ 3.E Report singleton/disjunctive precision change   ◄──── CHECKPOINT 3
```

Tasks 1.B and 1.C run in parallel where possible (1.B is a one-time env setup; 1.C is the inference body); 1.D depends on both.

---

## Task list

### Phase 1: Step-1 finishing — gLM operon caller proven on real data

#### Task 1.A: Sidecar scaffold + wrapper [DONE]

**Status:** Already landed. `benchmark/neural/run_glm_operon.py`, `gspa-predictors/.../GLMOperonPredictor.groovy`, `GLMOperonPredictorSpec.groovy`, `run_glm_operon.sbatch`, updated `.gitignore`. `./gradlew test` green; `--mode self-test` passes.

#### Task 1.B: Pixi env + gLM / ESM2 weights on ORIX

**Description:** Provision a reproducible Python environment on ORIX with `torch`, `numpy`, `gLM` (vendored from `y-hwang/gLM` at a pinned commit), and `fair-esm` for the ESM2-650M dependency. Cache weights under `/mnt/data/u/hohndor/`.

**Acceptance criteria:**
- [ ] `envs/glm.pixi.toml` (or equivalent) committed under `gspa/orix-workbench` or `gspa/envs/` describing the env.
- [ ] gLM repo cloned to `/mnt/data/u/hohndor/gLM/repo/` at a recorded commit SHA; checkpoint at `/mnt/data/u/hohndor/gLM/weights/`.
- [ ] ESM2-650M cache populated under `/mnt/data/u/hohndor/esm2/` (TORCH_HOME).
- [ ] On a `pi-hohndor` H200 node, `python -c "import glm; import esm; import torch; print(torch.cuda.is_available())"` returns `True`.

**Verification:**
- [ ] Submit a 5-minute interactive job: `srun --partition=pi-hohndor --gres=gpu:1 --pty bash` → `pixi shell` (or `conda activate gspa-glm`) → run the import check.
- [ ] Manual: weights directory `du -sh` < 10 GB and not empty.

**Dependencies:** None (parallel with 1.C scaffolding, but 1.C cannot be tested without this).

**Files likely touched:**
- `envs/glm.pixi.toml` (new) or `orix-workbench/envs/glm.pixi.toml`
- `tasks/notes/glm-env-pin.md` (new — records commit SHAs)

**Estimated scope:** S (1-2 files, but a few hours of cluster fiddling).

---

#### Task 1.C: Implement `real_run()` against the gLM API

**Description:** Replace the `NotImplementedError` in `benchmark/neural/run_glm_operon.py:real_run` with the actual gLM forward-pass: ESM2-650M per-protein embedding → gLM tokenization (esm_emb ⊕ intergenic-distance bin ⊕ strand) → forward over each contig in `context_window`-sized chunks → extract per-pair "operon boundary" probability and per-gene contextualized embeddings.

**Acceptance criteria:**
- [ ] `--mode real --weights /mnt/data/u/hohndor/gLM/weights --fasta tiny.faa --gff tiny.gff` exits 0 on a 5-protein fixture.
- [ ] `pair_break_prob` is non-`None` for adjacent same-strand same-contig pairs and `None` at hard breaks.
- [ ] `glm_embeddings.shape == (n_genes, d_ctx)` where `d_ctx` matches the gLM checkpoint config.
- [ ] `esm_embeddings.shape == (n_genes, 1280)`.
- [ ] No silent fall-through to mock — explicit error if weights cannot be loaded.

**Verification:**
- [ ] `python3 run_glm_operon.py --mode real --fasta test/tiny.faa --gff test/tiny.gff --weights /mnt/.../weights ...` — exits 0, schema valid.
- [ ] Sanity: outputs differ from `--mode mock` on the same input (proves model is doing something).

**Dependencies:** 1.B.

**Files likely touched:**
- `benchmark/neural/run_glm_operon.py` (`real_run` body only)
- `benchmark/neural/test/tiny.faa`, `test/tiny.gff` (new fixture, ~5 proteins)

**Estimated scope:** M (1 file body + tiny fixture; the work is reading gLM source).

---

#### Task 1.D: Single-genome sanity check on B. subtilis (ORIX)

**Description:** Run the full sidecar in `--mode real` on B. subtilis 168 and confirm operon distribution is biologically plausible (mean operon size 3–8 genes; well-known operons like `trpEDCFBA`, `hisGDCBHAFI` recovered).

**Acceptance criteria:**
- [ ] Sidecar runs to completion in < 15 min on 1× H200.
- [ ] `operons.tsv` schema valid (≥2 IDs/line, ≥ 100 operons emitted on B. subtilis ~4 k proteins).
- [ ] Mean confidence on emitted operons > 0.5; distribution non-degenerate (not all 1.0 or all 0.5).
- [ ] At least 3 known operons recovered (manually checked).

**Verification:**
- [ ] `wc -l operons.tsv` and `awk` for mean size; eyeball confidence distribution histogram.
- [ ] Cross-check `trpEDCFBA` (BSU22730–BSU22780) ends up co-operonic.

**Dependencies:** 1.C.

**Files likely touched:** None in repo. Output goes to `/mnt/data/u/hohndor/gspa-glm/sanity/bsubtilis/`.

**Estimated scope:** XS (driver only, no code).

---

#### Task 1.E: Add `--operon-caller {heuristic,glm}` switch to integrate runner

**Description:** Wire the operon-caller choice into `benchmark/run_integrate_full_priors.sh` (or its successor `bench10_phase10.sh`). When `glm`, the script reads operons from `/mnt/data/u/hohndor/gspa-glm/preds/<tag>/operons.tsv` instead of the heuristic `${ROOT}/operons/<tag>_operons.tsv`. RefSeq→UniProt remapping flows through the same `build_refseq_uniprot_map.py` step.

**Acceptance criteria:**
- [ ] Script accepts `--operon-caller heuristic|glm` (default `heuristic`).
- [ ] Both modes produce a complete per-genome integrated TSV on the same 13 genomes.
- [ ] Per-genome integrate logs explicitly record which operon source was used.
- [ ] Heuristic mode is byte-identical to current behaviour (regression check).

**Verification:**
- [ ] Run on B. subtilis with both flags; diff the integrated TSVs to confirm meaningful change in `glm` mode and zero change in `heuristic` mode (vs. baseline).
- [ ] `grep "operon source"` in the integrate log shows the correct value.

**Dependencies:** 1.D (needs at least one real-mode operons.tsv).

**Files likely touched:**
- `benchmark/run_integrate_full_priors.sh`
- `benchmark/bench10_phase10.sh` (if used in current 13-genome run)

**Estimated scope:** S.

---

#### Task 1.F: 13-genome manifest + array sbatch run

**Description:** Build a manifest TSV (tag, fasta path, gff path, output dir) for the 13 benchmark genomes and submit `run_glm_operon.sbatch` as a SLURM array. Then run the integrate step for both heuristic and glm modes.

**Acceptance criteria:**
- [ ] `manifest.tsv` covers all 13 benchmark genomes (3 original + 10 new) with valid paths.
- [ ] Array job completes for all 13 genomes; each output directory has `operons.tsv`, `operons_confidence.tsv`, `operons_centroids.npz`, `protein_embeddings.npz`.
- [ ] Heuristic and GLM integrate runs both complete on all 13 genomes; per-genome `integrated.tsv` files written.
- [ ] No genome failed silently (verify by file count).

**Verification:**
- [ ] `sacct` shows 13/13 array tasks COMPLETED with exit 0.
- [ ] `ls /mnt/data/u/hohndor/gspa-glm/preds/*/operons.tsv | wc -l` == 13.
- [ ] Both `${OUT}/heuristic/*_integrated.tsv` and `${OUT}/glm/*_integrated.tsv` exist with > 100 lines each.

**Dependencies:** 1.E.

**Files likely touched:**
- `benchmark/cross_genome/manifest_glm.tsv` or similar (new)
- Maybe small helper script `benchmark/build_glm_manifest.sh` (new)

**Estimated scope:** S (mostly driver; rerunning sbatch is the bulk of wall-clock).

---

#### Task 1.G: Report F-max + claims-fired delta + go/no-go

**Description:** Compute micro and CAFA F-max for both runs on all 13 genomes. Tally `GenomicContextPrior` fired-claims count from each integrate log. Generate a delta table and apply the gating bar.

**Acceptance criteria:**
- [ ] Per-genome table with columns: genome, micro F-max heuristic, micro F-max glm, micro Δ, CAFA F-max heuristic, CAFA F-max glm, CAFA Δ, `GenomicContextPrior` claims (heuristic), claims (glm), claims Δ, mean operon size (h/glm), mean confidence (glm).
- [ ] Mean micro F-max delta computed across 13 genomes.
- [ ] No-regression check: count of genomes with micro F-max delta < −0.01.
- [ ] Explicit GO / NO-GO verdict against the SPEC bar (mean micro Δ ≥ +0.005 AND no genome regresses by > 0.01).
- [ ] Result written to `benchmark/RESULTS.md` as a new section "FM-operon ablation (phase 1)".

**Verification:**
- [ ] Both `print_fmax.py` runs (heuristic + glm) executed.
- [ ] Numbers reproduce on a second invocation (deterministic).
- [ ] Human review: do the deltas point in the expected direction (more claims, F-max non-negative)?

**Dependencies:** 1.F.

**Files likely touched:**
- `benchmark/RESULTS.md` (append)
- `benchmark/glm_operon_ablation.py` (new — the reporting script)

**Estimated scope:** S.

---

### CHECKPOINT 1 — Phase-1 go/no-go

**Mandatory gate before any phase-2 work.**

- [ ] `./gradlew test` green.
- [ ] All 13 genomes have both heuristic and glm integrate outputs.
- [ ] Mean micro F-max delta ≥ +0.005 across the 13 genomes.
- [ ] No genome regresses by > 0.01 micro F-max.
- [ ] User signs off on RESULTS.md addendum.

If GO → phase 2.
If NO-GO → write a one-paragraph post-mortem in RESULTS.md, drop the line, archive the branch, **do not start phase 2**.

---

### Phase 2: GENOMIC_CONTEXT_FM evidence type into Phase 7 integrator

#### Task 2.A: Build gLM contextualized centroid catalog from the 500-genome KAUST panel

**Description:** For each GO term reachable through the panel's existing GO annotations (PGAP + GOA propagation per the panel's own annotation pipeline), and for which there are ≥ N panel proteins as members, run gLM on the panel genomes' GFFs, pull the per-protein contextualized embedding for each catalog member, L2-normalize, and average across members. Mirror the schema of `build_esm2_centroids.py`.

**Acceptance criteria:**
- [ ] `build_glm_centroids.py` emits NPZ + metadata TSV with the same schema as `build_esm2_centroids.py`.
- [ ] Catalog covers ≥ 3 000 GO terms (panel is smaller than Swiss-Prot, expect lower term count).
- [ ] Catalog NPZ < 200 MB.
- [ ] Leakage check: panel members overlapping with the 13-genome benchmark or the Empty-Quarter MAGs are excluded.
- [ ] Driver records the panel manifest version + ANI-95 cluster representatives used.

**Verification:**
- [ ] Spot-check 3 well-known GO terms: their catalog members include the expected enzymes (e.g. GO:0006633 fatty-acid biosynthesis → FAS subunits).
- [ ] `numpy.load(catalog).files` matches expected schema.
- [ ] Leakage check explicit assertion in the driver script.

**Dependencies:** Phase-1 GO ✓; KAUST 500-genome panel manifest from project phase-12 (memory: `project_phase12_panel_expansion`).

**Files likely touched:**
- `benchmark/glm/phase2/build_glm_centroids.py` (new, ~200 LOC)
- `benchmark/glm/phase2/build_glm_centroids.sbatch` (new)
- `benchmark/glm/phase2/panel_manifest.tsv` (new — pinned panel reference list)

**Estimated scope:** M.

---

#### Task 2.B: Per-protein gLM-centroid predictor sidecar

**Description:** Add a `--predictor glm-centroid` mode to `benchmark/neural/run_neural_predictors.py` that reads `protein_embeddings.npz` (already produced by phase-1's `run_glm_operon.py`), and for each protein finds the top-k catalog GO terms by cosine distance. Emit standard `(protein_id, term, score, annotation_type)` TSV.

**Acceptance criteria:**
- [ ] New predictor mode plumbed through the existing manifest interface.
- [ ] Output TSV is identical schema to `esm2-centroid`.
- [ ] No GPU required for the kNN step (centroid catalog + protein embeddings both on CPU).
- [ ] Top-k clip + min-score threshold parameters available.

**Verification:**
- [ ] On B. subtilis: predicted top-1 GO terms for `argG`, `purF`, `dnaA` are sensible.
- [ ] `head` of TSV has correct columns + sane scores in (0, 1].

**Dependencies:** 2.A; phase-1 `protein_embeddings.npz` available per genome.

**Files likely touched:**
- `benchmark/neural/run_neural_predictors.py` (one new function + dispatch row — sidecar code is shared infra, lives where the existing predictors live)
- `benchmark/glm/phase2/run_panel_glmcentroid.sbatch` (new — phase-2 driver)

**Estimated scope:** S.

---

#### Task 2.C: Wire as GENOMIC_LANGUAGE_MODEL claims

**Description:** Update `benchmark/02b_parse_predictors_to_claims.py` to ingest the glm-centroid TSV with `evidence_type = GENOMIC_LANGUAGE_MODEL`. Confirm the existing `EvidenceCombiner` correlation-group machinery routes these claims correctly (`ml_genomic` group, isolated from homology and other-ML claims).

**Acceptance criteria:**
- [ ] One CLI flag (e.g. `--glm-centroid <path>`) added to the parser.
- [ ] Emitted JSONL claims have `evidence_type: GENOMIC_LANGUAGE_MODEL`, `source: glm-centroid`.
- [ ] `EvidenceTypeSpec` gains a regression test for `ml_genomic` correlation grouping.
- [ ] No existing claim type is altered.

**Verification:**
- [ ] `./gradlew :gspa-core:test --tests EvidenceTypeSpec` green.
- [ ] `python3 02b_parse_predictors_to_claims.py ... --glm-centroid <x>` emits the new evidence type and `jq` confirms.

**Dependencies:** 2.B.

**Files likely touched:**
- `benchmark/02b_parse_predictors_to_claims.py`
- `gspa-core/src/test/groovy/gspa/integration/EvidenceTypeSpec.groovy` (extend existing)

**Estimated scope:** S.

---

#### Task 2.D: Calibration + small B. subtilis end-to-end check

**Description:** Register a default Platt calibration for `glm-centroid` (default `(4, -2)` likely fine; refine after empirical evaluation). Run the full integrate pipeline on B. subtilis with phase-1 operons + GENOMIC_CONTEXT_FM claims and confirm posteriors update sensibly.

**Acceptance criteria:**
- [ ] `CalibrationTable.loadDefaults` registers `glm-centroid`.
- [ ] B. subtilis integrated.tsv includes claims with provenance `glm-centroid`.
- [ ] Spot-check: at least one (protein, GO) pair posterior moves > 0.1 due to a glm-centroid claim.

**Verification:**
- [ ] `./gradlew :gspa-core:test --tests CalibrationTableSpec` green.
- [ ] Manual: `grep glm-centroid bsubtilis_integrated.tsv | wc -l` > 100.

**Dependencies:** 2.C.

**Files likely touched:**
- `gspa-core/src/main/groovy/gspa/integration/CalibrationTable.groovy`
- `gspa-core/src/test/groovy/gspa/integration/CalibrationTableSpec.groovy`

**Estimated scope:** S.

---

#### Task 2.E: 13-genome ablation, dark-matter slice

**Description:** Re-run the 13-genome benchmark with phase-1 gLM operons baseline vs. (phase-1 ops + phase-2 GENOMIC_CONTEXT_FM claims) treatment. Compute F-max delta on (a) all proteins, (b) the dark-matter subset (proteins with no DIAMOND/Pfam hit).

**Acceptance criteria:**
- [ ] Both runs complete on 13 genomes.
- [ ] Two delta tables: all-proteins, dark-matter.
- [ ] Dark-matter subset definition recorded (no DIAMOND hit AND no Pfam hit at a stated threshold).
- [ ] No regression on all-proteins (F-max delta ≥ 0).

**Verification:**
- [ ] `python3 print_fmax.py --subset dark` produces the dark-matter slice numbers.
- [ ] Both tables added to `benchmark/RESULTS.md`.

**Dependencies:** 2.D.

**Files likely touched:**
- `benchmark/print_fmax.py` (extend with `--subset` — shared utility)
- `benchmark/glm/phase2/ablation.sh` (new — drives both arms)
- `benchmark/glm/phase2/RESULTS.md` (new — phase-2-only writeup, mirrors top-level RESULTS.md once locked in)
- `benchmark/RESULTS.md` (append summary linking to phase-2 detail)

**Estimated scope:** M.

---

#### Task 2.F: Report + go/no-go on phase 3

**Description:** Aggregate phase-2 results, draft a one-page summary (delta on all proteins, delta on dark matter, claims-fired counts).

**Acceptance criteria:**
- [ ] Summary section in `benchmark/RESULTS.md`.
- [ ] Explicit GO / NO-GO verdict for phase 3.
- [ ] User sign-off recorded.

**Bar:** GO if mean dark-matter F-max delta ≥ +0.01 AND all-protein delta ≥ 0. Otherwise stop here — phase 2 is still useful for the methods paper, but phase 3's expected value goes down.

**Dependencies:** 2.E.

**Files likely touched:** `benchmark/RESULTS.md`.

**Estimated scope:** XS.

---

### CHECKPOINT 2 — Phase-2 go/no-go

- [ ] `./gradlew test` green.
- [ ] Mean dark-matter F-max delta ≥ +0.01.
- [ ] All-protein F-max delta ≥ 0.
- [ ] User signs off on RESULTS.md addendum.

If GO → phase 3.
If NO-GO → keep phase 2 as a permanent feature, archive phase 3 plan, drop the line.

---

### Phase 3: Augment `BF(O, P)` in DarkMatterSuggester

#### Task 3.A: Build known-pathway operon corpus

**Description:** Pick a curated set of pathway-annotated reference genomes (MetaCyc / KEGG-curated, e.g. ~50 organisms with high-quality pathway maps). Run phase-1 gLM operon detection on each, group operons by participating pathway via the existing pathway DB, store each operon's centroid embedding indexed by pathway ID.

**Acceptance criteria:**
- [ ] Corpus NPZ keyed `pathway_id -> array of centroids (n_op, d_ctx)`.
- [ ] At least 100 distinct pathways with ≥ 5 operons each.
- [ ] Source genome list recorded; **no overlap with the 13-genome benchmark or EQ MAGs**.
- [ ] Build is reproducible (driver script committed).

**Verification:**
- [ ] `np.load(corpus).files | head` shows 100+ pathway IDs.
- [ ] Leakage check identical to 2.A.

**Dependencies:** Phase-2 GO ✓.

**Files likely touched:**
- `benchmark/glm/phase3/build_pathway_operon_corpus.py` (new)
- `benchmark/glm/phase3/build_pathway_operon_corpus.sbatch` (new)

**Estimated scope:** M.

---

#### Task 3.B: Plumb corpus into IntegrationState / DarkMatterSuggester

**Description:** Extend `IntegrationState` with optional `pathwayOperonCentroids: Map<String, double[][]>`. Load via a CLI flag `--pathway-operon-corpus <npz>`. Pass through to `DarkMatterSuggester`.

**Acceptance criteria:**
- [ ] CLI flag wired.
- [ ] `IntegrationState` field populated when flag set; null-safe when absent.
- [ ] `DarkMatterSuggester` exposes `pathwayOperonCentroids` field; behaviour is unchanged when null.

**Verification:**
- [ ] `./gradlew :gspa-core:test` green (existing tests must not regress).
- [ ] New test: pass empty corpus → no behaviour change.

**Dependencies:** 3.A.

**Files likely touched:**
- `gspa-core/src/main/groovy/gspa/integration/IntegrationState.groovy`
- `gspa-core/src/main/groovy/gspa/integration/suggester/DarkMatterSuggester.groovy`
- `gspa-cli/src/main/groovy/gspa/cli/IntegrateCommand.groovy` (or wherever the flag lands)

**Estimated scope:** M.

---

#### Task 3.C: Augment `computeRefinedBayesFactor` with embedding-distance term

**Description:** Add a new term to the refined BF: for each candidate operon `O` and pathway `P`, compute the mean cosine similarity from `O`'s gLM centroid to `P`'s corpus centroids. Mix into the BF as a log-odds term, weighted by a new hyperparameter `embeddingWeight`. Behind an opt-in flag (`useEmbeddingDistance`, default off) so existing tests are stable.

**Acceptance criteria:**
- [ ] New hyperparameters: `embeddingWeight` (default 1.0), `embeddingTemperature` (default 0.1), `useEmbeddingDistance` (default false).
- [ ] Behaviour identical to current code when flag is false.
- [ ] When flag is true and corpus is null → fall back to current behaviour with a single log line.
- [ ] Unit tests cover both branches + an "embeddings shift BF" assertion.

**Verification:**
- [ ] `./gradlew :gspa-core:test --tests DarkMatter*` green.
- [ ] New test verifies that an operon whose centroid is near the pathway corpus gets a higher BF than an operon whose centroid is far.

**Dependencies:** 3.B.

**Files likely touched:**
- `gspa-core/src/main/groovy/gspa/integration/suggester/DarkMatterSuggester.groovy`
- New Spock test under `gspa-core/src/test/groovy/gspa/integration/suggester/`

**Estimated scope:** M.

---

#### Task 3.D: Dark-matter ablation

**Description:** Re-run the dark-matter benchmark on the 5 pathway-rich genomes (the gapseq-equipped ones from `RESULTS.md`: ecoli, ecolo157, paeruginosa, mjannaschii, synechocystis). Compare singleton/disjunctive precision and recall with vs. without embedding distance.

**Acceptance criteria:**
- [ ] Dark-matter suggestions counted with both flag values.
- [ ] Where ground truth exists (gapseq-validated genes), precision/recall computed.
- [ ] Manual eyeball: at least 5 anecdotal cases where embedding distance corrected a circular call.

**Verification:**
- [ ] Tables in `benchmark/RESULTS.md`.
- [ ] User reviews 5 anecdotal cases.

**Dependencies:** 3.C.

**Files likely touched:**
- `benchmark/glm/phase3/dark_matter_ablation.py` (new)
- `benchmark/glm/phase3/RESULTS.md` (new — phase-3-only writeup)
- `benchmark/RESULTS.md` (append summary linking to phase-3 detail)

**Estimated scope:** M.

---

#### Task 3.E: Phase-3 closeout report

**Description:** Final summary: ablation tables, anecdotes, hyperparameter sensitivity (sweep `embeddingWeight ∈ {0.5, 1, 2}`), recommendation for default flag value.

**Acceptance criteria:**
- [ ] Closeout section in `benchmark/RESULTS.md`.
- [ ] Default flag recommendation recorded.
- [ ] User sign-off.

**Dependencies:** 3.D.

**Files likely touched:** `benchmark/RESULTS.md`.

**Estimated scope:** XS.

---

### CHECKPOINT 3 — Phase-3 closeout

- [ ] `./gradlew test` green.
- [ ] Dark-matter precision delta non-negative on the 5 gapseq genomes.
- [ ] Default flag decision recorded.
- [ ] User signs off.

---

## Risks and mitigations

| Risk | Impact | Mitigation |
|---|---|---|
| gLM API at the pinned commit doesn't expose per-pair break probabilities directly | High — blocks 1.C | Fall back to attention-pattern segmentation (paper §3.2). Also flag `--mode mock` as a viable interim for benchmark plumbing. |
| ESM2-650M doesn't fit on 1 H200 alongside gLM | Med — slows 1.C | Run ESM2 in a separate pass (write embeddings to disk) before invoking gLM. H200 has 141 GB; should be fine. |
| Phase-1 F-max delta < +0.005 | High — kills the line | Drop the line per the spec, document in RESULTS.md. The plumbing remains for future improvement. |
| Catalog leakage in 2.A: Swiss-Prot members from a benchmark genome | Med — F-max inflated | Strict accession-list exclusion before catalog build. Verify with leakage assertion in 2.A. |
| Phase-2 calibration is wrong (centroids over-fire weak terms) | Med — false positives | Default Platt `(4, -2)` is conservative; tune against held-out PGAP comparison. |
| Phase-3 corpus too small to discriminate pathways | Med — null effect | Start with 50 genomes and ≥ 100 pathways; expand if signal is on the edge. |
| Compute on `pi-hohndor` becomes contended | Low — slows wall-clock only | Use `freecycle --gres=gpu:h200:1` for non-critical reruns. |

## Resolved decisions (2026-05-05)

1. **Gating bars** (1.G, 2.E, 3.D): accepted as written for now —
   mean micro F-max Δ ≥ +0.005 with no genome regressing by > 0.01
   for phase 1; dark-matter F-max Δ ≥ +0.01 with all-protein Δ ≥ 0
   for phase 2; dark-matter precision Δ ≥ 0 on the 5 gapseq genomes
   for phase 3. Revisit if any phase produces ambiguous numbers.

2. **Catalog reference set** (2.A): the **500-genome KAUST panel**
   (project memory: phase-12 panel expansion, Strategy C + ANI-95,
   approved 2026-04-20). **Not** Swiss-Prot. The catalog will draw
   GO labels from the panel's existing PGAP annotations / GOA
   propagation. Leakage check: exclude any panel members that overlap
   the 13-genome benchmark or downstream Empty-Quarter MAGs.

3. **Branch + folder layout**: stay on `phase11-crossgenome`; new
   per-phase subfolders under `benchmark/glm/phase{2,3}/` so the
   gated phases stay separable. Phase-1 artifacts keep their
   existing locations. ORIX outputs land under
   `/mnt/data/u/hohndor/gspa-glm/phase{1,2,3}/`.

## Still open

- **Whether step 2's calibration should be source-stratified** (2.D):
  one Platt curve for `glm-centroid` globally, or per-aspect
  (BP / MF / CC)? Decide during 2.D after looking at the empirical
  raw-score distribution per aspect.
- **gLM vs gLM2 upgrade path**: out of scope for this plan; revisit
  before 3.E so the closeout's default-flag recommendation considers
  whether to also recommend the gLM2 migration.

## Parallelization notes

- 1.B and 1.C scaffolding can be developed in parallel (different files, different concerns); 1.D needs both.
- 2.A and 2.B sit in series (2.B reads the catalog), but 2.B development can start with a synthetic catalog stub so it isn't blocked.
- Phase-3 is fully sequential; no useful parallelism inside.
- Across phases — each phase is gated, so cross-phase parallelism is forbidden by the design.
