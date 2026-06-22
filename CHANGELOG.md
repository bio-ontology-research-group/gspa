# Changelog

All notable changes to GSPA are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); this project
uses [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Per-version benchmark numbers, ablation tables, and the F-max protocol
notes live in `benchmark/RESULTS.md`. This file summarises the
user-visible deltas; for measured impact, follow the cross-references.

## [Unreleased]

### Added
- **DeepGO-PlusPlus full ablation + aggregator study** (`deepgo-plusplus/pipeline/ablation.py`;
  raw numbers `deepgo-plusplus/ablation_no_results.tsv`, deep dive in
  `RESULTS.md`). One reproducible pass (loads each component once, scores all
  configs in a single `cafaeval` call) reports per-component standalone,
  leave-one-out, and cumulative-build-up f_w plus a **logistic-regression vs
  gradient-boosted-trees** aggregator comparison. Findings: **`net` is the only
  large marginal contributor** (+0.049 LOO); the PLM/homology/domain channels are
  largely redundant for novel proteins and **`clean` *decreases* the score**;
  **leak-free XGBoost does not beat the linear stacker** (0.478 vs 0.483) and
  XGBoost + term-identity features **leaks** to a non-credible MF 0.963 — the
  bottleneck is signal diversity/generalisation, not aggregator capacity. This is
  why DG++ ships the linear, scores-only integrator.
- **DeepGO-PlusPlus-Light — no-GPU variant** (`eval_light.py` + `build_cnn_component.py`
  + `build_net_bridge.py` + `extract_sprot_fasta.py`; models
  `deepgo_plusplus_light_cpu.json`, `deepgo_plusplus_light.json`,
  `deepgo_plusplus_light_cnn.json`). Drops the GPU PLM heads (`mlp`/`prostt5`/
  `esm2_3b`) and uses only CPU components (DIAMOND, FoldSeek over AlphaFold-DB
  structures, InterProScan/-LR, STRING Net-KNN, optional BM25 `lit`) plus a CPU
  **1D-CNN over sequence** (`cnn`, DeepGOCNN-style, frequency-weighted BCE, trained
  on the 80,750 pre-t0 SwissProt proteins) meant to replace the PLM heads. Runs
  through the existing predictor unchanged (component list read from the JSON).
  Two CAFA6 findings (no-knowledge IA-weighted f_w; tables in `RESULTS.md`):
  **(1) the best no-GPU panel `diam+foldseek+interpro+net` (0.550) beats the full
  GPU model (0.532) on novel proteins** — the PLM heads are redundant with `net`;
  **(2) the 1D-CNN does not improve novel-protein f_w** (standalone 0.206; −0.012
  in panel) — a sequence model trained on pre-t0 data hits the same generalisation
  wall, so it ships as the separate coverage-first `_cnn` model (0.516, predicts for
  every protein incl. orphans) rather than in the default.
  **(3) `foldseek` is not unconditionally GPU-free** (it needs a query structure —
  a CPU AFDB lookup, but folding a novel sequence needs a GPU); the strictly-no-GPU
  panel drops it (`diam+interpro+net`, 0.544, −0.006). **(4) `net` only fires for
  STRING members**, so `build_net_bridge.py` adds a DIAMOND homology bridge (query →
  pre-t0 STRING-member homolog → neighbour labels) that takes the 356 no-STRING
  no-knowledge proteins from f_w 0 → 0.42 and lifts the structure-free panel to
  **0.564** — shipped as **`deepgo_plusplus_light_cpu.json`** (`diam,interpro,net_union`),
  the **recommended strictly-no-GPU model** (works for any sequence, incl. proteins
  not in STRING). Bridge is leak-safe (pre-t0 homolog DB; novel queries can't
  self-match).
- **`deepgo-plusplus` module + reproducible retraining pipeline**
  (`deepgo-plusplus/`). The learned-stacker predictor is consolidated into a
  self-contained module — the predictor was renamed **`cafa-baseline` →
  `deepgo-plusplus`** (the old id stays a working alias everywhere: sidecar
  runner, CLI flags `--cafa-baseline*`, and the `RUNNERS` registry). The new
  folder is **re-runnable at every UniProt / STRING release**: a `Makefile`
  drives the full DAG (UniProt index → CAFA6 ground truth → Net-KNN/literature
  components → train + freeze the integrator → `cafaeval`), inputs are declared
  in `config.mk`, dependencies are pinned in `pyproject.toml` (uv), and input
  release provenance is recorded in `VERSIONS.md`. Apply via the same canonical
  runner (`pipeline/apply_integrator.py` reuses `run_deepgo_plusplus` — one
  copy of the integration math). Frozen models shipped under
  `deepgo-plusplus/models/` (`deepgo_plusplus_integrator{,_net,_lit_net}.json`).
  **Comprehensive regression tests** (`deepgo-plusplus/tests/`, pure-CPU/offline
  pytest): ground-truth knowledge-class + t0 logic, the UniProt index parser,
  Net-KNN (vote/min-conf/corrupt-file skip), the literature **name-only leak
  guard** + sharding, integrator schema/stability/leak-config guard, the frozen
  apply (DAG propagation, min-score, the stable-sigmoid extreme-`z` regression),
  an end-to-end signal-recovery chain, and a shipped-model schema check.

- **DeepGO-PlusPlus learned-stacker GO predictor** (`DeepGoPlusPlusPredictor`,
  sidecar id `deepgo-plusplus`). A CAFA6-competitive predictor built with no new
  model architecture: it replaces naive max-merge of GSPA's component predictors
  (DIAMOND/BLAST-KNN, FoldSeek-KNN, CLEAN, InterPro, an ESM2 MLP head and a
  structure-aware ProstT5 head) with a **frozen per-aspect logistic-regression
  stacker**. On a faithful CAFA6 reconstruction (GOA snapshot,
  t0 = 2026-02-02, official `cafaeval` + `IA.tsv`) this recovers
  **novel-protein IA-weighted f_w from 0.359 → 0.483** (vs the 0.524
  first-place GOAlpha entry; our original entry was 0.377, rank 263/2177) —
  the gain is integration done right, not bigger models. ProstT5 is the one
  PLM that complements rather than cannibalises the MLP head; ESM2-3B and
  LR-InterPro proved redundant.

  At inference there is no ground truth, so the integrator is trained once and
  **frozen** to a 2 KB JSON
  (`deepgo-plusplus/pipeline/train_integrator.py --save-model`); an apply-only
  sidecar runner (`run_neural_predictors.py --predictor deepgo-plusplus`)
  combines precomputed per-component score TSVs, propagates each to GO
  ancestors (max) and emits `sigmoid(w·x + b)`. Applying the frozen model to
  the no-knowledge set scores 0.489, confirming the deployable artifact is
  faithful. New `DeepGoPlusPlusConfig` block in `GspaConfig`; wired into
  `AnnotationPipeline.createAllPredictors()`; `--deepgo-plusplus*` CLI flags
  (legacy `--cafa-baseline*` aliases retained); Spock coverage for flag
  serialisation, 4-column output parsing and fail-fast on missing config.
  Shipped model `deepgo-plusplus/models/deepgo_plusplus_integrator.json`;
  recipe and verified numbers in `deepgo-plusplus/README.md` + `RESULTS.md`.
  Off by default.

- **DeepGO-PlusPlus network + literature signals (Net-KNN, BM25 text-kNN).**
  Two more components for the stacker, reproducing the rest of GOAlpha's
  heterogeneous panel. **Net-KNN** (`net`) votes a query protein's STRING-v12
  PPI neighbours' pre-t0 GO labels (guilt-by-association); it is leak-free
  (STRING 2023 edges + `train_terms` labels both pre-t0) and the clear win —
  adding it lifts the official 3-class IA-weighted f_w on the CAFA6
  reconstruction **0.629 → 0.647** (no-knowledge 0.489 → 0.538), with gains
  across novel/limited proteins and only a small partial-knowledge dip.
  **Literature** (`lit`, BM25 text-kNN over SwissProt names) pushes
  no-knowledge higher (→ 0.553 combined) but lowers the 3-class mean by hurting
  partial-knowledge proteins and carries a name-leakage caveat, so it ships as
  an optional no-knowledge booster, not a default. New shipped models
  `deepgo_plusplus_integrator_net.json` (recommended) and
  `deepgo_plusplus_integrator_lit_net.json`; builders
  `pipeline/build_text_string_index.py`, `pipeline/build_net_component.py`,
  `pipeline/build_lit_component.py`, `benchmark/neural/run_net_ws.sh`. The
  Groovy predictor is unchanged (component list lives in the integrator JSON).

### Changed
- **Renamed `cafa-baseline` → `deepgo-plusplus`** across the predictor, config
  (`GspaConfig.neural.deepGoPlusPlus`), CLI flags, sidecar runner and shipped
  model filenames, and relocated the training scripts + models + docs into the
  new `deepgo-plusplus/` module. The `cafa-baseline` id and `--cafa-baseline*`
  flags continue to work as aliases, so existing pipelines and frozen models are
  unaffected.

### Fixed
- **DeepGO-PlusPlus sidecar: numerically stable logistic sigmoid.** The apply
  path computed `1/(1+exp(-z))` directly, which raised `OverflowError` for
  extreme component combinations (large negative `z`); switched to the
  branchless-by-sign stable form. The 6-component model was unaffected; the
  fix is required for the larger net/lit models. Now pinned by a regression
  test (`tests/test_apply_integrator.py::test_extreme_z_does_not_overflow`).

## [1.5.3] — 2026-05-11 — AntiFam pseudogene scanner

### Added
- **AntiFam pseudogene scanner** (`AntiFamPredictor`). Runs
  `hmmsearch --cut_ga` against the AntiFam HMM library (EBI/Pfam),
  which is the same pseudogene / spurious-ORF filter that Bakta uses.
  Hits become `Annotation` records with `type = PSEUDOGENE`,
  `source = 'antifam'`, and the AntiFam accession as the value, so
  downstream consumers (integrator, quality scorer, visualiser) can
  drop or down-weight claims attached to flagged proteins.
  Off by default — download the database from
  `https://ftp.ebi.ac.uk/pub/databases/Pfam/AntiFam/current/AntiFam.tar.gz`
  and enable via:
  ```yaml
  predictors:
    antifam:
      enabled: true
      database: /refs/AntiFam.hmm
  ```
  New `AnnotationType.PSEUDOGENE` enum value; new `AntiFamConfig`
  block in `GspaConfig`; wired into `AnnotationPipeline.createAllPredictors()`.
  Spock spec covers domtblout parsing, duplicate-domain
  de-duplication, the `--cut_ga` / `--domE` toggle, and the default
  output type. Validated end-to-end through Singularity on unimatrix01
  against the MR59-6 *Pontibacter* prokka FAA (4,372 proteins;
  0 canonical hits, 11 sub-threshold AntiFam matches under loose
  E-value gating). Bakta-as-replacement-gene-caller (parsing Bakta's
  GFF `pseudo=true` attributes) is a separate v1.6 task; this release
  adds the AntiFam capability that Bakta's pseudogene step depends on.
- Docker image `leechuck/gspa-cli:1.5.3` published; bundles the new
  predictor. `gspa-nf/modules/integrate.nf` pinned to the new tag.

## [1.5.2] — 2026-05-11 — container + documentation hygiene

Closes the documentation / containerisation gaps surfaced after 1.5.1
shipped (none of the new 1.5.x features were end-to-end-validated under
Docker, and the top-level README didn't mention any of them).

### Added
- **`leechuck/gspa-cli:1.5.2`** Docker image, built from
  `docker/Dockerfile` and published on Docker Hub. java 21 + python3 +
  DIAMOND + HMMER + BLAST+ + prodigal + barrnap + pyrodigal + cobra +
  embedded `gapsmith` (Rust). Sized for orchestrating `gspa annotate /
  integrate / visualize` end-to-end without depending on the heavier
  predictor-stack images. Replaces the unpublished `gspa/gspa:latest`
  placeholder previously referenced by `gspa-nf/modules/integrate.nf`.
- README v1.5.x section walks through `gspa visualize`,
  `benchmark/fetch_kegg_modules.py`, and pathway-source stacking via
  `--modules`. `gspa-nf/README.md` documents the six integrate-path
  processes (`BUILD_CLAIMS`, `SIDECAR_CLAIMS`, `MERGE_CLAIMS`,
  `OPERONS`, `INTEGRATE`, `VISUALIZE`) and the orchestration image.
- CI `docker-build` job (`.github/workflows/test.yml`): runs after
  `gradle-test`, builds `docker/Dockerfile` end-to-end, smoke-tests
  `gspa --version` and `gspa --help | grep visualize` inside the
  resulting image. Future Dockerfile drift now fails CI rather than
  being caught only at release time.

### Fixed
- `docker/Dockerfile` — `minced` package no longer in Ubuntu 24.04
  (Noble) repos; image was failing to build at the apt-install step.
  Dropped from the orchestration image; CRISPR detection remains
  available through the dedicated `quay.io/biocontainers/minced`
  container that `gspa-nf` already wires for the `MINCED` process.
- `gspa-nf/modules/integrate.nf` — `VISUALIZE` referenced
  `container 'gspa/gspa:latest'` which had never been published, so
  `--run_visualize` would have failed at the image-pull step. Now
  pinned to `leechuck/gspa-cli:1.5.2`. `MERGE_CLAIMS` was missing an
  explicit `container` directive (relied on the default); pinned to
  `python:3.12-slim` for parity with the other Python-only steps.



Discovered while building an interactive HTML browser for the MR59-6
*Pontibacter* tutorial: the visualisation surfaced a real product bug
in the integrator (only 5 of 10 tools contributing to posteriors)
plus several reporting gaps. All fixed end-to-end. The bug fix alone
shifts MR59-6 GAF coverage from 53.8% → 94.0% (15,938 → 75,130 GAF
rows; integrated annotations 84,215 → 311,440).

### Fixed
- `ClaimExtractor.SOURCE_TO_TYPE` was missing entries for `mdf`,
  `mdeepfri`, `proteinfer`, `clean`. Claims with these source names
  silently dropped at extract time (line 158: `if (type == null)
  return // unresolved claim; skip`), so 261k claims from those 3 tools
  never reached the Bayesian integrator. Added entries mapping all four
  to `EvidenceType.SEQUENCE_DEEPLEARNING`. Regression test in
  `ClaimExtractorSpec` loads a JSONL with each new source and asserts
  the claims survive.

### Added
- **GAEF report detail** — `quality_gspa.json` now ships per-pathway
  and per-process detail with human-readable names (not just bare GO
  ids), so callers can ask which essentials are missing or which
  pathway is incoherent without their own GO-ontology lookup.
  - `QualityReport.incoherentProcessPairs` — list of unsatisfied
    `(required, missing)` `has_part` pairs from the process-coherence
    check.
  - `QualityReport.incoherentPathways` — per-triggered-pathway
    completeness with present/missing GO terms.
  - `QualityReport.goLabels` — id → name lookup populated by
    QualityScorer / QualityPipeline whenever a `GoOntology` is
    available; consumed by `QualityReportWriter` so the JSON has both
    ids and names.
  - `Coherence.evaluate()` now wires the per-pair / per-pathway detail
    through `CoherenceResult` to the report.
  - `PathwayDatabase.computePathwayCoherenceDetailed()` returns the
    per-pathway result object alongside the existing aggregate score
    (backwards compatible).
  - `QualityReportWriter.buildReportMap` emits two new sections:
    `coherence.process_unsatisfied_pairs` and `coherence.pathway_detail`,
    plus `*_named` variants of the essential-functions lists.

- **Operon predictor ensemble + Noisy-OR** — `OperonEnsemble` runs three
  independent predictors per adjacent gene-pair (distance ≤ 300 bp,
  strict ≤ 50 bp, functional with shared GO BP terms ≤ 1000 bp) and
  combines per-predictor sensitivity θ via Noisy-OR. `Operon` now
  carries `supportSet`, `minPairPosterior`, `meanPairPosterior` so
  callers can rank operons by confidence. The original
  `OperonPredictor` (single rule) is unchanged.

- **`gspa annotate` persists operons by default** — `AnnotateCommand`
  re-detects operons after the predictor stack runs and writes
  `operons.tsv` (verbose), `protein_to_operon.tsv` (reverse index),
  and `operons_for_integrate.tsv` (the format `gspa integrate
  --operons` consumes) to the output dir. Cheap O(n log n) re-run
  that closes the gap between "operons were predicted" and "downstream
  steps can find them on disk".

- **`gspa visualize` subcommand** — emits a single self-contained HTML
  browser (~25 MB with embedded FASTA) for any GSPA workspace. Tabs:
  Proteins (virtualised search/filter table), Functions (confidence
  histogram + aspect donut + top GO terms), Genome browser (igv.js
  with CDS / operons / BGCs / AMR / localisation tracks), Operons
  (enrichment-based names + ensemble support + dominant pathway when
  significant), Pathways (per-KEGG-pathway coverage with per-reaction
  colouring + cross-references to operons), Special features (AMR +
  BGC tables), Quality (GAEF metrics + named missing essentials +
  incoherent process pairs + per-pathway completeness), Pipeline.
  Uses the bundled Python templater (`gspa-cli/src/main/resources/
  visualize/make_viz.py`) extracted at runtime. No external deps
  beyond `python3` and standard library.

- **Operon naming + pathway tags via hypergeometric enrichment** — the
  visualisation derives operon names by running a one-sided
  hypergeometric test on every GO BP term carried by ≥2 operon
  members (genome-wide BP-annotation frequency as the background).
  The smallest p-value with k≥2 wins; ties broken by global rarity
  (more specific first); name shown with `(k/M, p, fold)` qualifier.
  Falls back to the dominant non-hypothetical Prokka product when no
  term reaches p<0.05. On MR59-6 this recovers textbook bacterial
  operons (trp, F-ATPase, NUO, ribosomal-protein cluster, thiamine
  biosynthesis, …) at typical fold-enrichment 30–500×.

  Pathway tags use the same hypergeometric over pathway-membership
  in the genome background. Top-3 enriched pathways shown per operon;
  a pathway is labelled "dominant" only when k≥2 AND coverage≥25%
  AND p<0.05 — kills the noisy "1/M" tags that were misleading
  (single-member pathway hits are statistically meaningless on a
  large bacterial genome). Display caps coverage at 100% when the
  same enzyme slot is hit by paralogous operon members
  (e.g. NUO complex subunits all carrying GO:0008137).

- **gspa-nf integrator pipeline upgrades** —
  - `SIDECAR_CLAIMS` process: parses mDeepFRI / ProteInfer / CLEAN
    sidecar TSVs into a parallel claims.jsonl. With this + the
    ClaimExtractor fix, the Nextflow path now matches the JVM CLI's
    8-source claim coverage.
  - `MERGE_CLAIMS` process: concatenates builtin + sidecar claims.
  - `OPERONS` process: wraps the bundled 3-predictor operon ensemble
    so the Nextflow path produces the same `operons_for_integrate.tsv`
    that `gspa annotate` writes locally.
  - `INTEGRATE` process now consumes operons (not just claims) so
    `GenomicContextPrior` fires.
  - `VISUALIZE` process scaffold (gated on `params.run_visualize`,
    auto-trigger pending until `gspa evaluate` is wired in).
  - All processes documented inline in `gspa-nf/modules/integrate.nf`.

### Changed
- `Operon` model gains three optional fields (`supportSet`,
  `minPairPosterior`, `meanPairPosterior`); existing callers that
  construct `Operon` without them keep working (defaults are empty).

- **Pathway database supports stacking multiple sources**. New
  `PathwayLoader.loadPathwaysInto(db, file)` merges an additional
  pathway TSV into an existing `PathwayDatabase` so KEGG main maps,
  KEGG Modules, MetaCyc, BioCyc — and any future source in the same
  schema — can be layered. `gspa integrate --modules` and
  `gspa evaluate --modules` accept comma-separated TSV paths to stack.

- **KEGG Modules ingest** — new `benchmark/fetch_kegg_modules.py`
  fetches KEGG Modules from the public REST API (`/list/module` +
  `/link/ec/module`) and joins to GO via `ec2go.txt` to emit
  `kegg_modules.tsv` in the existing `kegg_pathways.tsv` schema. 520
  modules / 2,822 (module, EC) rows. Modules are smaller, more focused
  units (5–15 enzymes) than KEGG main pathways (50+), so per-genome
  enrichment is meaningful — on MR59-6, switching to main + Modules
  raised the count of operons with a significant dominant pathway
  from **20 → 58** while keeping the strict bar (k≥2, coverage≥25%,
  p<0.05) unchanged.

  The Pathways tab in the visualisation now classifies each entry by
  source (KEGG main / KEGG Module / other) with a per-source filter
  and badge.

### Fixed
- `PathwayGraph.computeCompleteness` used path-based scoring that
  collapsed to 1.0 whenever a pathway had no dependency edges (every
  single-node "path" reported 0 or 1; max → 1 if any enzyme present).
  KEGG Modules built from `link/ec/module` carry no reaction order,
  so all 274 module entries reported 100% coverage on MR59-6 — the
  metric was degenerate. Fixed: when the reaction graph has no edges,
  fall back to fraction-of-required-terms-covered, the honest
  aggregate. MR59-6 pathway coherence is now 95.7% → 55.8% with
  modules in the mix (314 of 397 triggered pathways are genuinely
  partial), and the per-pathway detail in the report finally calls
  out which enzymes are missing.

- `GoOntology.getLabel` and `buildLabelCache` rendered every GO label
  as the literal string `"true"`. Root cause: in Groovy, the property
  shorthand `((OWLLiteral) ann.value).literal` resolves to
  `OWLAnnotationValue.isLiteral()` (returns the boolean `true`), not
  `OWLLiteral.getLiteral()` (the lexical form). Coerced to `String`
  on the way into the label cache, every entry became `"true"`. The
  `goLabels` map shipped in `quality_gspa.json` therefore mapped every
  id to `"true"`, so the visualisation rendered "missingtrue
  GO:0090482" instead of "missing vitamin transmembrane transporter
  activity GO:0090482" for every essential function, every
  process-coherence pair, and every pathway term. Fix: explicit
  `getLiteral()` calls; defence-in-depth `resolveGoName` helper in the
  HTML viewer falls back to the top-level GO map when an inline name
  is `"true"`/`"false"` so already-generated quality JSONs render
  correctly without re-running `gspa evaluate`. Regression test in
  `GoOntologySpec` loads the bundled `go-tiny.obo` and asserts
  `getLabel('GO:0006259') == 'DNA metabolic process'`.

## [1.5.0] — 2026-05-07

### Added
- gspa-nf integrator parity: `--run_integrate` flag wires
  `BUILD_CLAIMS` + `INTEGRATE` processes after `MERGE_ANNOTATIONS`,
  producing per-(protein, function) posterior probabilities with the
  full Phase 7 prior stack from inside Nextflow. Closes the documented
  parity gap with the JVM CLI.
  (`gspa-nf/modules/integrate.nf`, `gspa-nf/README.md` "End-to-end
  with Phase 7 integration".)
- GPL-3.0-or-later `LICENSE` at repo root with dependency
  compatibility notes in `README.md` (covers OWL API LGPL, ELK / picocli
  / Jackson / Spock Apache-2.0, SAT4J via the LGPL leg, AGPL-3.0
  DeepEC via subprocess sidecar boundary).
- `VersionProvider` reads the GSPA version from a build-generated
  `version.properties` resource. Picocli's `--version` flag, the
  shadowJar artefact name, and the Gradle project version are now
  guaranteed in sync; future bumps require editing only the root
  `build.gradle.kts`.
- `CHANGELOG.md` (this file) and `.github/workflows/test.yml`.
  CI runs `./gradlew clean test` + `:gspa-cli:shadowJar` + version
  smoke on push to `main` and on PRs targeting `main`.
- Comparison vs metagenomic-deepFRI (mdF; Bezshapkin et al., bioRxiv
  2026-04-29). Sequence-only mdF F-max micro = 0.157, CAFA = 0.153
  on the 13-genome PGAP-comparison panel; GSPA C1 baseline averages
  0.842 / 0.868. **Mean GSPA / mdF ratio = 5.4× (micro), 5.6× (CAFA).**
  Frames the integrator's value vs. single-modality structure-aware
  prediction. Full tables in `benchmark/RESULTS.md` "v1.5.0 —
  comparison with metagenomic-deepFRI".
- mdF→GSPA-shape adapter `benchmark/parse_mdf_predictions.py`
  (BSD-3-licensed mdF v1.0 weights consumed via the canonical
  `results.tsv` schema; --self-test passes).
- Phase C / Phase D SLURM array drivers + scoring scripts under
  `benchmark/`: `phase10_retune.{sh,sbatch}`,
  `phase10_retune_array.{sh,sbatch}`, `mdf_array.{sh,sbatch}`,
  `score_phase_c.sh`, `score_mdf.sh`.

### Changed
- Phase 10 outer iterative loop (`--iterate-gapseq`) retune verdict
  for v1.5.0: **NO-GO for default-on**. Higher `qBase` (0.50 → 0.70 /
  0.75) narrows the regression vs. the C1 baseline (Δ −0.029 → −0.027
  micro F-max) but every tuned variant still regresses on every
  genome. q=0.70 and q=0.75 collapse to the same numbers because the
  default `qCap = 0.75` saturates both threshold paths. Default-off
  retained as in v1.4.x; flag remains opt-in for users who want it.
  Full data in `benchmark/RESULTS.md` "v1.5.0 — Phase 10 retune".

### Removed
- Phase 11 gLM operon caller and supporting benchmark artefacts.
  Gating verdict was NO-GO (mean micro F-max Δ = +0.0002, fails the
  +0.005 threshold). Full Phase 11 history preserved on the
  `parking/phase11-glm` branch for any future reroll (e.g. with
  gLM2). The `--operon-caller {heuristic,glm}` switch in
  `benchmark/run_integrate_full_priors.sh` is also gone.

## [1.4.1] — 2026-04-26

### Added
- ScanNet end-to-end validated on unimatrix01 (PPI interface sites
  via Apache-2 SIF + cwd-symlink workaround for /ScanNet read-only
  layout). Closes the v1.4 FOSS-protein-predictor line: 10/10 live.

See `benchmark/RESULTS.md` "v1.4 — Track A finishers" for the
hpylori-panel install gotchas and per-tool row counts.

## [1.4.0] — 2026-04-21

### Added
- TPpred3 (GPL-3.0, transit peptides), MusiteDeep (MIT, PTM sites),
  TMbed (Apache-2.0, TM regions via ProtT5) via per-tool Singularity
  images.
- Track B genomic-region ensemble: `build_genomic_ensemble.py` fuses
  geNomad / CheckV / PhiSpy / VirSorter2 / VIBRANT 6-column genomic
  region TSVs by reciprocal interval overlap.
- VirSorter2 + VIBRANT wired into main.nf as `params.run_virsorter2`
  / `params.run_vibrant` (FOSS biocontainers).

See `benchmark/RESULTS.md` "v1.4 — Track A finishers + viral
expansion + genomic-region ensemble".

## [1.3.0] — 2026-04-15

### Added
- Track A — DeepFRI (BSD-3, sequence-only GO via GCN) and DeepEC
  (AGPL-3, EC numbers) productionised after upstream API and pickle
  shimming. v1.3 hpylori panel: DeepFRI 10,591 GO rows; DeepEC 366 EC.
- Track B — phage / prophage genomic-region track with geNomad
  (Apache-2), CheckV (BSD-3), PhiSpy (BSD-3) using a 6-column
  genomic-region TSV shape. RDF / JSON-LD vocabulary extended with
  `gspa:Prophage`, `gspa:ViralContig`, etc.

See `benchmark/RESULTS.md` "v1.3 — Track A predictor fixes + Track B
phage / prophage track".

## [1.2.0] — 2026-04-08

### Added
- 10 OSI-licensed FOSS protein predictors with three new output
  shapes (region, term-extras, site): Metapredict v2, DeepSig v3,
  TMbed, TPpred 3, PSORTb 3.0, DeepFRI, DeepEC, DeepARG, MusiteDeep,
  ScanNet. Replaces license-walled DeepLoc 2 / DeepTMHMM / IUPred3 /
  TargetP 2 / SignalP 6 / NetPhos for production use.
- HTML, RDF/Turtle, JSON-LD report extended to all three shapes; SIO
  vocabulary added for `gspa:Region`, `gspa:Site`, `gspa:Prophage`
  and friends.
- Three new Docker images: `leechuck/gspa-region-stack:0.1`,
  `leechuck/gspa-tf-stack:0.1`, `leechuck/gspa-struct-stack:0.1`.

See `benchmark/RESULTS.md` "v1.2 — FOSS-only fast ML predictors".

## [1.1.0] — 2026-04-04

### Added
- Modern ML/DL function predictors via the neural sidecar pattern:
  ESM2-DeepGOPlus (frozen ESM2 + FC head), ProteInfer (CNN),
  ESM2-centroid (NPZ centroids over Swiss-Prot), CLEAN (ESM2 +
  contrastive head, EC).
- Multi-format reports (HTML / TTL / JSON-LD) replacing the v1.0
  text-only report.
- Phase 12 cross-genome layer scaffolding: `ReactionLocalContextSuggester`
  (M1, in-tree), `ReactionLocusCatalog` + `CrossGenomeReScorer` (M2,
  conditional-LR re-weighting with Jeffreys-prior CIs), `gspa.integration.ranker`
  package (M3 scaffold).

## [1.0.0] — 2026-03-25

### Added
- Initial release. Multi-module Gradle build (Java 21+, Groovy 4),
  five `gspa-cli` subcommands (annotate, evaluate, compare, report,
  integrate), Phase 7 evidence integrator (Noisy-OR with correlation
  groups, IterativeRefiner, 5 priors: Essentiality, Coherence,
  Consistency [SAT4J UNSAT-core], GapFilling, GenomicContext), GAEF
  quality metrics, Phase 8 one-shot DarkMatterSuggester. 22 predictor
  wrappers covering similarity (DIAMOND, MMseqs2), domains (Pfam,
  InterProScan), structure (FoldSeek), orthology (eggNOG-mapper),
  pathways (gapseq), localization (SignalP, DeepTMHMM), and several
  specialized tools (AMRFinder, dbCAN, antiSMASH, VFDB).
- gspa-nf Nextflow sibling (Docker / Singularity / SLURM profiles).
- 13-genome head-to-head benchmark vs PGAP, mean GSPA/PGAP F-max
  ratio = 1.93× micro, 1.96× CAFA. See `benchmark/RESULTS.md`
  "Main result" and "10-Genome PGAP Comparison".

[1.5.0]: https://github.com/bio-ontology-research-group/gspa/compare/v1.4.1...v1.5.0
[1.4.1]: https://github.com/bio-ontology-research-group/gspa/compare/v1.4.0...v1.4.1
[1.4.0]: https://github.com/bio-ontology-research-group/gspa/compare/v1.3.0...v1.4.0
[1.3.0]: https://github.com/bio-ontology-research-group/gspa/compare/v1.1.0...v1.3.0
[1.1.0]: https://github.com/bio-ontology-research-group/gspa/compare/v1.0.0...v1.1.0
[1.0.0]: https://github.com/bio-ontology-research-group/gspa/releases/tag/v1.0.0
