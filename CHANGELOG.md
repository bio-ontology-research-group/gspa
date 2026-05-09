# Changelog

All notable changes to GSPA are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); this project
uses [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Per-version benchmark numbers, ablation tables, and the F-max protocol
notes live in `benchmark/RESULTS.md`. This file summarises the
user-visible deltas; for measured impact, follow the cross-references.

## [Unreleased] — 1.5.x debug-surface fixes

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
  (with derived names + dominant pathway + ensemble support),
  Pathways (per-KEGG-pathway coverage with per-reaction colouring +
  cross-references to operons), Special features (AMR + BGC tables),
  Quality (GAEF metrics + named missing essentials + incoherent
  process pairs + per-pathway completeness), Pipeline. Uses the
  bundled Python templater (`gspa-cli/src/main/resources/visualize/
  make_viz.py`) extracted at runtime. No external deps beyond
  `python3` and standard library.

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
