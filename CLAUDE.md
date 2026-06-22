# GSPA — Genome-Scale Protein Annotation

## Build & Test

```bash
./gradlew build          # compile + test all modules
./gradlew clean test     # clean build with all tests
./gradlew :gspa-cli:run --args='--help'   # run CLI
```

Java 21+, Gradle 8.7 (wrapper included).

## Project Layout

Multi-module Gradle (Groovy 4 + Java):

- **gspa-core** — Data model, ontology (OWL API + ELK), quality metrics (SAT4J for consistency), I/O, config
- **gspa-predictors** — Predictor interfaces + tool wrappers (DIAMOND, InterProScan, FoldSeek, eggNOG, gapseq, operons, etc.), gene callers, crossfeeding analyzer, AnnotationPipeline orchestrator
- **gspa-cli** — picocli CLI: `annotate`, `evaluate`, `compare`, `report`. Each subcommand lives in its own file (`AnnotateCommand.groovy`, `EvaluateCommand.groovy`, `CompareCommand.groovy`, `ReportCommand.groovy`); `GspaMain.groovy` only wires them together.

Also in the repo but NOT part of the Gradle build:

- **gspa-nf/** — Nextflow sibling pipeline (`main.nf` + 5 `modules/*.nf`), runs the same external tools as the JVM CLI but via Singularity / Docker on HPC. Run with `nextflow run gspa-nf/main.nf`. Convenience Gradle tasks: `./gradlew nfHelp` (no-Nextflow quickstart) and `./gradlew nfLint` (parse check, skipped when Nextflow isn't on PATH). See `gspa-nf/README.md` for the full story and `gspa-nf/UNIMATRIX01.md` for the worked cluster example.
- **benchmark/** — Python evaluation harness + neural-predictor sidecar (`benchmark/neural/run_neural_predictors.py`). See `benchmark/README.md` for script layout.
- **deepgo-plusplus/** — Self-contained, reproducible retraining pipeline for the `deepgo-plusplus` predictor (the learned-stacker CAFA6 baseline; legacy sidecar/CLI alias `cafa-baseline`). `Makefile` rebuilds the frozen integrator from a UniProt/STRING/CAFA release; deps pinned via `uv` (`pyproject.toml`); regression suite in `deepgo-plusplus/tests/` (`uv run make test`). Inference still runs through the shared neural sidecar. See `deepgo-plusplus/README.md`.

The exclusions are recorded in `settings.gradle.kts` so they aren't mistaken for "forgotten" subprojects.

## Key Design Decisions

- **Consistency checking uses SAT4J** (not ELK) — taxon constraints encoded as propositional SAT with sibling disjointness. UNSAT core for violation explanation.
- **ELK** used only for completeness (subsumption) and coherence (has_part pair extraction).
- Essential function profiles are **runtime-configurable** (add/remove GO terms per run).
- All predictors implement `Predictor` interface; genome-level ones implement `GenomePredictor`.
- External tools wrapped via `AbstractToolPredictor` (command + parse pattern).

## Tests

Spock framework. Tests use synthetic data in `gspa-core/src/test/resources/test-genomes/` and `test-ontology/`.

## Conventions

- Groovy source, Spock tests, Gradle Kotlin DSL for build files
- Package root: `gspa.*`
- YAML config with hierarchical merging (defaults → kingdom preset → user file → CLI)

## Benchmark F-max conventions

`benchmark/benchmark_pgap_v2.py` reports **two F-max metrics per
genome**:

1. **Micro** (`fmax_with_ci`) — per-genome TP/FP/FN summed across all
   (protein, GO-term) pairs → one F1 → max over thresholds.
2. **CAFA** (`fmax_cafa_with_ci`) — CAFA III/IV protein-centric:
   per-protein precision/recall averaged across proteins → F1 from
   those averages → max over thresholds.

The result tables show both, side by side. See the "F-max definitions"
section of `benchmark/RESULTS.md` for the full procedures and rationale.
