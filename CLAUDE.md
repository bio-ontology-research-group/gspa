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
- **gspa-cli** — picocli CLI: `annotate`, `evaluate`, `compare`, `report`

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

## Benchmark F-max convention

All F-max numbers in `benchmark/RESULTS.md` are **per-genome
micro-averaged F-max** (TP/FP/FN summed across all (protein, GO-term)
pairs in the genome → one F1 → max over thresholds). This is NOT
CAFA's per-protein-then-averaged F-max. See the "F-max definition"
section of `benchmark/RESULTS.md` for the full procedure and rationale.
Implementation in `benchmark/benchmark_pgap_v2.py::fmax_with_ci()`.
