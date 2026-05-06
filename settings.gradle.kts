rootProject.name = "gspa"

include("gspa-core")
include("gspa-predictors")
include("gspa-cli")

// Intentionally NOT included:
//   gspa-nf/   — Nextflow pipeline, runs via `nextflow run gspa-nf/main.nf`
//                See gspa-nf/README.md.
//   benchmark/ — Python evaluation harness and neural sidecar.
//                See benchmark/README.md.
