# Implementation Plan: GSPA v1.5.0 release

Companion to `SPEC.md`. Status: **Draft, awaiting approval.**

The previous plan (FM-based operon understanding) is preserved on the
`parking/phase11-glm` branch (Phase 1 NO-GO verdict, 2026-05-05).

## Overview

Ship `v1.5.0` from `phase11-crossgenome` to `main` via six phases:

- **Phase A — Triage and yank** (gating prerequisite). Diagnose
  `./gradlew test`, yank Phase 11 gLM work to `parking/phase11-glm`,
  fix CLI version drift, add GPL-3.0 LICENSE.
- **Phase B — gspa-nf integrator parity**. New `gspa-nf/modules/integrate.nf`
  with `BUILD_CLAIMS` + `INTEGRATE` processes; opt-in via `--run_integrate`.
- **Phase C — Phase 10 retune** (parallel to B). Real gapseq output, qBase
  swept over {0.50, 0.70, 0.75}, decide default-on per F-max delta.
- **Phase D — mdF comparison**. Head-to-head vs metagenomic-deepFRI
  (Bezshapkin et al. 2026) on the 13-genome PGAP panel.
- **Phase E — CHANGELOG and CI**. Keep-a-Changelog format,
  `.github/workflows/test.yml`.
- **Phase F — Merge and tag**. Open PR `phase11-crossgenome → main`,
  annotated `v1.5.0` tag.

Detailed task list lives in `tasks/todo.md`.

## Out of scope (deferred to v1.6+)

- gLM2 evaluation as alternative to gLM
- mdF wired as a first-class GSPA predictor (`MdFPredictor.groovy`)
- Phase 11.2 / Phase 11.3 (GENOMIC_CONTEXT_FM evidence + BF(O,P)
  embedding term)
- KAUST 500-genome panel scaling validation
- Eukaryote support (Augustus integration)
- Phigaro, vConTACT3, eukaryote viral predictors
- `release.yml` / `docker.yml` GitHub Actions (Phase E optional sub-tasks)
