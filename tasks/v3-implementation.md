# GSPA v3.0.0 implementation and release evidence

Started 2026-10-26 from `7714377`; development branch `development/v3.0.0`.
The authoritative scope is all twelve items in `tasks/todo.md`. The release is
authorized after implementation, biological validation, and full benchmark reporting.
No item is complete merely because its scaffolding or synthetic test passes.

## Acceptance and evidence

| Item | Required evidence before completion | Status |
| --- | --- | --- |
| 1 | Production-engine benchmark runner; frozen splits, input hashes, matched baselines and all-case reports | Pending |
| 2 | Equivalent single/batch/web/Nextflow execution; stage status; coherent gene-keyed predictions | In progress: shared prediction state |
| 3 | Correct tool contracts and end-to-end biological parser/dispatch/failure tests | In progress: isolated Claude Code Sonnet 5.5 task |
| 4 | Trained context model and held-out calibration/precision evidence; revisable inference and null assignments | Pending |
| 5 | Reaction/interaction context wired into production candidates and attribution; biological network ablations | Pending |
| 6 | Cohort recurrence with evolutionary/missing-data controls; no singleton certainty bug; independent-clade evaluation | Pending |
| 7 | Weak-candidate retention; role-specific enzyme/complex attribution; same-operon decoy evaluation | Pending |
| 8 | Explicit media/biomass/assembly assumptions and alternative absence hypotheses; phenotype validation | Pending |
| 9 | Executed community models with transport/exchanges; established coculture/auxotroph validation | Pending |
| 10 | Auditable change ledger and truthful UI/statistics reconciled to predictions | Pending |
| 11 | Qwen on unimatrix evaluated against equally informed simple models on held-out cases | Pending |
| 12 | Evidence-guided experimental planning and executed metabolic perturbations; equal-budget controls, independent phenotype evaluation | Pending |
| Release | Complete matrix above; source/tests/docs/artifacts consistent; v3.0.0 tag and GitHub release with full results | Pending |

## Benchmark scope

Use established, traceable datasets; **AllTheBacteria is excluded**. Required
tracks cover exact EC attribution (Price with accession correction), prospective
and homology-controlled GO, reference genomes across domains, operon and complex
role discrimination, taxonomic/anatomical constraints, reaction masking, metabolic
growth/knockouts, BacDive phenotypes, incomplete genomes/communities, secondary
metabolism, and agent/condition-planning comparisons. Each track must disclose
experimental versus synthetic truth and missing coverage. Dataset preparation,
model training, calibration and final test labels are separate. Benchmark
execution must use production inference, not parallel scoring formulas.

## Current work and process handles

- Root: core prediction-state repair and integration; then context/learning and
  execution parity. Local Gradle tests are permitted; large jobs use compute nodes.
- Claude Code Sonnet 5.5 worktrees `/tmp/gspa-v3-predictors`,
  `/tmp/gspa-v3-benchmarks`, and `/tmp/gspa-v3-experiments`: predictor contracts,
  production benchmark harness, and executed metabolic experiment planning.
  Prompts and streaming logs are under `/tmp/gspa-v3-agents/`.
- agy attempts ended without edits because headless auto mode was unavailable.
  Claude Code with explicit file/test permissions is running successfully.
- Unimatrix SSH verified; `/data/hohndor/gspa` exists; no hohndor Slurm jobs at
  initial inspection. Qwen endpoint is configured in the existing agent sidecar.

## Validation log

- Audit baseline: 488 JVM tests passed, 2 skipped; seven targeted Python tests
  passed. These establish the starting suite only, not v3 correctness.
- Existing untracked research data and paper build files are retained unchanged.
- Full JVM suite passed after prediction identity synchronization and recurrence
  mass-preservation fixes. Biological regressions cover competing enzyme loci,
  evidence retraction, negative recurrence, and retained alternative candidates.
  These are implementation checks, not held-out benchmark results.
