"""Six follow-up experiments on Stage 2B and off-growth enzyme discovery (design: tasks/v3-stage2b-followups-design-2026-10-31.md).

Modules:
  common          shared loaders (layouts, KEGG modules, stored off-growth outputs, labels, receipts)
  assignment      competitive gene assignment rules (experiment 1) and the shared responsibility rule
  exp1_assignment re-scoring of the stored off-growth outputs under the assignment rules (experiment 1)
  exp3_panel      per-genome condition panels from transporter calls, re-scored on stored sole-source tests (experiment 3)
  exp2_targets    unblocking / dead-end tests as a gated evidence source for the Stage 2B update (experiment 2)
  stage2b_rerun   re-attribution of the stored Stage 2B cells with covariates and candidate limits (experiments 4, 5)
  price_score     Price scoring of sparse update records against stage 0 with the frozen metric functions
  exp6_fitness    knockout predictions on the fixed panel against Fitness Browser specific phenotypes (experiment 6)
"""
