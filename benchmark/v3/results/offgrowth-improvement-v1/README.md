# Off-growth enzyme discovery, second round: tier A results

Report: `docs/benchmarks/v3-offgrowth-improvement.md`.
Frozen plan: `tasks/v3-row12-offgrowth-improvement-2026-11-02.md`, sha256 in `PLAN_SHA256` (recorded before any
held-out score was read).

| File | What it holds |
|---|---|
| `PLAN_SHA256` | the frozen plan's hash and freeze date |
| `q_fixed_threshold_L0_six.{py,txt}` | the pre-registered fixed-precision search for the sequence-only comparator, run on the seen labels only |
| `tier-a-six/arms.tsv` | one row per group, arm, metric and label set, with bootstrap intervals and paired tests |
| `tier-a-six/cases.tsv` | one row per masked case and arm: gene, reaction and pathway recovery, the recovering test's fill size, the hidden protein's candidate rank |
| `tier-a-six/decisions.tsv.gz` | every distinct labelled decision of the unmasked runs, with its class |
| `tier-a-six/pathway_calls.tsv.gz` | every pathway call with its reactions and gene cluster |
| `tier-a-six/orphan_steps.tsv` | called reactions whose enzymes have no candidate above the score floor anywhere in the genome |
| `tier-a-six/cost.tsv` | per run and arm: admitted tests, linear programs, CPU seconds, fill sizes, stability, measured-panel fallback |
| `tier-a-six/summary.json` | the same numbers in one file, plus the additional-function list and the fill distribution |
| `tier-a-six/report.md` | the generated tables (captions only, no prose) |
| `tier-a-six/implementation-*.md` | the two implementation runs that made the evaluation fit in memory and fixed the measured-growth panel |
| `tier-a-six/supervisor_gate_check.{py,txt}` | independent check on real data that the call gate behaves as the plan prescribes |
| `tier-a-clean/`, `tier-b-clean/`, `tier-c-clean/` | the same tables for the three runs with our own CLEAN candidate tables (tier A re-run, and the two genome-level tiers) |
| `tier-a-clean/clean-inference-{report.md,receipt.json}` | how CLEAN was run on 47 proteomes on orix, with checkpoint hashes and per-tag counts |
| `tier-a-clean/operations-report-tiers-abc-clean.md` | the three waves and their scoring: job ids, node-hours, case counts, fallback counts, exact commands |
| `tier-a-clean/candidates-clean-checks.tsv` | per-proteome checks of the candidate tables built from the CLEAN output |
| `supervisor_number_check.{py,txt}` | re-derivation of every tier B and C number quoted in the report from the shipped tables |

Simulation outputs (361 MB of per-genome tars) stay on Shaheen at
`/scratch/project/k10018/gspa-v3/out/offgrowth2-v1/genomes/`, with their own sha256 files.
