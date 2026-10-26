# v3 condition-planner evaluation (equal budget, measured cases, actual Qwen)

Implementation: `benchmark/v3/run_planner_benchmark.py` (runner, scoring, statistics), `benchmark/v3/make_planner_report.py`
(renders `REPORT.md`), `benchmark/v3/tests/test_planner_benchmark.py` (26 tests), `benchmark/v3/planner_evaluation.sbatch`
(prepared for the root, not submitted). Results: `benchmark/v3/results/planner-evaluation/` (`protocol.json`, `sweep.json`,
`deterministic.jsonl`, `qwen.shard*.jsonl`, `qwen.shard*_transcripts.jsonl.gz`, `report.json`, `REPORT.md`).
The production engine `scripts/metabolic_experiments.py` is **unmodified**; every arm runs through it.

## What the planner decides, and what that can and cannot show

Cases: the 61 measured Keio knockout / carbon-source rows of the Maranas-lab 2026 auxotroph validation (44 loci, 24
substrates), protocol `manuscript_regulatory` (ATPM 8.39, O2 20, carbon 10, substrate-specific regulatory closures), the
same model and manifests as `run_metabolic_benchmark.py`. Prepared data: `docs/benchmarks/v3-metabolic-data.md`.

**The measured set cannot define a measured discovery task.** All 61 rows were chosen because iML1515 predicted the deletion
abolishes growth on that substrate, and for most loci exactly one substrate was tested. A condition menu that lists only
the measured substrate of each gene, or any shortlist built from which rows exist, would hand the planner the answer.
So the menu is the full locus x substrate grid and the planner never sees which pairs are measured:

* menu per gene: for each of the 24 substrates a matched wild-type control (`carbon_source:<exchange>`) and a deletion under
  that carbon source (`knockout:H:<locus>:<exchange>@<exchange>`), 48 items;
* label-free evidence per (gene, substrate): a candidate record (`transport`) and, from a model-structure graph distance
  (no LP, undirected, high-degree metabolites excluded, no labels), `network_near` (<= 4 reaction hops) and `network_reachable`
  (<= 10). Thresholds were fixed before any run. Which substrate a gene was measured on is not an input.

Therefore two experiments, never mixed:

* **Experiment A, model-derived condition discovery (development).** Success = the arm itself ran a matched control that
  grows and a deletion that does not (both completed solves). The reference is an exhaustive evaluator sweep of the same model
  (48 LPs per gene, outside any planner budget). A is circular with respect to function: it says whether a planner finds
  model-predicted conditional essentiality cheaply. It is not measured function recovery and no gene function is recovered.
* **Experiment B, measured phenotype comparison on visited conditions.** For visited (gene, substrate) pairs that have a
  measured label, the model's deletion prediction is compared with the measurement after the run. Because the 61 rows were
  chosen by a model-based design (the design model and the planner model coincide), this estimates precision of predicted
  auxotrophy on whatever the arm visited. It supports no recall or specificity claim and is not independent of the model.
  Partial growth is a separate class and is never counted as confirmation. The exhaustive reference for all 61 cases is reported
  with the same rules, and cross-checked against `run_metabolic_benchmark.py` (same engine) and the independent COBRApy replay
  (`check_metabolic_cobra.py`) that establish solver mechanics.

## Frozen protocol

Recorded in `protocol.json` with the SHA-256 of cases, engine and runner before any planner is run.

| item | value |
|---|---|
| budgets | 8 and 16 experiments after the baseline; LP cap = budget + 1 for every arm (`max_solver_calls`), 1 LP per experiment |
| focus / FVA / product LPs | disabled (`max_focus_reactions = 0`, `max_baseline_focus_reactions = 0`): no forced-flux checks, no product objectives |
| arms | `fixed_paired`, `random_item` (10 seeds), `random_paired` (10 seeds), `evidence_item`, `evidence_paired`, `qwen` |
| random seeds | 0-9, frozen |
| Qwen | `qwen3.8-27b`, `http://unimatrix01.kaust.edu.sa:8000/v1`, temperature 0, seed 0, thinking off, max 300 output tokens, 2 invalid choices then evidence fallback (production engine loop), timeout 600 s |
| statistics | bootstrap over genes (10,000, seed 20261026), exact McNemar on discordant genes, Wilson intervals for proportions |

Arms:

* `fixed_paired`: control then knockout for the first budget/2 substrates in alphabetical exchange order (gene-independent).
* `random_item`: the engine's `random` strategy over all 48 items. It is not pair-aware.
* `random_paired`: seeded random condition order, each visited as control + knockout.
* `evidence_item`: the engine's `evidence` strategy unmodified. Ties break by item id, so every control precedes every
  knockout and the arm spends its budget on controls. This is production behaviour, reported as is.
* `evidence_paired`: the same `evidence_key`, applied per condition, visited as control + knockout.
* `qwen`: the engine's own selection loop (`run_experiments(strategy="qwen")`), with the transport (`UsageChat`, same
  payload as `make_qwen_chat` plus token accounting) and the prompt text replaced. The prompt adds, per menu item, the content of the
  cited evidence records (type and hop count) and the goal; the evidence arms order by exactly these records.
  Qwen only chooses an id; nothing it returns is executed.

Paired arms are an explicit `panel` of the `fixed` strategy; item-level arms use the engine strategies. The 24 closures are
passed as `condition_constraints` and apply only when their substrate is open, as in the registered cases.

## Rules that keep the comparison honest

* **No truth in menus, evidence, hypotheses or prompts.** `public_view` keeps only the locus list, the substrate panel and the
  closures. Labels are joined after each run (`measured_visits`). Tests flip every label and measured value and check menus,
  evidence, prompts and run records are identical, and that dropping a measured pair does not change any gene's menu.
* **Information parity.** Each Qwen run audits its first prompt: identical menu ids and identical evidence content as the
  manifest the evidence arms use (`prompt_parity`, per run in the results). The audit function is tested to fail when evidence
  is hidden or items are dropped. Parity is of evidence and menu; the goal sentence is part of the task definition, and the
  deterministic arms realise the goal through pairing and the evidence key. This is stated, not hidden.
* **Equal budget.** All arms: same cap, same per-experiment cost. Unused LPs (early stop, skipped items) are reported, never reallocated.
  Qwen requests do not consume solver calls (engine policy); they are reported separately as API requests, tokens and seconds, as is wall time of every arm.
* **Qwen failures are visible.** Endpoint errors, malformed or non-menu choices, fallback activation and early stops are counted per gene.
  A Qwen discovery requires both its control and its knockout to be Qwen-selected; fallback-selected items never count. The hybrid
  success (any source), the Qwen-attributable success and the success among genes with zero invalid choices and no fallback are all reported.
* **Circular objectives.** Planner menus reject product (`bgc`) evidence and any focus check (a growth-floor forced flux). A growth-coupled
  product is flagged by the engine as `secretion_unavoidable_in_model`; tests show it is never independent evidence
  (`independent_proof_of_function = false`, no probability, no promotion) and never a planner discovery.
* **Every case reported.** All 44 genes and all 61 registered cases are kept with eligibility and reason; genes without any model-predicted
  essential condition stay in the tables and are only excluded from the Experiment A success denominator.

## Qwen endpoint behaviour observed

A smoke run on one gene before the evaluation showed 87-115 s per request regardless of prompt size (a 23-token prompt took 87 s; a 6.6k-token
menu prompt about 112 s). The timeout was therefore raised from the engine's 180 s to 600 s before any evaluated Qwen request
(recorded under `revisions` in `protocol.json`). Four shard workers keep concurrent load at four requests. No service setting was touched and no GPU job was run.

## Secondary-metabolism feasibility audit

MIBiG 4.0 (606 ready rows after the registered eligibility rule) records characterised production as **positive observations only**:
no negatives, no titres, no medium, no mapping of compounds to model metabolites. It therefore supports route-existence or recall
checks and cannot support an empirical medium-effect claim, a precision, or a specificity. No negative was synthesised from absence.

A traceable condition-dependent source was fetched: Westhoff et al., "Competition sensing changes antibiotic production in Streptomyces"
(Dryad 10.5061/dryad.wpzgmsbkw, Zenodo record 4397644, CC0; files and SHA-256 in
`results/planner-evaluation/secondary_condition_dependent_audit.json`, produced by `audit_condition_dependent_production.py`).
It has the same Streptomyces panel under a low and a high resource regime (pairwise growth inhibition matrices), carbon-source
utilisation (Biolog) and BiG-SCAPE BGC families per strain. It measures **bioactivity**, not a named compound: no compound identity, no
titre, no curated genome-scale models for these soil isolates, and no-inhibition is not evidence of non-production. It can anchor descriptive
statements about resource dependence of antibiotic activity; it cannot validate a flux-based medium effect on a specific product, so none
was attempted. Rejected after inspection: Zenodo 11065062 (F. verticillioides OSMAC, 4.8 GB raw UHPLC-HRMS, no model, no processed table).
No production dataset with measured negatives, named compounds and explicit media together was found.

## Empirical Results (`benchmark/v3/results/planner-evaluation/REPORT.md`)

All 44 registered loci (`41` with $\ge 1$ model-essential condition among the 24 substrates; `3` without: `b3564`, `b3599`, `b3605`) completed across all 6 arms under both budgets ($B = 8$, LP cap $9$; and $B = 16$, LP cap $17$). Every arm faced the identical menu and context per gene (`88` gene-budget cells checked, `0` mismatches, `prompt_parity_all_identical = True` at both $B=8$ and $B=16$, `over_cap_runs = 0`).

### 1. Experiment A — Model-Derived Condition Discovery ($n = 41$ Eligible Genes)

| Budget | Arm | Runs | Genes Scored | Success Rate | 95% CI (Bootstrap) | Mean Yield | Mean LPs | LP Cap | Wall Time (s) |
|---:|---|---:|---:|---:|---|---:|---:|---:|---:|
| **$B=8$** | **`qwen` (`qwen3.8-27b`, strict Qwen-attributable)** | 44 | 41 | **85.4%** (`35/41`) | `[73%, 95%]` | **81.7%** | 8.30 | 9 | 28,572.8 |
| **$B=8$** | `evidence_paired` | 44 | 41 | `39.0%` (`16/41`) | `[24%, 54%]` | `35.4%` | 9.00 | 9 | 16.6 |
| **$B=8$** | `random_paired` (10 seeds) | 440 | 41 | `19.5%` | `[16%, 24%]` | `17.4%` | 9.00 | 9 | 222.9 |
| **$B=8$** | `evidence_item` | 44 | 41 | `17.1%` (`7/41`) | `[7%, 29%]` | `17.1%` | 9.00 | 9 | 16.5 |
| **$B=8$** | `fixed_paired` | 44 | 41 | `12.2%` (`5/41`) | `[2%, 22%]` | `8.5%` | 9.00 | 9 | 18.5 |
| **$B=8$** | `random_item` (10 seeds) | 440 | 41 | `1.2%` | `[0%, 2%]` | `0.7%` | 9.00 | 9 | 204.9 |
| **$B=16$** | **`qwen` (`qwen3.8-27b`, strict Qwen-attributable)** | 44 | 41 | **80.5%** (`33/41`) | `[68%, 93%]` | **78.3%** | 15.41 | 17 | 37,338.1 |
| **$B=16$** | `evidence_paired` | 44 | 41 | `75.6%` (`31/41`) | `[61%, 88%]` | `71.5%` | 17.00 | 17 | 26.1 |
| **$B=16$** | `evidence_item` | 44 | 41 | `51.2%` (`21/41`) | `[37%, 66%]` | `51.2%` | 17.00 | 17 | 25.8 |
| **$B=16$** | `random_paired` (10 seeds) | 440 | 41 | `39.8%` | `[33%, 47%]` | `36.7%` | 17.00 | 17 | 276.4 |
| **$B=16$** | `fixed_paired` | 44 | 41 | `31.7%` (`13/41`) | `[17%, 46%]` | `27.6%` | 17.00 | 17 | 24.5 |
| **$B=16$** | `random_item` (10 seeds) | 440 | 41 | `8.0%` | `[5%, 11%]` | `6.9%` | 17.00 | 17 | 256.6 |

- **Pairwise comparisons vs. `qwen`:**
  - At **$B=8$**, `qwen` significantly outperforms every deterministic comparator: `+46.3%` over `evidence_paired` (`20` Qwen-only vs. `1` comparator-only, exact McNemar $p = 0.0000$), `+68.3%` over `evidence_item` (`28` vs. `0`, $p = 0.0000$), and `+73.2%` over `fixed_paired` (`31` vs. `1`, $p = 0.0000$). Among the `24` clean genes with zero invalid choices or fallbacks at $B=8$, `qwen` achieves **100.0%** (`1.000`) discovery.
  - At **$B=16$**, `qwen` achieves **80.5%** strict Qwen-attributable discovery (**87.8%** hybrid discovery including fallback steps; **84.6%** on clean genes), compared to **75.6%** for `evidence_paired` (`+4.9%`, `8` Qwen-only vs. `6` `evidence_paired`-only, exact McNemar $p = 0.7905$), **51.2%** for `evidence_item` (`+29.3%`, $p = 0.0075$), and **31.7%** for `fixed_paired` (`+48.8%`, $p = 0.0001$).
  - Why does `qwen`'s strict Qwen-attributable rate dip slightly from `85.4%` at $B=8$ to `80.5%` at $B=16$ while its hybrid rate rises from `85.4%` to `87.8%`? Under the longer 16-step horizon, more runs encounter $\ge 2$ endpoint timeouts or invalid menu choices (`28` genes activate fallback at $B=16$ vs. `9` at $B=8$), and our strict accounting rule refuses to credit `qwen` for any condition where either the control or knockout was chosen by `evidence_fallback`.

### 2. Experiment B — Measured Phenotype Comparison on Visited Conditions

| Budget | Arm | Visited Measured Pairs | With Growing Control | Binary Evaluable | Binary Agree | Agreement (95% Wilson) | Predicted No-Growth | Confirmed No-Growth | Confirmed (95% Wilson) |
|---:|---|---:|---:|---:|---:|---|---:|---:|---|
| **$B=8$** | **`qwen`** | **39** | **33** | **29** | **21** | **72.4%** `[54%, 85%]` | 38 | **21** | `55.3%` `[40%, 70%]` |
| **$B=8$** | `evidence_paired` | 16 | 16 | 14 | 10 | `71.4%` `[45%, 88%]` | 14 | 10 | `71.4%` `[45%, 88%]` |
| **$B=8$** | `evidence_item` | 7 | 7 | 7 | 6 | `85.7%` `[49%, 97%]` | 7 | 6 | `85.7%` `[49%, 97%]` |
| **$B=8$** | `fixed_paired` | 8 | 8 | 6 | 3 | `50.0%` `[19%, 81%]` | 6 | 3 | `50.0%` `[19%, 81%]` |
| **$B=8$** | `random_paired` (10 seeds) | 11.2 | 11.2 | 7.5 | 4.4 | `58.7%` (pooled) | 7.6 | 4.1 | `53.9%` (pooled) |
| **$B=8$** | `random_item` (10 seeds) | 8.6 | 0.8 | 6.0 | 3.2 | `53.3%` (pooled) | 5.3 | 3.1 | `58.5%` (pooled) |
| **$B=16$** | **`qwen`** | **36** | **34** | **29** | **22** | **75.9%** `[58%, 88%]` | 34 | **22** | `64.7%` `[48%, 79%]` |
| **$B=16$** | `evidence_paired` | 34 | 34 | 28 | 17 | `60.7%` `[42%, 76%]` | 29 | 17 | `58.6%` `[41%, 74%]` |
| **$B=16$** | `evidence_item` | 21 | 21 | 17 | 13 | `76.5%` `[53%, 90%]` | 18 | 13 | `72.2%` `[49%, 88%]` |
| **$B=16$** | `fixed_paired` | 17 | 17 | 15 | 10 | `66.7%` `[42%, 85%]` | 14 | 10 | `71.4%` `[45%, 88%]` |
| **$B=16$** | `random_paired` (10 seeds) | 22.9 | 22.9 | 16.2 | 9.0 | `55.6%` (pooled) | 15.8 | 8.3 | `52.5%` (pooled) |
| **$B=16$** | `random_item` (10 seeds) | 20.1 | 5.6 | 13.5 | 6.8 | `50.4%` (pooled) | 13.3 | 6.5 | `48.9%` (pooled) |

- **Exhaustive 61-case reference (no planner budget):** `61` eligible (`8` growth, `36` no_growth, `17` partial); binary agreement `26/44` (`59.1%`, 95% Wilson `[44%, 72%]`); confirmed no-growth `25/44` (`56.8%`). Cross-check against `auxotroph61-manuscript_regulatory-report.json`: `max_abs_knockout_growth_difference = 0.0`, `0` differing phenotype calls.

### 3. Artifact Digests (SHA-256)

- `benchmark/v3/results/planner-evaluation/protocol.json`: `29f2c8e8c6f1efc51e4131e4dd1f97920a79ce0d1dfb18478902bbb2976634ab`
- `benchmark/v3/results/planner-evaluation/sweep.json`: `9304be72aefd4f967de799daa62b1224dbbf5115a7d696df22e8b0ff10104a7e`
- `benchmark/v3/results/planner-evaluation/deterministic.jsonl`: `8336bc6018a2587ded3eddfe99bb26e4ea225904263d7a48bfc14684983f77cb`
- `benchmark/v3/results/planner-evaluation/report.json`: `9be97ec5fcd647576369f5d4e9afee56a9532bbdb8e37c27d3a380dc805f0d41`
- `benchmark/v3/results/planner-evaluation/REPORT.md`: `a3e756e575de6983b08e414ba8c21078001aebf753140b8a10b4e79bf55abb00`

## Reproduce

```bash
PY=/tmp/gspa-v3-testenv/bin/python
D=/tmp/gspa-v3-metabolic-data
$PY -I -m unittest discover -s benchmark/v3/tests -p "test_planner_benchmark.py"
$PY -I benchmark/v3/run_planner_benchmark.py --cases $D/derived/cases_manuscript_regulatory.jsonl --model $D/raw/iML1515.xml \
    --out benchmark/v3/results/planner-evaluation --stages sweep,deterministic          # 2 CPUs, minutes
$PY -I benchmark/v3/run_planner_benchmark.py ... --stages qwen --shard 0/4              # actual endpoint calls; repeat for 1/4..3/4
$PY -I benchmark/v3/run_planner_benchmark.py ... --stages report --phenotype-reference <run_metabolic_benchmark report.json>
$PY -I benchmark/v3/make_planner_report.py --dir benchmark/v3/results/planner-evaluation
```

The key is read only by `ME.load_api_key()` (environment or `~/.config/borg-llm/api-key`); nothing in this package prints or stores it.
