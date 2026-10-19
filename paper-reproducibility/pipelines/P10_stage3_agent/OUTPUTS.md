# P10 outputs

All files are written to `results/` of the repository (benchmark v2; the v1 `results/pillar4_*`
files are left untouched). Remote intermediates live under
`/data/hohndor/gspa/proteomes/pillar4_repro_2026-10-19/` on unimatrix01.

| File | Producer | Content |
| --- | --- | --- |
| `results/pillar4_v2_benchmark.jsonl` | `build_pillar4_dataset_v2.py assemble` | Benchmark v2: extractor records, de-identified, plus scoring fields (`expected_*`, `ground_truth_label`, `acceptable_protein_ids`, `benchmark_*`, `jar_target_label`) that `run_pillar4_v2_agent.py` strips before the model sees a record |
| `results/pillar4_v2_id_mapping.tsv` | `assemble` | Opaque token (P####, O####) -> genome, original UniProt / operon id; anomaly id -> genome, case id, construction variant |
| `results/pillar4_v2_leak_checks.tsv` | `assemble` | Model-visible schema identity TP vs negative, context field sets, score distributions (AUC, Mann-Whitney), identifier hygiene, 20-label check (jar label vs OBO name), selection notes |
| `results/pillar4_v2_resolutions_<cfg>_<seed>.jsonl` | `run_pillar4_v2_agent.py` | Agent decisions per anomaly (c0_42, c2_42, c3_42, c3_43, c3_44, c4_42) |
| `results/pillar4_v2_transcripts_<cfg>_<seed>.jsonl` | `run_pillar4_v2_agent.py` | Full transcripts: run header (model, sampling, tool schema, sidecar sha256), every record sent, every chat-completion request (all messages incl. tool results) and raw response with UTC timestamps, final resolution |
| `results/pillar4_v2_lro_resolutions_<cfg>_<seed>.jsonl` | `run_pillar4_v2_agent.py --deidentify` | Same six configurations on the unchanged LRO-50 cases (protein ids mapped back) |
| `results/pillar4_v2_transcripts_lro50_<cfg>_<seed>.jsonl` (+ `.idmap.json`) | idem | LRO-50 transcripts (model saw opaque tokens; idmap = token -> accession per case) |
| `results/pillar4_v2_evaluation.tsv` | `evaluate_pillar4_v2.py` | Columns of `pillar4_agent_evaluation.tsv` + `benchmark_version=v2`; agent configs, deterministic, threshold, supervised (LR/RF/GBDT/MLP x 3 pre-tool / 8 post-tool features x with/without guardrail), single-feature baselines |
| `results/pillar4_v2_per_category.tsv` | idem | Correct / lenient-correct / promotions / false promotions per category and label for every method |
| `results/pillar4_v2_ml_baselines.tsv` | idem | Baselines with benchmark-v2 metrics and LRO-50 Hit@1 / MRR (single-feature LRO rows, RLGC order, supervised LRO) |
| `results/pillar4_v2_multiseed_summary.tsv` | idem | c3 seeds 42 (T=0.0), 43 and 44 (T=0.25) and mean +/- population SD, both benchmarks |
| `results/pillar4_v2_mcnemar.tsv` | idem | Exact McNemar (two-sided binomial) of c3 seed 42 against every other method on both benchmarks |
| `results/pillar4_v2_lro_reranking.tsv` | idem | Per-case LRO-50 re-ranking for c3 seed 42 |
| `results/pillar4_v2_trigger_rates.tsv` | `trigger_rates_v2.py aggregate` | 88 proteomes: n_proteins, n_triggered, pct, anomalies per category (production cap 150 and uncapped), per-domain and overall summary rows |
