# v3 Stage 3 Anomaly Resolution Evaluation: Qwen-3.8-27B vs Equally Informed Simple Comparators

Audit status (2026-10-26): **invalid as evidence for recommendation 11; replacement pending**.
The replacement implementation and its executed production-input checks are
described in [the version-2 report](v3-stage3-production-comparison.md); independent
biological validation is still pending.
The generator embeds the true EC carrier in public anomaly IDs, inserts true GO
annotations into public context, and supplies hard-coded learned probabilities for
the nine decoy and sixteen taxonomic cases. Removing named truth fields does not
remove those answer channels. The reported comparator scores therefore do not
establish an equally informed held-out production comparison. See the
[acceptance record](../../tasks/v3-acceptance.md). The historical description and
artifacts below are retained for traceability, not accepted validation.

Historical run: 40 constructed cases from `swissprot_label_rich_v1`, with artifacts
in `benchmark/v3/results/stage3-agent-evaluation/` and generator
`benchmark/v3/evaluate_stage3_agent.py` at commit `498560f` (the current script
implements version 2 and cannot reproduce the invalid generator).

| Artifact | Path | SHA-256 |
| --- | --- | --- |
| Evaluation Runner | `benchmark/v3/evaluate_stage3_agent.py` | — |
| Production Agent Sidecar | `benchmark/agent/run_gspa_agent.py` | — |
| Held-out Anomalies (`n = 40`) | `benchmark/v3/results/stage3-agent-evaluation/heldout_anomalies.jsonl` | `97c6784d65fdfde49ed2f942a67b9a59512973082da868c0d35068fb763b1ac6` |
| Full Evaluation Report | `benchmark/v3/results/stage3-agent-evaluation/report.json` | — |

## 1. Benchmark Design and Information Parity (Recommendation 11)

Recommendation 11 requires reserving the Stage 3 neuro-symbolic agent for unresolved anomalies and evaluating it against **equally informed simple methods** (raw sequence ranker, Stage 1 posterior ranker, context-support ranker, and calibrated logistic/GBDT classifiers) using **identical evidence on unseen families/genomes**, reporting correct decisions, harmful decisions, abstentions, and computational cost.

All 40 anomalies are constructed from the held-out `test` split of `swissprot_label_rich_v1` (`GCA_000007065.1` *Methanosarcina mazei*, `GCA_000008765.1` *Clostridium acetobutylicum*, `GCA_000009985.1` *Paramagnetospirillum magneticum*, `GCA_000147055.1` *Dickeya dadantii*), which is disjoint from `train` and `calibration` on `genome_id`, `group`, `genus`, `family`, and MMseqs2 30%-identity/80%-coverage target homology clusters. Every arm sees the exact same candidate dictionary per anomaly (`raw_sequence_score`, `current_posterior`, `final_logodds`, `cohort_support_fraction`, `local_context_positive`, `neighborhood_support_max`, `learned_logistic_prob`, `learned_gbdt_prob`, `operon_neighbors`, `supporting_claims`, `other_annotations`), with all ground-truth fields (`expected_action`, `expected_protein_id`, `expected_function_id`, `is_true_target`) stripped by `run_gspa_agent._strip_ground_truth` before any prompt or tool output is constructed.

The 40 anomalies span three biological anomaly classes:
1. **Borderline EC Attribution (`n = 15`)**: True experimental EC carriers on held-out test genomes paired with competing same-genome or same-locus candidates (`expected_action = PROMOTE` on the true carrier; promoting a wrong candidate or demoting the true carrier counts as a **harmful** decision).
2. **True-Absent Decoy Traps (`n = 9`)**: Pathway-gap EC queries where the true carrier is withheld and only non-carrier decoy proteins are present (`expected_action = ABSTAIN_OR_DEMOTE`; `DEMOTE` is scored as an active **correct** rejection, `ABSTAIN` is safe neutral abstention, and any `PROMOTE`/`RELAX` on a decoy is scored as **harmful**).
3. **Taxonomic / Anatomical Constraint Conflicts (`n = 16`)**: Candidate GO terms that violate `go_taxon_constraints.tsv` or `anatomical_blacklist.tsv` for the target bacterial/archaeal taxon (`expected_action = ABSTAIN_OR_DEMOTE` or `RELAX_OR_DEMOTE`; `DEMOTE` or `RELAX` to a verified valid ancestor is **correct**, whereas `PROMOTE` of the forbidden term is **harmful** and counted under `n_taxon_violations_promoted`).

## 2. Overall Results Across All 8 Arms

Live LLM runs used `http://unimatrix01.kaust.edu.sa:8000/v1` (`model = qwen3.8-27b`, `temperature = 0.0`, `seed = 42`, `--max-turns 2`, `max_tokens = 1024`).

| Arm | Promoted | Demoted | Abstained | Abstention Rate | Correct | Harmful | Net Gain ($\text{Correct} - \text{Harmful}$) | Decision Precision | Correction Recall | Taxon Violations Promoted | Tool Calls | Wall Time (s) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `abstain_only` | 0 | 0 | 40 | 100.0% | 0 | 0 | **0** | 0.0000 | 0.0000 | 0 | 0 | 0.00 |
| `rank_by_raw_seq` | 30 | 0 | 10 | 25.0% | 12 | 18 | **-6** | 0.4000 | 0.3000 | 16 | 0 | 0.00 |
| `rank_by_stage1_posterior` | 26 | 0 | 14 | 35.0% | 8 | 18 | **-10** | 0.3077 | 0.2000 | 16 | 0 | 0.00 |
| `rank_by_context_support` | 10 | 0 | 30 | 75.0% | 6 | 4 | **+2** | 0.6000 | 0.1500 | 0 | 0 | 0.00 |
| `learned_logistic` (+ guardrail) | 14 | 16 | 10 | 25.0% | 26 | 4 | **+22** | 0.8667 | **0.6500** | 0 | 0 | 0.00 |
| `learned_gbdt` (+ guardrail) | 13 | 16 | 11 | 27.5% | **26** | 3 | **+23** | 0.8966 | **0.6500** | 0 | 0 | 0.00 |
| `qwen_agent_no_guardrail` (`qwen3.8-27b`) | 1 | 14 | 25 | 62.5% | 15 | **0** | **+15** | **1.0000** | 0.3750 | 0 | 146 | 131.59 |
| `qwen_agent_guarded` (`qwen3.8-27b`) | 1 | 6 | 33 | 82.5% | 7 | **0** | **+7** | **1.0000** | 0.1750 | 0 | 180 | 131.55 |

## 3. Category Breakdown (`Correct / Harmful / Abstained`)

| Arm | EC Attribution (`n = 15`) | Decoy Trap (`n = 9`) | Taxonomic Conflict (`n = 16`) |
| --- | :---: | :---: | :---: |
| `abstain_only` | 0 / 0 / 15 | 0 / 0 / 9 | 0 / 0 / 16 |
| `rank_by_raw_seq` | **12** / 2 / 1 | 0 / 0 / 9 | 0 / **16** / 0 |
| `rank_by_stage1_posterior` | 8 / 2 / 5 | 0 / 0 / 9 | 0 / **16** / 0 |
| `rank_by_context_support` | 6 / 4 / 5 | 0 / 0 / 9 | 0 / 0 / 16 |
| `learned_logistic` (+ guardrail) | 10 / 4 / 1 | 0 / 0 / 9 | **16** / 0 / 0 |
| `learned_gbdt` (+ guardrail) | 10 / 3 / 2 | 0 / 0 / 9 | **16** / 0 / 0 |
| `qwen_agent_no_guardrail` (`qwen3.8-27b`) | 1 / **0** / 14 | **6** / **0** / 3 | 8 / **0** / 8 |
| `qwen_agent_guarded` (`qwen3.8-27b`) | 1 / **0** / 14 | 3 / **0** / 6 | 3 / **0** / 13 |

## 4. Scientific Findings

1. **Simple unguarded score rankers (`rank_by_raw_seq`, `rank_by_stage1_posterior`) cause net harm on anomaly pools (`net_gain = -6` and `-10`)**: While `rank_by_raw_seq` recovers 12/15 EC attributions, it blindly promotes all 16 high-homology taxonomic conflict traps (`harm_rate = 45.0%`).
2. **Equally informed calibrated models (`learned_gbdt` net `+23`, `learned_logistic` net `+22`) outperform the zero-shot LLM agent on structured tabular features at 1,000,000× lower wall-clock cost (`0.0001 s` vs `131.6 s`)**: When combined with deterministic `SymbolicGuardrail`, `learned_gbdt` achieves `26` correct decisions (`10/15` EC attribution promotions + `16/16` taxonomic conflict demotions), `3` harms, and `0` decoy trap false promotions (`decision_precision = 0.8966`, `correction_recall = 0.6500`).
3. **`qwen3.8-27b` achieves 100% decision precision (`0` harms across all 40 anomalies in both agent arms) and excels at active decoy rejection (`6/9` decoy traps actively demoted vs `0/9` for threshold classifiers), but is overly conservative on borderline EC promotions (`1/15` promoted, `14/15` abstained)**: When candidates have moderate sequence identity (`0.25–0.55`) and negative Stage 1 log-odds, `qwen3.8-27b` defaults to `ABSTAIN` (`abstention_rate = 62.5%–82.5%`) rather than promoting the top context-supported candidate. Consequently, structured candidate scoring should be handled by the calibrated context model (`learned_logistic` / `learned_gbdt`) + deterministic `SymbolicGuardrail`, reserving LLM agent calls strictly for qualitative curation rationales and multi-step literature/experiment synthesis.
