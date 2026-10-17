# Guarded agentic completion

## Question

When does tool-mediated language-model reasoning complement deterministic probabilistic and stoichiometric inference, how does it compare against supervised tabular rerankers, and does a deterministic symbolic guardrail eliminate cross-clade hallucinations on structured anomalies and ambiguous metabolic gaps?

Paper: Section 2.8, Table 7 (in latex `tab:agent`), Figure 5a to 5d, and Supplementary Table S14 (`paper/main.tex`).

## Datasets

### De-identified 120-anomaly benchmark

120 post-Stage-2 structured anomalies extracted from seven prokaryotic and eukaryotic reference genomes (*E. coli*, *B. subtilis*, *P. aeruginosa*, *M. tuberculosis*, *H. pylori*, *S. cerevisiae*, and *D. melanogaster*). All protein identifiers and gene symbols are masked, and candidate order is hash-shuffled to prevent model memorization:

- **54 true-positive promotable gaps:** Genuine pathway gaps and incomplete complex subunits where true candidate genes exist in the candidate pool.
- **36 score-matched decoys:** Non-cognate loci with matching domain roles placed as negative decoys where the agent should abstain.
- **30 borderline ontology conflicts:** 15 genuine cross-clade conflicts and 15 adversarial traps constructed to bait language models into hallucinating vertebrate anatomical terms (*heart morphogenesis*, *neural crest migration*) in bacterial hosts.

Source data: `results/pillar4_agent_evaluation.tsv`, `results/pillar4_agent_benchmark_120.jsonl`.

### 50 ambiguous metabolic gaps

50 unique-protein ambiguous leave-reaction-out gaps in *E. coli* K-12 where Stage 1 RLGC yields multiple candidate loci with nearly identical operon distances and domain scores. Source data: `results/pillar4_lro_reranking.tsv`.

### Multi-seed replication cohort

Three random seeds (seeds 42, 43, and 44) evaluated on the full guarded multi-turn configuration to quantify stochastic variance in tool selection, precision, recall, and candidate ranking. Source data: `results/pillar4_multiseed_summary.tsv`.

## Protocol

### Configurations evaluated

1. **Stage 1 + Stage 2 Deterministic Only:** Deterministic inference baseline; abstains on all unresolved anomalies.
2. **Posterior Threshold (\(\ge 0.35\)) + Guardrail:** Promotes candidates whose Stage 2 posterior exceeds 0.35, filtered by the symbolic guardrail.
3. **Logistic Regression (3 Pre-Tool Scalars):** Supervised classifier evaluated under 5-fold cross-validation using three pre-tool features (sequence odds, RLGC Bayes factor, pathway score).
4. **GBDT (3 Pre-Tool Scalars):** Gradient-boosted decision tree using the same three pre-tool features.
5. **Stage 3 Fast 1-Step Mode:** Single-turn language-model reasoning (maximum 2 turns, single tool invocation).
6. **Stage 3 Multi-Turn Without Guardrail:** Multi-turn tool calling using Qwen 3.8 27B without symbolic verification.
7. **Full Guarded Multi-Turn Stage 3 Agent (3-Seed Mean):** Multi-turn tool calling using Qwen 3.8 27B verified by the deterministic Symbolic Guardrail.
8. **Supervised Post-Tool Classifiers (5-Fold CV + Guardrail):** Logistic Regression, Random Forest, GBDT, and Multilayer Perceptron (MLP) trained on eight numerical features extracted by the agent's inspection tools (domain overlap, operon co-directionality, synteny score, substrate distance, complex partner count, etc.).

### Metrics

- **Binary metrics:** Precision, Recall, \(F_1\) score, and Resolution Accuracy on the 120 anomalies.
- **Action distribution:** True-positive promotions, false-positive promotions, abstentions, ancestral relaxations, and demotions.
- **Violation rate (%):** Percentage of accepted proposals violating taxonomic or anatomical axioms.
- **LRO-50 Hit@1 (%) and MRR:** Top-1 accuracy and mean reciprocal rank across the 50 ambiguous metabolic gaps.
- **Latency and efficiency:** Concurrency-normalized wall-clock seconds per anomaly at batch concurrency 8, raw HTTP latency, and mean tool calls per anomaly.

## Baselines

Deterministic Stages 1 and 2 alone, fixed posterior thresholding (\(\ge 0.35\)), pre-tool supervised classifiers (Logistic Regression, GBDT), fast 1-step mode, unguarded multi-turn agent, and post-tool supervised classifiers.

## Reproduce

The evaluation script `benchmark/unified/pillar4/evaluate_pillar4.py` re-evaluates all configurations from pre-recorded trajectory archives:

```bash
uv run --with numpy --with scipy --with scikit-learn \
    python benchmark/unified/pillar4/evaluate_pillar4.py \
        --skip-live-runs --output-dir results/
```

To run the live agent over an OpenAI-compatible vLLM endpoint:

```bash
uv run --with requests --with numpy \
    python benchmark/unified/pillar4/evaluate_pillar4.py \
        --endpoint http://localhost:8000/v1 \
        --model Qwen/Qwen2.5-32B-Instruct \
        --seeds 42,43,44 \
        --output-dir results/
```

## Result tables

- `results/pillar4_agent_evaluation.tsv`: columns `configuration`, `n_anomalies`, `tp_promoted`, `fp_promoted`, `abstentions`, `relaxations`, `demotions`, `precision`, `recall`, `f1`, `accuracy`, `violations_accepted`, `violation_rate_pct`, `mean_latency_sec`, `mean_tool_calls`. ([GitHub](https://github.com/bio-ontology-research-group/gspa/blob/main/results/pillar4_agent_evaluation.tsv))
- `results/pillar4_ml_baselines.tsv`: side-by-side comparison of pre-tool and post-tool supervised classifiers on the 120 anomalies and 50 LRO gaps. ([GitHub](https://github.com/bio-ontology-research-group/gspa/blob/main/results/pillar4_ml_baselines.tsv))
- `results/pillar4_multiseed_summary.tsv`: per-seed metrics and standard deviations across seeds 42, 43, and 44. ([GitHub](https://github.com/bio-ontology-research-group/gspa/blob/main/results/pillar4_multiseed_summary.tsv))
- `results/pillar4_lro_reranking.tsv`: per-case ranks, reciprocal ranks, tool calls, and agent rationales for all 50 ambiguous metabolic gaps. ([GitHub](https://github.com/bio-ontology-research-group/gspa/blob/main/results/pillar4_lro_reranking.tsv))

### 120-anomaly evaluation results

Values transcribed directly from `results/pillar4_agent_evaluation.tsv`:

| Configuration | Anomalies | TP Prom. | FP Prom. | Abstain | Relax. | Precision | Recall | \(F_1\) Score | Accuracy | Accepted Violations | Violation Rate | Latency / Tools |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Stage 1 + 2 Deterministic Only | 120 | 0 | 0 | 120 | 0 | 0.0000 | 0.0000 | 0.0000 | 0.5500 | 0 | **0.0%** | 0.00 s / 0.00 |
| Posterior Threshold (\(\ge 0.35\)) + Guard | 120 | 29 | 22 | 39 | 30 | 0.5686 | 0.5370 | 0.5524 | 0.7083 | 0 | **0.0%** | < 0.01 s / 0.00 |
| Logistic Regression (3 Pre-Tool) | 120 | 36 | 6 | 78 | 0 | 0.8571 | 0.6667 | 0.7500 | 0.8000 | 3 | 2.5% | < 0.01 s / 0.00 |
| GBDT (3 Pre-Tool) | 120 | 41 | 4 | 75 | 0 | 0.9111 | 0.7593 | 0.8283 | 0.8583 | 0 | **0.0%** | < 0.01 s / 0.00 |
| Stage 3 Fast 1-Step Mode | 120 | 9 | 0 | 90 | 4 | **1.0000** | 0.1667 | 0.2857 | 0.6250 | 0 | **0.0%** | 19.80 s / 5.01 |
| Stage 3 Multi-Turn Without Guardrail | 120 | 43 | 4 | 35 | 0 | 0.9149 | 0.7963 | 0.8515 | 0.8833 | 1 | 0.83% | 29.75 s / 5.07 |
| **Stage 3 Qwen 27B + Guardrail** | **120** | **46** | **1** | **33** | **15** | **0.9787** | **0.8519** | **0.9109** | **0.9250** | **0** | **0.0%** | **32.53 s / 5.74** |
| Log. Reg. (8 Tool Features + Guard) | 120 | 47 | 1 | 71 | 1 | 0.9792 | 0.8704 | 0.9216 | 0.9333 | 0 | **0.0%** | < 0.01 s / 0.00 |
| Random Forest (8 Tool Features + Guard) | 120 | 52 | 4 | 62 | 2 | 0.9286 | 0.9630 | 0.9455 | 0.9500 | 0 | **0.0%** | < 0.01 s / 0.00 |
| GBDT (8 Tool Features + Guard) | 120 | 49 | 1 | 70 | 0 | 0.9800 | 0.9074 | 0.9423 | 0.9500 | 0 | **0.0%** | < 0.01 s / 0.00 |
| MLP (8 Tool Features + Guard) | 120 | 52 | 1 | 67 | 0 | **0.9811** | **0.9630** | **0.9720** | **0.9750** | **0** | **0.0%** | < 0.01 s / 0.00 |

### Multi-seed replication across random seeds

Values transcribed directly from `results/pillar4_multiseed_summary.tsv`:

| Seed | Precision | Recall | \(F_1\) Score | Accuracy | Violation Rate | LRO-50 Hit@1 (%) | LRO-50 MRR |
|---|---|---|---|---|---|---|---|
| seed_42 | 0.9787 | 0.8519 | 0.9109 | 0.9250 | 0.00% | 82.00% | 0.8740 |
| seed_43 | 0.9592 | 0.8704 | 0.9126 | 0.9333 | 0.00% | 84.00% | 0.8890 |
| seed_44 | 0.9787 | 0.8519 | 0.9109 | 0.9250 | 0.00% | 82.00% | 0.8747 |
| **Mean \(\pm\) s.d.** | **0.9722 \(\pm\) 0.0092** | **0.8581 \(\pm\) 0.0087** | **0.9115 \(\pm\) 0.0008** | **0.9278 \(\pm\) 0.0039** | **0.00 \(\pm\) 0.00%** | **82.67 \(\pm\) 0.94%** | **0.8792 \(\pm\) 0.0069** |

### Ambiguous metabolic gap re-ranking (LRO-50)

Values transcribed from `results/pillar4_ml_baselines.tsv` and `results/pillar4_lro_reranking.tsv`:

| Method | Feature Scope | LRO-50 Hit@1 (%) | LRO-50 MRR | Anomaly 120 \(F_1\) | Anomaly 120 Accuracy |
|---|---|---|---|---|---|
| Posterior Threshold (\(\ge 0.35\)) | None | 34.0% | 0.5700 | 0.5524 | 0.7083 |
| Stage 1 RLGC Operon Kernel Alone | Operon Distance | 34.0% | 0.6020 | n/a | n/a |
| Logistic Regression (5-Fold CV) | 3 Pre-Tool Scalars | 94.0% | 0.9667 | 0.7500 | 0.8000 |
| GBDT (5-Fold CV) | 3 Pre-Tool Scalars | 90.0% | 0.9383 | 0.8283 | 0.8583 |
| **Stage 3 Guarded Agent (3-Seed Mean)** | **Tool-Mediated Interactive** | **82.67%** | **0.8792** | **0.9115** | **0.9278** |
| Logistic Regression (5-Fold CV) | 8 Post-Tool Features + Guard | 96.0% | 0.9767 | 0.9216 | 0.9333 |
| GBDT (5-Fold CV) | 8 Post-Tool Features + Guard | 92.0% | 0.9500 | 0.9423 | 0.9500 |
| MLP (5-Fold CV) | 8 Post-Tool Features + Guard | 92.0% | 0.9507 | 0.9720 | 0.9750 |

## Main findings

1. **Selective invocation bounds compute cost:** Across all 88 multi-kingdom reference proteomes, deterministic anomaly triggers select a mean of only 1.84% of proteome loci (fewer than 2.5% in every clade), ensuring language-model inference remains tightly bounded to ambiguous cases.
2. **Guarded reasoning outperforms pre-tool ML:** On the de-identified 120-anomaly benchmark, the zero-shot guarded agent achieves \(F_1 = 0.9115 \pm 0.0008\) (accuracy 0.9278), outperforming supervised classifiers restricted to pre-tool summary scalars (Logistic Regression \(F_1 = 0.7500\), GBDT \(F_1 = 0.8283\)).
3. **Necessity of the deterministic Symbolic Guardrail:** Ablating the guardrail admits cross-clade violations (0.83% violation rate, accepting an adversarial trap) and drops \(F_1\) to 0.8515. The guarded agent achieves an absolute 0.0% violation rate across all three random seeds, guaranteeing that accepted proposals satisfy all taxonomic and anatomical axioms.
4. **Post-tool classifiers versus zero-shot agent:** Classifiers trained on the eight post-tool numerical features achieve up to \(F_1 = 0.9720\) on binary anomaly promotion at sub-millisecond latency. However, on borderline ontology conflicts, the zero-shot agent proposes 13.7 valid ancestral GO relaxations, whereas fixed post-tool classifiers can only choose among pre-defined binary actions. Supervised post-tool classifiers are optimal when labeled anomaly data exist; zero-shot guarded reasoning provides label-free operation and multi-class semantic ancestor relaxation.
5. **Resolving ambiguous metabolic gaps:** On 50 ambiguous multi-candidate LRO gaps where Stage 1 RLGC ranks the true enzyme at Rank 1 in only 34.00% of cases (MRR 0.6020), Stage 3 guarded agent reasoning lifts Hit@1 to 82.67 \(\pm\) 0.94% (MRR 0.8792; \(p = 2.8 \times 10^{-7}\), exact McNemar test).
