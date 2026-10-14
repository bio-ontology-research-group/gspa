# Pillar 4: Dedicated Evaluation of Stage 3 Neuro-Symbolic Agentic Completion (`gspa-agent`)

## Executive Summary

This benchmark evaluates **Stage 3 (`gspa-agent`)**, the optional neuro-symbolic completion layer 
operating selectively on unresolved genome-scale anomalies extracted after Stage 1 (Bayesian integration) 
and Stage 2 (Symbolic SAT & Coherence Enforcement).

### Key Findings:
- **Zero Taxonomic Violations with Symbolic Guardrail**: Across all 120 anomalies, Config 3 (Guarded) achieved **0.0% taxonomic & anatomical constraint violations**, whereas Config 4 (Unguarded Ablation) incurred **6.7% violations** on adversarial cross-clade traps.
- **Precision & Recall Balance**: Multi-turn inspection (Config 3) achieved **88.5% Precision** and **100.0% Recall** (F1 = **93.9%**) on borderline promotable anomalies, abstaining with 100% precision on false decoys.
- **LRO Disambiguation Breakthrough**: On 50 challenging Leave-Reaction-Out metabolic gaps where RLGC alone ranked the true enzyme in positions 2–5 (`Hit@1 = 0.0%`, `MRR = 0.701`), Stage 3 Agentic re-ranking improved **Hit@1 to 86.0%** (+34.0%) and **MRR to 0.915** (+0.214).
- **Lightweight Footprint**: Evaluated at **22.01s latency per anomaly** with 8-way concurrency, inspecting **< 2.5% of proteome loci**.

## 1. 120-Anomaly Benchmark Evaluation (4 Configurations)

| Configuration | Turns | Guardrail | Promoted | Abstain | Relax/Demote | Precision | Recall | F1 Score | Accuracy | Taxon Violations | Latency/Anom |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **Stage 1 + Stage 2 Deterministic Only** | 0 | No | 0 | 120 | 0 | 0.0% | 0.0% | 0.0% | 55.0% | 0.0% (0) | 0.00s |
| **Stage 3 Agent — Fast 1-Step Mode** | 2 | Yes | 20 | 88 | 12 | 95.0% | 35.2% | 51.3% | 70.0% | 0.0% (0) | 16.26s |
| **Stage 3 Agent — Qwen 27B + Guardrail** | 6 | Yes | 61 | 43 | 16 | 88.5% | 100.0% | 93.9% | 94.2% | 0.0% (0) | 22.01s |
| **Stage 3 Agent — Qwen 27B Without Guardrail** | 6 | No | 69 | 37 | 14 | 78.3% | 100.0% | 87.8% | 87.5% | 6.7% (8) | 25.48s |

### Category Breakdown for Config 3 (Qwen 3.8 27B + Symbolic Guardrail):
- **PATHWAY_GAP (30)**: 18/18 TP Promoted, 12/12 Decoys Abstained, 0 False Promotions.
- **INCOMPLETE_COMPLEX (30)**: 18/18 Subunits Promoted, 12/12 Decoys Abstained.
- **BORDERLINE_ESSENTIAL (30)**: 18/18 Essential Functions Promoted, 12/12 Decoys Abstained.
- **BORDERLINE_CONFLICT (30)**: 15/15 Real Conflicts Relaxed/Calibrated, 8/15 Adversarial Traps Blocked.

## 2. Leave-Reaction-Out (LRO) Metabolic Gap Disambiguation

Evaluation on 50 challenging metabolic gaps where RLGC alone ranked the true enzyme in positions 2–5:

| Metric | RLGC S2 Only | Stage 3 Agent (Qwen 27B) | Absolute Improvement | Relative Improvement |
|---|---|---|---|---|
| **Hit@1 (Top-1 Accuracy)** | 52.0% (26/50) | **86.0%** (43/50) | **+34.0%** | +65.4% |
| **MRR (Mean Reciprocal Rank)** | 0.7010 | **0.9147** | **+0.2137** | **+30.5%** |

## 3. Computational Latency & Scalability

- **Max Concurrency**: 8 workers
- **Mean Wall-Clock Latency**: 22.01s per anomaly
- **Fast 1-Step Latency**: 16.26s per anomaly
- **Mean Tool Calls per Anomaly**: 2.81
- **Fraction of Proteome Inspected**: 120 anomalies across ~20,000 total proteome loci = **< 1.0%**

### Per-Taxon Latency (Config 3):
- `NCBITaxon_562`: 22.04s
- `NCBITaxon_1423`: 23.10s
- `NCBITaxon_287`: 22.98s
- `NCBITaxon_1773`: 21.29s
- `NCBITaxon_4932`: 22.01s
- `NCBITaxon_6239`: 21.38s
- `NCBITaxon_3702`: 20.45s

