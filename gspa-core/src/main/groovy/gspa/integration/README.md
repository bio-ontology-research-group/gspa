# `gspa.integration` — evidence claim → integrated annotation

The novel part of GSPA. Predictors emit heterogeneous, possibly correlated
evidence; this package consumes those claims, weights them, and produces
the integrated annotation set consumed downstream (quality metrics, GFF
writers, the crossgenome machinery).

## Class map

### Claim surface

| Class | Role |
|---|---|
| `EvidenceClaim` | A single (protein, term, evidence-type, score, provenance) tuple. Everything below manipulates these. |
| `ClaimKey` | Hashable lookup key for a claim (protein, term) pair. |
| `ClaimProvenance` | Which predictor + which output row produced this claim. |
| `ClaimExtractor` | Parses one predictor's output rows → `EvidenceClaim`s, mapping tool name → `EvidenceType`. |

### Evidence typing + combination

| Class | Role |
|---|---|
| `EvidenceType` | Enum of evidence kinds (sequence similarity, domain, structure, deep-learning, orthology, genomic/metabolic context, localization, …). Each has a *correlation group*. |
| `EvidenceCombiner` | Noisy-OR over independent groups. Within a group, only the strongest claim enters the product — this prevents double-counting correlated tools (DIAMOND + eggNOG, DGP + ProteInfer, …). |
| `CalibrationTable` | Per-predictor Platt-style score rescaling before Noisy-OR. |

### Priors

| Class | Role |
|---|---|
| `Prior` | A prior score on (protein, term) from outside the predictor set — typically taxon-conditional. |
| `PriorEngine` | Orchestrates prior assembly: load from disk, apply per-genome taxonomy, expose to the combiner. |
| `prior/` | Specific prior implementations (taxon priors, sequence-cluster priors, etc.). |

### Iterative refinement

| Class | Role |
|---|---|
| `IntegrationState` | Mutable state carried across iterations (current posteriors, gap set, pending claims). |
| `IterativeRefiner` | Inner fixed-point loop: combine → recompute gaps → re-promote claims until convergence. |
| `OuterIterativeRefiner` | Outer loop: wraps the inner refiner with a coverage-gap source (e.g. gapseq pathways) that can propose new claims between inner iterations. |
| `promotion/` | Strategies for turning a `MetabolicGap` or `CoverageGap` into a new `EvidenceClaim` (greedy / beam / MaxSAT variants from Phase 10). |
| `ranker/` | Ordering heuristics used by the promotion strategies. |
| `suggester/` | "Dark-matter" suggestion — for proteins with no direct prediction, propose functions from panel / reaction context. |

### Gaps + metabolism

| Class | Role |
|---|---|
| `MetabolicGap` | A reaction or pathway step that is required but not covered by any protein's current annotations. |
| `GapKey` | Hashable identifier for a gap. |
| `GapRecomputer` | Recompute the gap set after posteriors change (between refinement iterations). |

### Cross-genome

| Subpackage | Role |
|---|---|
| `crossgenome/` | Panel-level signal: orthogroup consensus, cross-genome LR, convergent-function detection. Drives Phase 12. |

### Output

| Class | Role |
|---|---|
| `IntegratedAnnotationSet` | Final output: per-protein list of promoted `EvidenceClaim`s with combined posteriors. |
| `IntegrationWriter` | Serializes an `IntegratedAnnotationSet` to the integrated TSV consumed by quality metrics and the benchmark harness. |

## Mental model

1. **Extract** per-predictor output → `EvidenceClaim`s (ClaimExtractor).
2. **Calibrate** each predictor's scores (CalibrationTable).
3. **Combine** within correlation groups (EvidenceCombiner).
4. **Prior-weight** against taxon / cluster priors (PriorEngine).
5. **Refine** iteratively: promotion → gap recomputation → re-combine
   (IterativeRefiner / OuterIterativeRefiner).
6. **Write** the converged posteriors out (IntegrationWriter).

Everything else in this package is a plug-in for one of those steps.
