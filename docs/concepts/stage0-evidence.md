# Per-protein evidence (Stage 0)

Stage 0 builds the per-protein starting distribution that all later stages refine. It converts the raw scores of one or more function predictors into calibrated probabilities, weights them by source reliability, collapses sources that share underlying evidence, and pools the remaining independent groups with Noisy-OR. The output is an anchored log-odds \( L^{(0)}(p, f) \) for every candidate pair of protein \( p \) and function \( f \).

Stage 0 is the foundation of GSPA. Every Stage 1 update is added to \( L^{(0)} \), so a protein with strong sequence or structural evidence keeps its high-confidence assignment, and the set of candidate functions that Stage 0 proposes defines the space in which later stages can act. In the published reference genomes about 85% of annotated proteins already had \( P_0 \ge 0.90 \); genome-scale inference concentrates on the remainder.

## Inputs

Stage 0 consumes *evidence claims*: one claim per (protein, function, source). [`gspa annotate`](../reference/cli-annotate.md) produces claims from the predictors it runs; [`gspa integrate`](../reference/cli-integrate.md) reads them from a JSONL file given with `--claims` (format in [Input and output formats](../reference/file-formats.md)). Each claim carries:

- `protein_id`, `function_type` (`GO` or `EC`), `function_id`, `source`, `raw_score` (required);
- `evidence_type`, `go_aspect`, `metadata`, and `calibrated_prob` (optional). If `evidence_type` is absent it is inferred from the source name; claims whose type cannot be resolved are dropped. If `calibrated_prob` is present it is used as is and the calibration step below is skipped for that claim.

Stage 0 can run on a single predictor (for example one enzyme predictor such as CLEAN, MAPred, DeepProZyme, EnzBERT, GraphEC, or TopEC) or on several complementary channels. The published configurations use up to seven channels:

| Channel | Example tools | Evidence type | Correlation group |
|---|---|---|---|
| Sequence alignment | DIAMOND against Swiss-Prot | `SEQUENCE_SIMILARITY` | homology |
| Profile HMMs and domain signatures | HMMER/Pfam with Pfam2GO, InterProScan | `SEQUENCE_DOMAIN` | homology |
| Protein language model nearest neighbours | ESM-2, ProtT5 \(k\)-NN transfer | `PROTEIN_LM_EMBEDDING` | ml_protein_seq |
| Neural sequence classifiers | ProteInfer, DeepGOWeb, DeepGO-SE | `SEQUENCE_DEEPLEARNING` | ml_protein_seq |
| Structural similarity | FoldSeek, ProstT5 3Di tokens | `STRUCTURE_SIMILARITY` | structure |

Other evidence types (orthology, genomic and metabolic context, genomic language models, motifs, localization, region- and genome-level predictors, and domain-specific tools such as AMRFinder, dbCAN, antiSMASH, and VFDB) enter the same machinery. The [DeepGO-PlusPlus scenario](../scenarios/deepgo-plusplus.md) shows how to use a stacked GO predictor as the base.

## Model

### Calibration

A raw score \( s(c) \in [0, 1] \) of claim \( c \) from source \( e(c) \) is mapped to a probability by Platt scaling and clipped away from 0 and 1:

\[
\tilde q(c) = \operatorname{clip}\bigl(\sigma(a_{e(c)}\, s(c) + b_{e(c)}),\ \varepsilon,\ 1 - \varepsilon\bigr), \qquad \varepsilon = 10^{-4},
\]

where \( \sigma \) is the logistic function. The default coefficients \( (a_e, b_e) \) shipped in the code are:

| Source | \( a_e \) | \( b_e \) |
|---|---|---|
| `diamond`, `mmseqs2` | 6.0 | -3.0 |
| `hmmer`, `pfam`, `antismash` | 4.0 | -1.5 |
| `interproscan`, `foldseek`, `dbcan` | 5.0 | -2.0 |
| `eggnog-mapper`, `eggnog` | 3.5 | -1.3 |
| `signalp`, `deeptmhmm` | 3.0 | -1.0 |
| `operon` | 3.0 | -1.2 |
| `gapseq` | 4.0 | -1.8 |
| `amrfinder` | 8.0 | -4.0 |
| any other source | 4.0 | -2.0 |

Sources not listed (for example language-model or neural predictors) use the fallback curve unless their claims carry a precomputed `calibrated_prob`, which is how the benchmark pipelines supply calibrated scores for those channels. Implementation: `gspa-core/.../integration/CalibrationTable.groovy`.

### Reliability weighting and correlation groups

Each calibrated probability is multiplied by a reliability weight \( r_{e} \in (0, 1] \) of its evidence type. Claims for the same pair \( (p, f) \) are then partitioned into *correlation groups* \( \mathcal{G}(p, f) \): evidence types that draw on overlapping information (for example DIAMOND, Pfam, and eggNOG all reflect sequence homology) share a group, so that the same signal is not counted several times as if it were independent.

Within each group only the strongest weighted claim is kept:

\[
q_g(p, f) = \operatorname{clip}\Bigl(\max_{c \in g} r_{e(c)}\,\tilde q(c) \;+\; \beta_{\text{so}}\, q_{g,(2)}(p, f),\ 10^{-3},\ 0.999\Bigr),
\]

where \( q_{g,(2)} \) is the second-strongest weighted claim in the group and \( \beta_{\text{so}} \) is an optional second-opinion bonus that applies only when the two claims come from different sources and different databases. It is 0 (disabled) by default.

Default reliability weights (`gspa-core/.../integration/EvidenceCombiner.groovy`) and groups (`EvidenceType.groovy`):

| Evidence type | \( r_e \) | Group |
|---|---|---|
| `SEQUENCE_SIMILARITY` | 0.70 | homology |
| `SEQUENCE_DOMAIN` | 0.75 | homology |
| `ORTHOLOGY` | 0.70 | homology |
| `SEQUENCE_DEEPLEARNING` | 0.65 | ml_protein_seq |
| `PROTEIN_LM_EMBEDDING` | 0.65 | ml_protein_seq |
| `STRUCTURE_SIMILARITY` | 0.80 | structure |
| `STRUCTURE_DEEPLEARNING` | 0.75 | structure |
| `GENOMIC_CONTEXT` | 0.45 | context |
| `METABOLIC_CONTEXT` | 0.55 | context |
| `GENOMIC_LANGUAGE_MODEL` | 0.55 | ml_genomic |
| `DARK_MATTER`, `REACTION_LOCAL_CONTEXT`, `CROSS_GENOME_TRANSFER`, `ML_RANKER` | 0.60 | inferred_context |
| `SEQUENCE_MOTIF` | 0.50 | motif |
| `LOCALIZATION` | 0.50 | localization |
| `SEQUENCE_REGION_ML` | 0.60 | region_features |
| `GENOMIC_REGION_ML` | 0.60 | viral |
| `DOMAIN_SPECIFIC_AMR` | 0.85 | own group |
| `DOMAIN_SPECIFIC_CAZY` | 0.80 | own group |
| `DOMAIN_SPECIFIC_BGC`, `DOMAIN_SPECIFIC_VF` | 0.75 | own group (one each) |

### Noisy-OR pooling and the anchored log-odds

Independent groups are pooled with Noisy-OR, and the result is converted to clipped log-odds:

\[
P_0(p, f) = \operatorname{clip}\Bigl(1 - \prod_{g \in \mathcal{G}(p, f)} \bigl(1 - q_g(p, f)\bigr),\ 10^{-6},\ 1 - 10^{-6}\Bigr),
\qquad
L^{(0)}(p, f) = \operatorname{clip}\Bigl(\log \tfrac{P_0(p, f)}{1 - P_0(p, f)},\ -12,\ 12\Bigr).
\]

\( L^{(0)} \) is the *anchor*: it is computed once and is never overwritten by later stages. Stage 1 adds prior shifts to it at every iteration.

## Hyperparameters

| Parameter | Default | How to change |
|---|---|---|
| Platt coefficients \( (a_e, b_e) \) | table above | supply `calibrated_prob` in the claims file |
| Reliability \( r_e \) | table above | `reliability` object in the `--theta` JSON, keyed by evidence-type name |
| Second-opinion bonus \( \beta_{\text{so}} \) | 0.0 | `second_opinion_bonus` in the `--theta` JSON |
| Log-odds clip | \([-12, 12]\) | fixed |

See [Hyperparameters (theta)](../reference/theta.md) for the JSON schema.

## Outputs and audit trail

For every pair \( (p, f) \) Stage 0 records the list of supporting claims (with source, raw score, and calibrated probability) and the likelihood log-odds \( L^{(0)} \). Both are carried through all later stages: in the output TSV of `gspa integrate` as the columns `likelihood_logodds` and `n_supporting`, and in the optional `--provenance` JSON as `likelihood_logodds` and `supporting_sources`. Running `gspa integrate --enable-priors none` returns the Stage 0 posterior alone, which is the baseline row in the published ablations.

## Limitations

- Stage 0 defines the candidate set. A function that no predictor proposes for a protein cannot be recovered by Stage 1 priors on that protein; the leave-reaction-out benchmarks report recovery both unconditionally and within the Stage 0 candidate pool for this reason (see [Benchmarks](../benchmarks/index.md)).
- Taking the maximum within a correlation group is conservative: two genuinely independent tools assigned to the same group contribute only once.
- The default Platt coefficients are provisional curves for the listed tools, not learned per dataset. For well-calibrated results with new predictors, supply calibrated probabilities in the claims file or override the reliability weights.
