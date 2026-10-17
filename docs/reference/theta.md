# Hyperparameters (theta)

The `--theta <thetaFile>` JSON file supplies reliability weights, prior strengths,
and coordinate ascent convergence parameters to `gspa integrate`. When `--theta`
is omitted, GSPA applies verified default values calibrated across benchmark
panels.

## Mathematical framework

GSPA evidence integration operates in two stages: Stage 0 per-protein Noisy-OR
pooling and Stage 1 genome-scale Bayesian iterative refinement.

### Stage 0 Noisy-OR pooling

Input claims from diverse predictors for a given protein $p$ and function $f$
are grouped into correlation groups according to `evidenceType.correlationGroup()`.
This prevents non-independent evidence sources (such as DIAMOND and eggNOG,
which query overlapping sequence homology databases) from artificially inflating
confidence.

Within each correlation group $g$, GSPA retains the single strongest claim:

$$s_g = \max_{c \in g} \left( w_{\text{type}(c)} \cdot p_c \right)$$

where $w_{\text{type}(c)} \in [0, 1]$ is the reliability hyperparameter for that
evidence type, and $p_c$ is the calibrated probability of claim $c$.

If `second_opinion_bonus` $> 0$ and the top two claims within group $g$
originate from genuinely distinct databases, a secondary support increment is
added:

$$s_g \leftarrow \min\left(0.999, s_g + \text{second\_opinion\_bonus} \cdot s_{\text{second}}\right)$$

Across independent correlation groups, probabilities are pooled via Noisy-OR:

$$p_{\text{pooled}} = 1 - \prod_{g} (1 - s_g)$$

The pooled probability is converted into log-odds space and clipped to numerical
stability bounds $[-12.0, 12.0]$:

$$L_{\text{lik}} = \ln \left( \frac{p_{\text{pooled}}}{1 - p_{\text{pooled}}} \right)$$

### Stage 1 coordinate ascent refinement

Refinement proceeds iteratively over the whole genome. At iteration 0, posteriors
are initialized from likelihood alone: $L^{(0)} = L_{\text{lik}}$.

At iteration $t$, the total prior log-odds boost for protein $p$ and function key
$k$ is evaluated against the current proteome state:

$$L_{\text{pri}}^{(t)} = \sum_{m} \lambda_m \cdot \text{Prior}_m\left(p, k, \text{State}^{(t-1)}\right)$$

where $\lambda_m$ is the prior weight from `prior_weights`.

The unrelaxed target log-odds is:

$$L_{\text{target}}^{(t)} = \text{clip}\left(L_{\text{lik}} + L_{\text{pri}}^{(t)}\right)$$

To guarantee numerical stability and prevent limit cycles when non-monotone priors
(such as consistency violations) flip status, Jacobi under-relaxation is applied
with damping parameter $d \in (0, 1]$:

$$L^{(t)} = (1 - d) \cdot L^{(t-1)} + d \cdot L_{\text{target}}^{(t)}$$

Iterations terminate when the mean absolute change in posterior probability
across all keys falls below threshold $\epsilon$:

$$\frac{1}{|K|} \sum_{k \in K} \left| \sigma(L_k^{(t)}) - \sigma(L_k^{(t-1)}) \right| < \epsilon$$

where $\sigma(x) = \frac{1}{1 + e^{-x}}$. If divergence is detected (two successive
iterations with increasing delta), the state rolls back to the last stable snapshot.

## JSON schema and parameter reference

A `theta.json` file accepts the top-level keys detailed below.

### Evidence reliability (`reliability`)

Map from `EvidenceType` string identifier to reliability weight $w \in [0, 1]$.
Controls the maximum contribution of a single predictor channel before Noisy-OR
combination.

| Key | Default | Biological rationale |
|---|---|---|
| `SEQUENCE_SIMILARITY` | `0.70` | Alignment against curated databases (DIAMOND, MMseqs2). |
| `SEQUENCE_DOMAIN` | `0.75` | Profile HMM signature matches (InterProScan, Pfam). |
| `STRUCTURE_SIMILARITY` | `0.80` | Tertiary structure alignment (FoldSeek). High independence from sequence. |
| `STRUCTURE_DEEPLEARNING` | `0.75` | Structural neural predictors (DeepFRI, GraphGOSeq). |
| `ORTHOLOGY` | `0.70` | Phylogenetic orthology transfers (eggNOG-mapper). |
| `SEQUENCE_DEEPLEARNING` | `0.65` | Sequence-based neural models (DeepGO-Plus, ProteInfer). |
| `PROTEIN_LM_EMBEDDING` | `0.65` | Protein language model embedding similarity (SaProt, GOPredSim). |
| `GENOMIC_CONTEXT` | `0.45` | Operon co-occurrence without direct homology. Conservative weight. |
| `METABOLIC_CONTEXT` | `0.55` | Metabolic pathway requirement from draft reconstructions. |
| `GENOMIC_LANGUAGE_MODEL` | `0.55` | Nucleotide language model scores over operon neighborhoods. |
| `LOCALIZATION` | `0.50` | Subcellular compartment filters (SignalP, DeepTMHMM). |
| `DOMAIN_SPECIFIC_AMR` | `0.85` | Curated antimicrobial resistance markers (AMRFinderPlus). |
| `DOMAIN_SPECIFIC_CAZY` | `0.80` | Curated carbohydrate-active enzyme models (dbCAN). |
| `DOMAIN_SPECIFIC_BGC` | `0.75` | Biosynthetic gene cluster functional units (antiSMASH). |
| `DOMAIN_SPECIFIC_VF` | `0.75` | Virulence factor database hits (VFDB). |
| `SEQUENCE_MOTIF` | `0.50` | Short linear motifs (PROSITE, ELM). |
| *Other unlisted types* | `0.60` | Fallback default for custom or newly registered evidence types. |

### Correlated group collapse

Evidence types are mapped to correlation groups to prevent double-counting:

| Correlation Group | Member Evidence Types |
|---|---|
| `homology` | `SEQUENCE_SIMILARITY`, `SEQUENCE_DOMAIN`, `ORTHOLOGY` |
| `ml_protein_seq` | `SEQUENCE_DEEPLEARNING`, `PROTEIN_LM_EMBEDDING` |
| `structure` | `STRUCTURE_SIMILARITY`, `STRUCTURE_DEEPLEARNING` |
| `context` | `GENOMIC_CONTEXT`, `METABOLIC_CONTEXT` |
| `ml_genomic` | `GENOMIC_LANGUAGE_MODEL` |
| `inferred_context` | `DARK_MATTER`, `REACTION_LOCAL_CONTEXT`, `CROSS_GENOME_TRANSFER`, `ML_RANKER` |
| `region_features` | `SEQUENCE_REGION_ML` |
| `viral` | `GENOMIC_REGION_ML` |
| `motif` | `SEQUENCE_MOTIF` |
| `localization` | `LOCALIZATION` |
| `domain_*` | `DOMAIN_SPECIFIC_AMR`, `DOMAIN_SPECIFIC_CAZY`, `DOMAIN_SPECIFIC_BGC`, `DOMAIN_SPECIFIC_VF` |

### Coordinate ascent parameters

| Key | Type | Default | Mathematical role |
|---|---|---|---|
| `max_iter` | integer | `6` | Maximum number of coordinate ascent iterations. |
| `epsilon` | float | `0.005` | Convergence threshold on mean absolute posterior change across all keys. |
| `damping` | float | `0.5` | Jacobi under-relaxation factor $d$ applied to log-odds updates. |
| `second_opinion_bonus` | float | `0.0` | Fractional boost added when distinct databases agree within the same correlation group. |

### Prior engine weights (`prior_weights`)

Map specifying scalar multiplier $\lambda_m$ applied to each active prior.
Defaults to `1.0` for all priors when omitted.

| Key | Default | Associated prior class |
|---|---|---|
| `essentiality` | `1.0` | `gspa.integration.prior.EssentialityPrior` |
| `coherence` | `1.0` | `gspa.integration.prior.CoherencePrior` |
| `consistency` | `1.0` | `gspa.integration.prior.ConsistencyPrior` |
| `anatomical` | `1.0` | `gspa.integration.prior.AnatomicalPrior` |
| `gap_filling` | `1.0` | `gspa.integration.prior.GapFillingPrior` |
| `genomic_context` | `1.0` | `gspa.integration.prior.GenomicContextPrior` |
| `homology_transfer` | `1.0` | `gspa.integration.prior.HomologyTransferPrior` |

### Individual prior hyperparameters

#### Essentiality prior

- `alpha_ess` (float, default `1.5`): Positive log-odds boost applied to candidate
  proteins whose functions match an uncovered domain-essential core function.
  Capped at $2 \times \alpha_{\text{ess}}$ per claim.

#### Coherence prior

- `alpha_coh` (float, default `1.0`): Positive log-odds boost applied to candidate
  proteins completing an unsatisfied biological process `has_part` dependency
  or metabolic pathway, weighted by the fraction of the module already present.
- `alpha_complex_boost` (float, default `0.85`): Weight factor boosting candidate
  partner subunits that complete an obligate heteromeric complex.
- `alpha_complex_penalty` (float, default `0.75`): Weight factor penalizing an
  isolated singleton subunit of an obligate heteromeric complex when no candidate
  partners exist.

#### Consistency prior

- `alpha_cons` (float, default `0.68`): Soft log-odds penalty subtracted from
  claims participating in a taxonomic constraint conflict (identified by SAT4J
  UNSAT cores).
- `consistency_hard_filter` (boolean, default `false`): When true, converts the
  soft penalty into an absolute filter (subtracts 1000.0 log-odds), removing the
  annotation entirely.

#### Anatomical prior

- `alpha_anat_penalty` (float, default `0.68`): Soft log-odds penalty subtracted
  from candidate GO claims whose target anatomical structure cannot occur in the
  clade of the organism.
- `alpha_anat_boost` (float, default `0.65`): Positive log-odds boost for sub-threshold
  candidate claims acting on an anatomical structure that is actively supported
  by multiple high-confidence proteins in the genome.
- `anatomical_hard_filter` (boolean, default `false`): When true, subtracts 1000.0
  log-odds from anatomically forbidden claims.

#### Gap filling prior

- `alpha_gap` (float, default `1.5`): Positive log-odds boost applied to candidate
  proteins matching a metabolic gap identified by draft metabolic reconstruction
  (e.g. gapseq). Guesses generated by gapseq itself receive a discounted boost of
  $0.7 \times \alpha_{\text{gap}}$.

#### Genomic context prior

- `alpha_ctx` (float, default `1.0`): Positive log-odds boost applied to candidate
  proteins in an operon where high-confidence members establish strong consensus
  for a shared metabolic pathway.
- `alpha_gap_ctx` (float, default `1.5`): Additional multiplier applied when the
  candidate function simultaneously satisfies an operon pathway consensus and
  fills an identified metabolic network gap.

#### Cross-genome homology transfer prior

- `alpha_homology` (float, default `1.5`): Maximum positive log-odds boost added to
  a protein when its orthogroup consensus strongly supports the function.
- `homology_min_consensus` (float, default `0.5`): Minimum consensus posterior
  probability required across orthologs before the prior fires.
- `homology_min_delta` (float, default `0.05`): Minimum positive difference between
  orthogroup consensus and individual posterior required to trigger a boost.
- `homology_max_boost` (float, default `2.0`): Absolute ceiling on the log-odds
  boost from homology transfer.

### Dark matter suggester (`dark_matter`)

Configures parameters for nominating new candidate functions from unannotated
operon open reading frames:

| Key | Type | Default | Description |
|---|---|---|---|
| `bf_min` | float | `10.0` | Minimum hypergeometric Bayes factor required to nominate an operon. |
| `gamma_in_p` | float | `50.0` | Likelihood elevation factor for functions belonging to the pathway. |
| `coverage_threshold` | float | `0.9` | Cumulative posterior credible set threshold for disjunctive sets. |

## Complete example `theta.json`

The following example configuration specifies custom reliabilities and prior weights:

```json
{
  "reliability": {
    "SEQUENCE_SIMILARITY": 0.72,
    "SEQUENCE_DOMAIN": 0.78,
    "STRUCTURE_SIMILARITY": 0.82,
    "STRUCTURE_DEEPLEARNING": 0.75,
    "ORTHOLOGY": 0.70,
    "SEQUENCE_DEEPLEARNING": 0.65,
    "PROTEIN_LM_EMBEDDING": 0.65,
    "GENOMIC_CONTEXT": 0.45,
    "METABOLIC_CONTEXT": 0.55,
    "GENOMIC_LANGUAGE_MODEL": 0.55,
    "LOCALIZATION": 0.50,
    "DOMAIN_SPECIFIC_AMR": 0.85,
    "DOMAIN_SPECIFIC_CAZY": 0.80,
    "DOMAIN_SPECIFIC_BGC": 0.75,
    "DOMAIN_SPECIFIC_VF": 0.75,
    "SEQUENCE_MOTIF": 0.50
  },
  "second_opinion_bonus": 0.10,
  "max_iter": 6,
  "epsilon": 0.005,
  "damping": 0.50,
  "prior_weights": {
    "essentiality": 1.0,
    "coherence": 1.0,
    "consistency": 1.0,
    "anatomical": 1.0,
    "gap_filling": 1.0,
    "genomic_context": 1.0,
    "homology_transfer": 1.0
  },
  "alpha_ess": 1.5,
  "alpha_coh": 1.0,
  "alpha_complex_boost": 0.85,
  "alpha_complex_penalty": 0.75,
  "alpha_cons": 0.68,
  "consistency_hard_filter": false,
  "alpha_anat_penalty": 0.68,
  "alpha_anat_boost": 0.65,
  "anatomical_hard_filter": false,
  "alpha_gap": 1.5,
  "alpha_ctx": 1.0,
  "alpha_gap_ctx": 1.5,
  "alpha_homology": 1.5,
  "homology_min_consensus": 0.5,
  "homology_min_delta": 0.05,
  "homology_max_boost": 2.0,
  "dark_matter": {
    "bf_min": 10.0,
    "gamma_in_p": 50.0,
    "coverage_threshold": 0.90
  }
}
```
