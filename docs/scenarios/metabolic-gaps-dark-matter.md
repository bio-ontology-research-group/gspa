# Metabolic gaps and dark matter

Use this workflow to identify and fill missing enzymatic steps in metabolic
pathways. In genome-scale metabolic reconstructions, pathways often contain
unassigned reactions (gaps) even when physiological context or neighboring
enzymes indicate the pathway is present. GSPA handles metabolic gaps at two
levels:

1. **Gap-filling prior (`gap_filling`)**: boosts existing low-confidence sequence
   evidence for proteins matching missing pathway reactions.
2. **Dark-matter and RLGC suggesters**: scan unannotated or hypothetical loci
   located in genomic proximity to pathway anchors, evaluating whether genomic
   co-localization and operon structure support assigning the missing catalytic
   function.
3. **Outer fixed-point loop**: iteratively runs Bayesian integration, scores
   candidate dark-matter loci, promotes the top candidates above a calibrated
   decision threshold, pins them as evidence floors, and repeats until network
   closure or fixed-point convergence.

## What you need

| Input | Required | Description and generation |
|---|---|---|
| `claims.jsonl` | yes | Initial predictor claims file |
| `gaps.jsonl` | yes | Metabolic gaps extracted from gapseq or pathway databases |
| `layout.tsv` | for RLGC | Locus coordinates (`protein_id contig start end strand`) built by `make_layout.py` |
| Reaction graph | for RLGC | ModelSEED reaction database (`seed_reactions.tsv`) |
| EC aliases | for RLGC | ModelSEED EC-to-reaction mapping (`seed_Enzyme_Class_Reactions_Aliases_unique.tsv`) |
| Currency metabolites | for RLGC | High-degree hub metabolites to exclude (`diffusion_mets.tsv`) |
| `go.owl` | recommended | Gene Ontology release file |
| `ec2go` | recommended | Mapping from EC numbers to GO molecular function terms |

## Step 1: extract metabolic gaps from gapseq

If you run `gapseq find` or `gapseq doall`, gapseq produces `-all-Pathways.tbl`
and `-all-Reactions.tbl`. Use `benchmark/parse_gapseq_gaps.py` to extract
unresolved reaction steps in predicted pathways:

```bash
uv run python benchmark/parse_gapseq_gaps.py \
  --pathways-tbl gapseq_out/MYISO-all-Pathways.tbl \
  --reactions-tbl gapseq_out/MYISO-all-Reactions.tbl \
  --ec2go ec2go.txt \
  --tag MYISO \
  --out gaps.jsonl
```

Expected output:

```text
  MYISO: 142 triggered pathways, 1850 reaction rows
  MYISO: wrote 87 gaps (64 with EC) -> gaps.jsonl
```

Each line in `gaps.jsonl` is a JSON record specifying the pathway, reaction,
EC number, and mapped GO term:

```json
{"pathway_id": "MetaCyc:PWY-6163", "reaction_id": "RXN-11213", "ec_number": "2.7.7.60", "go_term": "GO:0046415", "gapseq_guessed": false}
```

## Step 2: prepare the genome layout

The Reaction-Local Genomic Context (RLGC) suggester measures the physical
clustering of candidate genes relative to established pathway anchors on the
chromosome. Extract CDS coordinates from your GFF3 file using
`benchmark/panel/make_layout.py`:

```bash
uv run python benchmark/panel/make_layout.py \
  --gff prokka_out/MYISO.gff \
  --map refseq_to_uniprot.tsv \
  --out genome_layout.tsv
```

The resulting TSV contains five columns:

```text
protein_id  contig    start   end     strand
MYISO_00001 contig_1  100     1150    +
MYISO_00002 contig_1  1200    2450    +
MYISO_00003 contig_1  2500    3100    -
```

## Step 3: run integration with RLGC and dark-matter suggesters

To run a single integration pass with RLGC candidate generation:

```bash
gspa integrate \
  --claims claims.jsonl \
  --out integrated.tsv \
  --gaps gaps.jsonl \
  --operons operons.tsv \
  --taxon bacteria \
  --rlc-suggester \
  --reaction-graph reference_data/seed_reactions.tsv \
  --reaction-ec-aliases reference_data/seed_Enzyme_Class_Reactions_Aliases_unique.tsv \
  --diffusion-mets reference_data/diffusion_mets.tsv \
  --genome-layout genome_layout.tsv \
  --suggestions-out suggestions.tsv \
  --features-out rlgc_features.tsv \
  --lite
```

Key RLGC parameters:

- `--rlc-radius-k`: reaction graph breadth-first search radius (default: `2`).
  Controls how many reaction hops away from the gap anchor genes are searched.
- `--rlc-alpha`: decay factor across reaction graph hops (default: `0.5`).
- `--rlc-kernel-bandwidth`: Gaussian kernel bandwidth in base pairs for genomic
  density smoothing (default: `5000.0`).
- `--rlc-anchor-threshold`: minimum posterior probability required for a locus
  to serve as an RLGC pathway anchor (default: `0.3`).
- `--rlc-currency-pct`: degree percentile threshold above which common hub
  metabolites (such as ATP, water, H+) are treated as currency metabolites and
  pruned from graph traversal (default: `99.0`).
- `--features-out`: exports the numerical feature vector computed for each
  candidate locus (including distance to anchor, operon co-membership, and
  topological reaction distance).

## Step 4: run the outer fixed-point promotion loop

In complex metabolic pathways, assigning one missing enzyme may reveal
contextual support for an adjacent gap. The outer loop executes iterative cycles
of Bayesian integration, dark-matter candidate scoring, and promotion:

```bash
gspa integrate \
  --claims claims.jsonl \
  --out integrated_final.tsv \
  --gaps gaps.jsonl \
  --operons operons.tsv \
  --taxon bacteria \
  --dark-matter \
  --suggestions-out dark_matter_suggestions.tsv \
  --iterate-gapseq \
  --max-gapseq-iter 5 \
  --promotion-strategy greedy \
  --gapseq-q-base 0.50 \
  --gapseq-q-step 0.05 \
  --gapseq-q-cap 0.75 \
  --gapseq-tau-cover 0.50 \
  --gapseq-pin-promotions true \
  --lite
```

### Outer loop flags

- `--iterate-gapseq`: enables the outer fixed-point loop (requires `--dark-matter`).
- `--max-gapseq-iter`: maximum number of outer iterations (default: `5`).
- `--gapseq-q-base`: base threshold $q_0$ on suggestion probability required for
  promotion in the first outer round (default: `0.5`).
- `--gapseq-q-step`: threshold increment added per outer round (default: `0.05`),
  ensuring later iterations require higher confidence.
- `--gapseq-q-cap`: maximum upper bound for the rising promotion threshold
  (default: `0.75`).
- `--gapseq-tau-cover`: posterior probability threshold for evaluating whether a
  metabolic gap is satisfied and covered (default: `0.5`).
- `--gapseq-pin-promotions`: pins promoted candidates as posterior probability
  floors in subsequent iterations (default: `true`), preventing rollback of
  established assignments.

### Promotion strategies

Configure how candidates are selected from suggestions via `--promotion-strategy`:

| Strategy | Flag value | Behavior |
|---|---|---|
| All above threshold | `default` | Promotes every candidate exceeding the current iteration threshold $q_t$. |
| Greedy rank | `greedy` | Sorts candidates by posterior confidence and commits the highest-ranked conflict-free candidates. |
| Weighted MaxSAT | `maxsat` | Solves an exact SAT4J MaxSAT optimization balancing candidate log-likelihood against pathway conflicts. Set `--maxsat-coherence-bonus` (e.g. `0.5`) to reward joint assignments within the same pathway. |
| Beam search | `beam` | Explores a search tree of joint assignments with `--beam-width` (default: `5`) and `--beam-candidates-per-gap` (default: `3`). |

### Real output

```text
  Loaded 4120 claims
  Asserted organism taxon: NCBITaxon_2
  Reaction graph: 34125 reactions, 82 currency metabolites, 4812 ECs
  Genome layout: 3840 loci across 1 contig(s)
  Running outer loop (maxIter=5, qBase=0.5, qStep=0.05, pin=true, strategy=greedy, refinedBf=true)
  Outer loop iter 1: 87 gaps -> 14 promoted
  Outer loop iter 2: 73 gaps -> 6 promoted
  Outer loop iter 3: 67 gaps -> 1 promoted
  Outer loop iter 4: 66 gaps -> 0 promoted (fixed point reached)
  Outer loop: iter=4, fixedPoint=true, cascade=false, promoted_per_iter=[14, 6, 1, 0], gaps_per_iter=[87, 73, 67, 66]
  Wrote: integrated_final.tsv
  Suggestions: dark_matter_suggestions.tsv
```

## Step 5: interpret dark-matter outputs

The suggestions TSV (`--suggestions-out`) records candidate assignments evaluated
for missing gaps:

```text
kind       pathway_id          reaction_id  function_id  go_aspect  operon_id  bayes_factor  suggestion_score  protein_ids   q_values  provenance
singleton  MetaCyc:PWY-6163    RXN-11213    GO:0046415   MF         op_0012    14.2851       0.6842            MYISO_01428   0.6842    dark_matter_rlgc
```

- `kind`: `singleton` (single unassigned gene candidate) or `disjunctive`
  (alternative candidate genes in an operon or cluster).
- `bayes_factor`: contextual Bayes factor derived from operon co-expression,
  synteny, and genomic distance to pathway anchors.
- `suggestion_score`: calibrated marginal probability that the candidate locus
  fills the specified reaction.
- `provenance`: records the source of the suggestion (`dark_matter_rlgc` or
  `dark_matter_operon`).

In `integrated_final.tsv`, promoted loci show an explicit prior contribution in
the `priors_fired` column:

```text
protein_id   type  function_id  go_aspect  posterior_prob  likelihood_logodds  final_logodds  n_supporting  priors_fired
MYISO_01428  GO    GO:0046415   MF         0.7124          -1.1200             0.9074         1             gap_filling:1.500;dark_matter_promote:0.527
```

## Pitfalls

- **Missing EC mappings**: if your gaps file uses EC numbers that lack direct GO
  terms in `ec2go`, those gaps cannot be matched against GO evidence claims.
  Verify that `ec2go` is loaded when running `benchmark/parse_gapseq_gaps.py`.
- **Coordinate mismatch in layout**: locus tags in `genome_layout.tsv` must match
  the protein IDs in `claims.jsonl` exactly. If you use Prokka or Bakta GFFs, ensure
  `make_layout.py` maps the CDS attributes to the same ID namespace.
- **Unbounded outer loops**: if $q_{\text{base}}$ is set too low (for example $< 0.3$),
  spurious hypothetical proteins may be promoted and trigger cascades of false
  positive pathway completions. Maintain `--gapseq-q-base` at or above `0.50`.
