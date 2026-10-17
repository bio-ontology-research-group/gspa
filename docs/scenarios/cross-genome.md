# Cross-genome evidence

Use this workflow to share functional evidence across multiple related genomes
in a clade, genus, or pangenome cohort. When annotating a single isolate in
isolation, lineage-specific or divergent genes may lack direct matches in public
databases. However, when multiple genomes from the same group are analysed
together, shared genomic context and conservation patterns provide strong
corroborating signals.

GSPA supports two complementary cross-genome mechanisms:

1. **Orthogroup consensus (`homology_transfer`)**: transfers high-confidence
   functional consensus across orthologous protein clusters.
2. **Cross-genome reaction-locus catalog re-scoring**: computes an empirical
   conditional likelihood ratio (LR) statistic across a multi-genome panel,
   measuring how consistently an orthogroup occurs in the genomic neighborhood
   of a metabolic reaction across genomes.

## What you need

| Input | Required | Description |
|---|---|---|
| Proteome FASTAs | yes | Protein FASTA files for all genomes in the cohort |
| MMseqs2 | yes | High-speed sequence clustering tool (`mmseqs`) |
| Baseline integrated TSVs | yes | Initial per-genome `integrated.tsv` from Stage 1 |
| Genome layouts | for catalog | Locus layout TSVs for each genome (from `make_layout.py`) |
| Reaction database | for catalog | ModelSEED reaction graph and EC aliases |

## Step 1: cluster cohort proteins into orthogroups

Concatenate all protein sequences from your cohort into a single FASTA file,
ensuring sequence headers contain unique protein IDs:

```bash
cat genomes/*.faa > all_cohort_proteins.faa
```

Cluster proteins using MMseqs2:

```bash
mmseqs easy-cluster all_cohort_proteins.faa cohort_clu /tmp/mmseqs_tmp \
  --min-seq-id 0.5 \
  -c 0.8 \
  --cov-mode 0 \
  --threads 16
```

MMseqs2 writes `cohort_clu_cluster.tsv` containing `cluster_rep_id <TAB> member_id`.
Convert this into GSPA's orthogroups format (`protein_id <TAB> cluster_id`):

```bash
uv run python benchmark/orthogroups_to_map.py \
  cohort_clu_cluster.tsv \
  orthogroups.tsv
```

The resulting `orthogroups.tsv` maps every protein ID to its representative
cluster ID:

```text
MYISO_00001    REP_00124
MYISO_00002    REP_00891
OTHER_00412    REP_00124
```

## Step 2: build orthogroup consensus

Calculate the consensus functional annotation for each orthogroup across the
baseline single-genome integration runs using `benchmark/build_orthogroup_consensus.py`:

```bash
uv run python benchmark/build_orthogroup_consensus.py \
  --orthogroups orthogroups.tsv \
  --integrated runs/genome1_integrated.tsv \
  --integrated runs/genome2_integrated.tsv \
  --integrated runs/genome3_integrated.tsv \
  --min-cluster-members 2 \
  --min-consensus 0.40 \
  --out cluster_consensus.tsv
```

The output file `cluster_consensus.tsv` records the mean posterior probability
per function for each cluster:

```text
cluster_id  type  function_id  consensus_prob
REP_00124   GO    GO:0004022   0.8840
REP_00124   EC    1.1.1.1      0.8412
REP_00891   GO    GO:0005524   0.9120
```

### Running integration with homology transfer

To integrate an individual genome using the orthogroup consensus prior:

```bash
gspa integrate \
  --claims claims/genome1.jsonl \
  --out integrated_cross.tsv \
  --orthogroups orthogroups.tsv \
  --cluster-consensus cluster_consensus.tsv \
  --enable-priors essentiality,coherence,consistency,homology_transfer \
  --taxon bacteria \
  --lite
```

Proteins in `genome1` receive a prior boost in `priors_fired` corresponding to
the consensus confidence of their orthogroup:

```text
protein_id   type function_id posterior_prob likelihood_logodds final_logodds n_supporting priors_fired
MYISO_00001  GO   GO:0004022  0.8624         -0.1200            1.8340        1            homology_transfer:1.954
```

## Step 3: compile a reaction-locus catalog

When predicting metabolic gap assignments, genomic co-localization in a single
genome can sometimes be accidental. The Reaction-Locus Catalog measures the
cross-genome conservation of physical synteny between an orthogroup $C$ and
reaction $R$.

### Mathematical definition

For every pair of orthogroup $C$ and reaction $R$ across a panel of genomes,
we count:

- $n_{\text{sig,with}}$: number of genomes where reaction $R$'s anchor signature
  is present and a member of $C$ lies in the genomic window around the anchor.
- $n_{\text{sig,total}}$: number of genomes where $R$'s anchor signature is present.
- $n_{\text{base,with}}$: number of genomes containing any member of $C$.
- $n_{\text{base,total}}$: total number of genomes in the panel.

The conditional likelihood ratio is:

$$\mathrm{LR}(C, R) = \frac{p(C \text{ in window} \mid R \text{ present})}{p(C \text{ present} \mid \text{baseline})} = \frac{n_{\text{sig,with}} / n_{\text{sig,total}}}{n_{\text{base,with}} / n_{\text{base,total}}}$$

A Jeffreys-prior beta-binomial 90% credible interval is computed on $\log \mathrm{LR}$.
Entries whose credible interval spans $\mathrm{LR} = 1.0$ (no significant enrichment)
are filtered out.

### Compiling the catalog

Build the catalog with `benchmark/cross_genome/build_catalog.py`:

```bash
uv run python benchmark/cross_genome/build_catalog.py \
  --orthogroups orthogroups.tsv \
  --reaction-graph reference_data/seed_reactions.tsv \
  --reaction-ec-aliases reference_data/seed_Enzyme_Class_Reactions_Aliases_unique.tsv \
  --diffusion-mets reference_data/diffusion_mets.tsv \
  --panel-manifest panel_manifest.tsv \
  --out rxn_locus_catalog.tsv
```

The resulting `rxn_locus_catalog.tsv` has the following schema:

```text
# panel_size=48
orthogroup_id  reaction_id  n_sig_with  n_sig_total  n_base_with  n_base_total
REP_00124      rxn00148     38          40           42           48
REP_00891      rxn01200     12          14           15           48
```

## Step 4: re-score RLGC suggestions with catalog evidence

When running the Reaction-Local Genomic Context (RLGC) suggester during
`gspa integrate`, pass `--rxn-locus-catalog`:

```bash
gspa integrate \
  --claims claims.jsonl \
  --out integrated_catalog.tsv \
  --gaps gaps.jsonl \
  --genome-layout genome_layout.tsv \
  --rlc-suggester \
  --reaction-graph reference_data/seed_reactions.tsv \
  --reaction-ec-aliases reference_data/seed_Enzyme_Class_Reactions_Aliases_unique.tsv \
  --diffusion-mets reference_data/diffusion_mets.tsv \
  --rxn-locus-catalog rxn_locus_catalog.tsv \
  --cg-lambda 1.0 \
  --cg-min-support 3 \
  --cg-require-credible true \
  --suggestions-out catalog_suggestions.tsv \
  --lite
```

### Parameters

- `--rxn-locus-catalog`: path to `rxn_locus_catalog.tsv`.
- `--cg-lambda`: exponent applied to the cross-genome LR in the posterior update (default: `1.0`).
- `--cg-min-support`: minimum $n_{\text{sig,total}}$ required for a pair $(C, R)$
  to be trusted (default: `3`).
- `--cg-require-credible`: drops LRs whose 90% credible interval includes 1.0
  (default: `true`).

### Real output

```text
  Reaction graph: 34125 reactions, 82 currency metabolites, 4812 ECs
  Genome layout: 3840 loci across 1 contig(s)
  RLGC suggester emitted 42 suggestions
  ReactionLocusCatalog: 18450 (C, R) entries; panel=48
  Cross-genome rescoring complete; 42 suggestions after rescore
  Suggestions: catalog_suggestions.tsv
```

## Pitfalls

- **Inconsistent protein identifiers**: sequence headers across cohort FASTA
  files must not collide. Use distinct prefixes per genome (such as `GENOME1_0001`
  and `GENOME2_0001`).
- **Insufficient panel size**: compiling a Reaction-Locus Catalog on fewer than
  5 to 10 genomes results in wide credible intervals that are discarded by
  `--cg-require-credible`. For small datasets, rely primarily on the
  `homology_transfer` prior.
- **Over-transferring specific functions**: if cluster members share an enzyme
  superfamily fold but have diverged in substrate specificity, single-member
  outliers can bias the consensus. Keep `--min-consensus` at or above `0.40`.
