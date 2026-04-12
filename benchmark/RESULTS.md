# GSPA 9-Genome Benchmark — Full Results

## Genomes

| Tag | Organism | Assembly | Domain | Proteins |
|---|---|---|---|---|
| ecoli | E. coli K-12 MG1655 | GCF_000005845.2 | Bacteria | 4,403 |
| hpylori | H. pylori 26695 | GCF_000008525.1 | Bacteria | 1,554 |
| mgenitalium | M. genitalium G37 | GCF_000027325.1 | Bacteria | 483 |
| mjannaschii | M. jannaschii DSM 2661 | GCF_000091665.1 | Archaea | 1,787 |
| ecolo157 | E. coli O157:H7 Sakai | GCF_000008865.2 | Bacteria | 5,056 |
| bsubtilis | B. subtilis 168 | GCF_000009045.1 | Bacteria | 4,288 |
| mtb | M. tuberculosis H37Rv | GCF_000195955.2 | Bacteria | 3,997 |
| synechocystis | Synechocystis sp. PCC 6803 | GCF_000009725.1 | Bacteria | 3,508 |
| paeruginosa | P. aeruginosa PAO1 | GCF_000006765.1 | Bacteria | 5,563 |

## Method

- **Leave-9-out DIAMOND + HMMER/Pfam** against Swiss-Prot (556k
  proteins after excluding all 9 genomes' accessions).
- **Noisy-OR integration** with correlation-group collapse (DIAMOND +
  Pfam in the homology group).
- **Four Bayesian priors** (EssentialityPrior, CoherencePrior,
  GapFillingPrior, GenomicContextPrior) using:
  - KEGG pathway DB (169 pathways, 5,262 reactions, 4,060 with GO terms)
  - Per-genome operons from intergenic-distance clustering on the
    RefSeq GFF (300 bp threshold)
  - Per-genome gapseq metabolic gaps (where available)
  - Essential function profiles (32 bacterial / 25 archaeal GO terms)
- **ConsistencyPrior** disabled in these results (requires per-genome
  NCBI taxonomy lineage; without it, the SAT checker over-penalizes).
- **Phase 8 DarkMatterSuggester** for metabolic-gap recovery via
  operon context.
- **Bootstrap F-max** (200 resamples, 95% CI) against dual ground truth:
  experimental-only GOA and full GOA (all evidence including IEA).

## F-max — Full-GOA truth (all evidence)

| Genome | GSPA (baseline) | GSPA (+priors) | PGAP | GSPA/PGAP |
|---|---|---|---|---|
| ecoli | 0.670 [0.662, 0.679] | **0.670** [0.662, 0.678] | — | — |
| ecolo157 | 0.835 [0.827, 0.845] | **0.835** [0.826, 0.844] | — | — |
| bsubtilis | 0.673 [0.659, 0.689] | **0.674** [0.659, 0.690] | — | — |
| mtb | 0.716 [0.705, 0.726] | **0.715** [0.705, 0.726] | — | — |
| synechocystis | 0.614 [0.599, 0.633] | **0.614** [0.600, 0.635] | — | — |
| paeruginosa | 0.601 [0.588, 0.616] | **0.598** [0.584, 0.613] | — | — |
| hpylori | 0.754 [0.730, 0.775] | **0.753** [0.730, 0.775] | 0.316 [0.298, 0.336] | **2.4×** |
| mgenitalium | 0.913 [0.897, 0.930] | **0.912** [0.897, 0.929] | 0.469 [0.446, 0.492] | **1.9×** |
| mjannaschii | 0.641 [0.625, 0.668] | **0.639** [0.624, 0.667] | 0.285 [0.267, 0.303] | **2.2×** |

### Coverage and IC-recall improvement from priors

F-max is a precision-recall trade-off that stays flat when priors push
borderline claims across the 0.5 threshold (gaining recall but also
some false positives). The priors' real contribution shows in
**coverage** (fraction of truth proteins with ≥1 prediction) and
**IC-weighted recall** (correctly-predicted terms weighted by
information content):

| Genome | Cov baseline | Cov +priors | ΔCov | icR baseline | icR +priors | ΔicR |
|---|---|---|---|---|---|---|
| ecoli | 0.740 | 0.744 | +0.004 | 0.511 | 0.513 | +0.002 |
| ecolo157 | 0.831 | 0.834 | +0.003 | 0.782 | 0.785 | +0.003 |
| bsubtilis | 0.618 | 0.625 | +0.007 | 0.480 | 0.489 | +0.009 |
| mtb | 0.720 | 0.721 | +0.001 | 0.591 | 0.592 | +0.001 |
| synechocystis | 0.618 | 0.633 | +0.015 | 0.432 | 0.443 | +0.011 |
| paeruginosa | 0.290 | 0.312 | +0.022 | 0.346 | 0.353 | +0.007 |
| hpylori | 0.528 | 0.535 | +0.007 | 0.607 | 0.609 | +0.002 |
| mgenitalium | 0.652 | 0.702 | **+0.050** | 0.659 | 0.677 | **+0.018** |
| mjannaschii | 0.208 | 0.276 | **+0.068** | 0.282 | 0.309 | **+0.027** |

mjannaschii (Archaea, the most phylogenetically divergent genome)
benefits the most: coverage jumps 33% relative. This matches the
expectation that priors are most valuable where homology evidence is
weakest.

### Prior activity per genome

| Genome | Essential (uncov.) | Coherence (pw-missing) | GapFill (fns) | Context (boosts) | Gaps available |
|---|---|---|---|---|---|
| ecoli | 8 | 1,764 | 219 | 802 | 368 |
| ecolo157 | 8 | 1,775 | 223 | 706 | 386 |
| bsubtilis | 14 | 1,906 | — | 459 | — |
| mtb | 10 | 2,042 | — | 644 | — |
| synechocystis | 15 | 1,888 | 190 | 260 | 277 |
| paeruginosa | 13 | 1,914 | 271 | 583 | 427 |
| hpylori | 14 | 1,740 | — | 204 | — |
| mgenitalium | 17 | 1,305 | — | 55 | — |
| mjannaschii | 10 | 1,618 | 135 | 120 | 124 |

Synechocystis required the NCBI `gene_refseq_uniprotkb_collab.gz`
mapping file because UniProt's `xref_refseq` API returned zero results
for this proteome (see "Reference data requirements" below).

Dashes (—) indicate gapseq re-runs still in progress for those
genomes; GapFillingPrior runs but finds no gaps to boost.

## Dark Matter Suggester (Phase 8)

The DarkMatterSuggester runs after integration converges and assigns
metabolic-gap functions to specific proteins via operon context. For
each gapseq-identified gap (pathway P missing reaction R with target
GO term f_R):

1. Score every operon with a Bayes factor BF(O, P) — "does this
   operon participate in pathway P?" — using the current posteriors
   of its members.
2. Within passing operons (BF ≥ 10), compute per-protein log-odds
   for f_R: L_likelihood + L_operon + commitment penalty.
3. Softmax over operon members → q(p).
4. If q(top) > 0.5: singleton suggestion ("this protein does it").
   Otherwise: disjunctive suggestion over the smallest top-k
   whose cumulative q ≥ 0.9.

| Genome | Gaps | Singleton | Disjunctive | Total |
|---|---|---|---|---|
| ecoli | 368 | **913** | 2,903 | 3,816 |
| ecolo157 | 386 | **1,095** | 3,516 | 4,611 |
| paeruginosa | 427 | **1,083** | 3,835 | 4,918 |
| mjannaschii | 124 | **193** | 377 | 570 |
| synechocystis | 277 | **314** | 1,015 | 1,329 |

Singleton suggestions are high-confidence: "protein X in operon Y
fills the gap for reaction R in pathway P." Each carries full
provenance: the Bayes factor, the per-protein log-odds decomposition,
and the softmax probability.

## Reference data requirements

### RefSeq → UniProt ID mapping

Operons are derived from RefSeq GFFs (protein_id = `WP_*` or `NP_*`)
but claims use UniProt accessions. The mapping between the two is
essential for GenomicContextPrior and DarkMatterSuggester.

**Problem**: the UniProt REST API's `xref_refseq` field is empty for
some proteomes (notably Synechocystis sp. PCC 6803), returning zero
mappings.

**Solution**: use NCBI's dedicated collaboration file
[`gene_refseq_uniprotkb_collab.gz`](https://ftp.ncbi.nlm.nih.gov/refseq/uniprotkb/gene_refseq_uniprotkb_collab.gz)
(~1.2 GB, 176M rows). Each row maps a RefSeq protein accession to a
UniProt accession with the method (`identical` or `similar`).
`build_refseq_uniprot_map.py` scans this file once and extracts
per-genome TSV mappings.

For Synechocystis, this recovered 2,913 mappings (2,905 identical)
where the API returned zero — fixing the operon/claim mismatch and
enabling GenomicContextPrior (260 boosts) and DarkMatterSuggester
(314 singleton + 1,015 disjunctive = 1,329 suggestions).

### Required downloads for a full GSPA benchmark run

| File | Source | Size | Purpose |
|---|---|---|---|
| `go.owl` | [GO consortium](http://purl.obolibrary.org/obo/go.owl) | 130 MB | GO hierarchy, class labels |
| `go-computed-taxon-constraints.obo` | [GO imports](http://current.geneontology.org/ontology/imports/go-computed-taxon-constraints.obo) | 2 MB | 13k only_in + 6k never_in constraints |
| `ec2go` | [GO external2go](http://current.geneontology.org/ontology/external2go/ec2go) | 350 KB | EC → GO mapping |
| `gene_refseq_uniprotkb_collab.gz` | [NCBI RefSeq](https://ftp.ncbi.nlm.nih.gov/refseq/uniprotkb/gene_refseq_uniprotkb_collab.gz) | 1.2 GB | RefSeq ↔ UniProt mapping |
| `goa_uniprot_all.gaf.gz` | [GOA FTP](https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz) | 15 GB | Ground truth (experimental + IEA) |
| `uniprot_sprot.fasta.gz` | [UniProt](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz) | 93 MB | Leave-N-out reference DB |
| KEGG ec-pathway link | [KEGG REST](https://rest.kegg.jp/link/pathway/ec) | 1 MB | Built into `kegg_pathways.tsv` via `build_kegg_pathway_tsv.py` |
| gapseq (conda) | `conda install -c bioconda gapseq` | 2 GB | Metabolic model + gap detection |

### ConsistencyPrior requires taxonomy lineage

The ConsistencyPrior loads 13,265 `only_in_taxon` + 5,896
`never_in_taxon` constraints from `go-computed-taxon-constraints.obo`.
However, without a taxonomy hierarchy file (`--taxonomy`) telling the
SAT solver "this genome is taxon X, which is a child of Y, which is a
child of Z," the solver treats every `only_in_taxon` constraint as a
potential violation — because it can't prove the genome IS in the
required taxon.

This caused F-max drops of up to 0.051 (ecolo157) when
ConsistencyPrior was enabled without taxonomy context. The fix: gate
ConsistencyPrior on `--taxonomy` being present. The SAT4J
infrastructure works correctly (314 conflicting GO terms identified in
E. coli, driven down to 47 through iterative penalization); it just
needs proper per-genome taxon context.

### Gapseq re-runs pending

Four genomes (hpylori, bsubtilis, mtb, mgenitalium) had their gapseq
`-all-Reactions.tbl` files corrupted to zero bytes by a GlusterFS
`sed -i` rename race condition. Re-runs are in progress from `/tmp`
(local disk). Until they complete, these genomes lack GapFillingPrior
and DarkMatterSuggester data.

## Ablation summary

| Config | Mean F-max (full-GOA) | Description |
|---|---|---|
| DIAMOND only | 0.652 | Single predictor, sequence similarity |
| Pfam only | 0.116 | Single predictor, domain-based |
| Combined | **0.684** | Noisy-OR integration, no priors |
| +4 priors | **0.682** | + essentiality, coherence, gap-filling, context |

Priors maintain F-max (±0.003) while consistently improving coverage
and IC-recall. The main F-max signal comes from DIAMOND homology;
Pfam adds +0.03-0.07 on larger genomes; priors add coverage depth
rather than threshold-crossing mass.
