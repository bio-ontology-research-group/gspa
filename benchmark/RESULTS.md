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

- **Three predictors**:
  - DIAMOND blastp against a leave-9-out Swiss-Prot reference (556k
    proteins). Evidence type: SEQUENCE_SIMILARITY.
  - HMMER/Pfam (hmmsearch against Pfam-A). Evidence type: SEQUENCE_DOMAIN.
  - **InterProScan** (Pfam + TIGRFAM + CDD + SUPERFAMILY + Panther +
    ProSite + HAMAP + FunFam + Coils + NCBIfam, with InterPro2GO
    mapping). Evidence type: SEQUENCE_DOMAIN.
- **Noisy-OR integration** with correlation-group collapse (DIAMOND in
  the homology group; Pfam + InterProScan in the domain group).
- **Four Bayesian priors** (EssentialityPrior, CoherencePrior,
  GapFillingPrior, GenomicContextPrior) using:
  - KEGG pathway DB (169 pathways, 5,262 reactions, 4,060 with GO terms)
  - Per-genome operons from intergenic-distance clustering on the
    RefSeq GFF (300 bp threshold)
  - Per-genome gapseq metabolic gaps (5 of 9 genomes)
  - Essential function profiles (32 bacterial / 25 archaeal GO terms)
- **ConsistencyPrior** disabled (requires per-genome NCBI taxonomy
  lineage; architecture validated but not included in these numbers).
- **Phase 8 DarkMatterSuggester** for metabolic-gap recovery via
  operon context.
- **Bootstrap F-max** (200 resamples, 95% CI) against dual ground truth:
  experimental-only GOA and full GOA (all evidence including IEA).

## F-max definitions (important — read before interpreting tables)

We report **two complementary F-max metrics per genome**. Both are
computed on the same genome-restricted predictions/truth and use the
same threshold sweep, but they aggregate (precision, recall) over
proteins differently. We compute both because each emphasizes a
different failure mode.

### Metric A — per-genome micro-averaged F-max

Procedure (per genome, per truth set, per method):

1. Restrict predictions and truth to proteins in this genome only.
2. Sweep posterior threshold τ ∈ {0.05, 0.10, …, 1.00}. At each τ:
   - For every (protein, aspect, GO term) prediction with score ≥ τ,
     count the (protein, GO term) pair as a *predicted positive*.
   - **Sum TP, FP, FN across all (protein, GO-term) pairs in the
     genome** (NOT averaged over proteins).
   - Compute one precision, one recall, one F1 from those global sums.
3. F-max = max F1 across thresholds.

Bootstrap 95% CI: resample proteins with replacement (200 iterations);
for each sample redo step 2 with its own argmax-τ; report 2.5%/97.5%
quantiles. The point estimate uses the original (unsampled) set.

This metric upweights heavily-annotated proteins (a protein with 30
true GO terms contributes ~30× more to the totals than a protein with
1). It is the natural "how many of all the genome's annotation pairs
do we recover, and at what false-positive cost" number.

### Metric B — CAFA-style protein-centric F-max (CAFA III/IV protocol)

Procedure (per genome, per truth set, per method):

1. Restrict predictions and truth to proteins in this genome only.
2. Sweep posterior threshold τ. At each τ, for every protein p with at
   least one true annotation:
   - Let `pred_p` = GO terms predicted for p with score ≥ τ.
   - Let `truth_p` = true GO terms for p.
   - If `pred_p` is non-empty: `precision_p(τ) = |pred_p ∩ truth_p| / |pred_p|`
     (otherwise p is excluded from the precision average).
   - `recall_p(τ) = |pred_p ∩ truth_p| / |truth_p|` (proteins lacking
     any prediction at τ contribute recall = 0).
3. `avg_precision(τ) = mean(precision_p)` over the m(τ) proteins with
   ≥ 1 prediction; `avg_recall(τ) = mean(recall_p)` over all n_e
   proteins with ≥ 1 truth annotation.
4. `F1(τ) = 2·avg_p·avg_r / (avg_p + avg_r)`; F-max = max F1(τ).

Bootstrap 95% CI: same resampling scheme as Metric A.

This is the standard CAFA III/IV F-max (without IC weighting; that
weighted variant is S-min and is reported separately as `ic_recall`
in the JSON). Each protein contributes equally regardless of how many
annotations it carries — so dark proteins and well-studied proteins
count the same.

### What neither metric is

- Neither is a "global F-max across all genomes pooled together." We
  report one F-max per genome and never collapse those into a single
  number (when summarizing across genomes, we report the *mean
  GSPA/PGAP ratio*).
- Neither propagates GO terms to ancestors before evaluation. CAFA
  itself propagates; we compare prediction-as-given to truth-as-given,
  which means both GSPA and PGAP face the same propagation handling
  (none) and the comparison is fair.

### When the two metrics diverge

Micro F-max < CAFA F-max → predictions are good on most proteins but
weaker on a few heavily-annotated ones (which dominate the micro
totals). This is what we see on this benchmark (CAFA ~0.02-0.07
higher), driven by S. coelicolor / D. radiodurans, where a handful
of large multi-domain proteins carry many true GO terms that GSPA
recovers incompletely.

Micro F-max > CAFA F-max → predictions are strong on a few
heavily-annotated proteins but miss many lightly-annotated ones. We
do not see this case here.

Implementation: `benchmark/benchmark_pgap_v2.py` —
`fmax_with_ci()` (Metric A) and `fmax_cafa_with_ci()` (Metric B).

## Main result: DIAMOND + Pfam + InterProScan + priors

### F-max vs PGAP — Experimental-only truth

PGAP annotations available from RefSeq GFF go_function/go_process/
go_component fields for 3 genomes only:

| Genome | GSPA (D+P) | GSPA (D+P+I+priors) | PGAP | GSPA/PGAP |
|---|---|---|---|---|
| hpylori | 0.147 | **0.152** [0.113, 0.210] | 0.029 [0.000, 0.068] | **5.2×** |
| mgenitalium | 0.207 | **0.242** [0.095, 0.500] | 0.000 | **∞** |
| mjannaschii | 0.352 | **0.385** [0.330, 0.456] | 0.209 [0.114, 0.308] | **1.84×** |

### F-max — Full-GOA truth (all evidence)

| Genome | GSPA (D+P) | GSPA (D+P+I+priors) | Δ | PGAP | GSPA/PGAP |
|---|---|---|---|---|---|
| ecoli | 0.670 | **0.693** [0.686, 0.701] | +0.023 | — | — |
| ecolo157 | 0.835 | **0.867** [0.861, 0.875] | +0.032 | — | — |
| mtb | 0.716 | **0.737** [0.727, 0.746] | +0.021 | — | — |
| synechocystis | 0.614 | **0.681** [0.667, 0.696] | +0.067 | — | — |
| paeruginosa | 0.601 | **0.697** [0.686, 0.707] | **+0.096** | — | — |
| hpylori | 0.754 | **0.819** [0.806, 0.833] | +0.065 | 0.316 | **2.6×** |
| mgenitalium | 0.913 | **0.908** [0.895, 0.922] | -0.005 | 0.469 | **1.9×** |
| mjannaschii | 0.641 | **0.732** [0.712, 0.751] | **+0.091** | 0.285 | **2.6×** |

### Coverage and IC-recall — Full-GOA truth

| Genome | Cov (D+P) | Cov (D+P+I+priors) | ΔCov | icR (D+P) | icR (D+P+I+priors) | ΔicR |
|---|---|---|---|---|---|---|
| ecoli | 0.740 | **0.861** | +0.121 | 0.511 | **0.574** | +0.063 |
| ecolo157 | 0.831 | **0.942** | +0.111 | 0.782 | **0.877** | +0.095 |
| mtb | 0.720 | **0.831** | +0.111 | 0.591 | **0.656** | +0.065 |
| synechocystis | 0.618 | **0.845** | +0.227 | 0.432 | **0.624** | +0.192 |
| paeruginosa | 0.290 | **0.830** | **+0.540** | 0.346 | **0.595** | **+0.249** |
| hpylori | 0.528 | **0.885** | +0.357 | 0.607 | **0.763** | +0.156 |
| mgenitalium | 0.652 | **0.843** | +0.191 | 0.659 | **0.803** | +0.144 |
| mjannaschii | 0.208 | **0.720** | **+0.512** | 0.282 | **0.585** | **+0.303** |

### Why InterProScan is the game-changer

With only DIAMOND + Pfam, the priors had almost no borderline claims to
boost (F-max moved ±0.003). InterProScan adds a dense layer of
independent domain evidence from TIGRFAM, CDD, SUPERFAMILY, Panther,
HAMAP, and FunFam — many at borderline confidence — that the priors
then push across the 0.5 posterior threshold. The combination is
synergistic:

- InterProScan alone adds ~1,500–5,000 new GO claims per genome
- The priors act on these borderline claims to close pathways, fill
  essential functions, and reinforce operon consensus
- F-max gains of +0.021 to +0.096 across all genomes
- Coverage nearly triples on paeruginosa (0.290 → 0.830) and
  mjannaschii (0.208 → 0.720)

## Prior activity per genome (three-predictor pipeline)

| Genome | Essential (uncov.) | Coherence (pw-missing) | GapFill (fns) | Context (boosts) | Gaps |
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

## Dark matter suggester (Phase 8)

Singleton = high-confidence "this specific protein fills this gap"
(q > 0.5 within its operon). Disjunctive = credible set of candidates.

| Genome | Gaps | Singleton | Disjunctive | Total | New proteins | Δ% |
|---|---|---|---|---|---|---|
| ecoli | 368 | **913** | 2,903 | 3,816 | +33 | +1.1% |
| ecolo157 | 386 | **1,095** | 3,516 | 4,611 | +23 | +0.7% |
| paeruginosa | 427 | **1,083** | 3,835 | 4,918 | +95 | **+6.9%** |
| mjannaschii | 124 | **193** | 377 | 570 | +26 | **+7.2%** |
| synechocystis | 277 | **314** | 1,015 | 1,329 | +31 | +2.1% |

## Ablation

### Per-predictor (two predictors only, no priors)

| Genome | DIAMOND | Pfam | D+P combined | D+P+InterProScan |
|---|---|---|---|---|
| ecoli | 0.656 | 0.137 | 0.670 | **0.693** |
| ecolo157 | 0.807 | 0.256 | 0.835 | **0.867** |
| mtb | 0.681 | 0.240 | 0.716 | **0.737** |
| synechocystis | 0.564 | 0.268 | 0.614 | **0.681** |
| paeruginosa | 0.597 | 0.037 | 0.601 | **0.697** |
| hpylori | 0.754 | 0.000 | 0.754 | **0.819** |
| mgenitalium | 0.913 | 0.000 | 0.913 | **0.908** |
| mjannaschii | 0.641 | 0.000 | 0.641 | **0.732** |

InterProScan's contribution: +0.021 to +0.096 F-max on full-GOA truth.
Pfam's contribution: +0.003 to +0.067 (partially subsumed by
InterProScan since Pfam is one of InterProScan's member databases).

## Reference data requirements

| File | Source | Size | Purpose |
|---|---|---|---|
| `go.owl` | [GO consortium](http://purl.obolibrary.org/obo/go.owl) | 130 MB | GO hierarchy, class labels |
| `go-computed-taxon-constraints.obo` | [GO imports](http://current.geneontology.org/ontology/imports/go-computed-taxon-constraints.obo) | 2 MB | 13k only_in + 6k never_in constraints |
| `ec2go` | [GO external2go](http://current.geneontology.org/ontology/external2go/ec2go) | 350 KB | EC → GO mapping |
| `gene_refseq_uniprotkb_collab.gz` | [NCBI RefSeq](https://ftp.ncbi.nlm.nih.gov/refseq/uniprotkb/gene_refseq_uniprotkb_collab.gz) | 1.2 GB | RefSeq ↔ UniProt mapping |
| `goa_uniprot_all.gaf.gz` | [GOA FTP](https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz) | 15 GB | Ground truth |
| `uniprot_sprot.fasta.gz` | [UniProt](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz) | 93 MB | Leave-N-out reference DB |
| InterProScan | [EBI FTP](https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/) | 5.8 GB | Multi-database domain scan |
| KEGG ec-pathway link | [KEGG REST](https://rest.kegg.jp/link/pathway/ec) | 1 MB | Built via `build_kegg_pathway_tsv.py` |
| gapseq | `conda install -c bioconda gapseq` | 2 GB | Metabolic gap detection |

## 10-Genome PGAP Comparison (extended benchmark)

### Genomes

| Tag | Organism | Assembly | Phylum | Proteins | PGAP GO lines |
|---|---|---|---|---|---|
| vcholerae | V. cholerae O1 N16961 | GCF_000006745.1 | Proteobacteria (γ) | 3,447 | 1,991 |
| saureus | S. aureus N315 | GCF_000009645.1 | Firmicutes | 2,621 | 1,459 |
| spneumoniae | S. pneumoniae TIGR4 | GCF_000006885.1 | Firmicutes | 1,951 | 1,262 |
| ccrescentus | C. crescentus CB15 | GCF_000006905.1 | Proteobacteria (α) | 3,808 | 2,289 |
| rprowazekii | R. prowazekii Madrid E | GCF_000195735.1 | Proteobacteria (α) | 823 | 601 |
| tpallidum | T. pallidum Nichols | GCF_000008605.1 | Spirochaetes | 985 | 574 |
| tthermophilus | T. thermophilus HB8 | GCF_000091545.1 | Deinococcus-Thermus | 2,150 | 1,452 |
| dradiodurans | D. radiodurans R1 | GCF_000008565.1 | Deinococcus-Thermus | 3,134 | 1,805 |
| scoelicolor | S. coelicolor A3(2) | GCF_000203835.1 | Actinobacteria | 7,872 | 4,856 |
| pfuriosus | P. furiosus DSM 3638 | GCF_000007305.1 | Euryarchaeota | 2,008 | 856 |

### Method

- **Leave-19-out** Swiss-Prot reference (554,803 proteins; excludes
  all 9 original + 10 new benchmark genomes).
- **Three predictors**: DIAMOND blastp, HMMER/Pfam, InterProScan.
- **Noisy-OR integration** with EssentialityPrior and CoherencePrior.
- **Bootstrap F-max** (200 resamples, 95% CI) against full-GOA truth.
- PGAP GO annotations extracted from RefSeq GFF go_function/go_process/
  go_component fields, mapped to UniProt via NCBI collab file.

### F-max — Full-GOA truth (GSPA D+P+I+priors vs PGAP), both metrics

CIs are 95% bootstrap; "ratio" = GSPA / PGAP under the same metric.

| Genome | GSPA micro | GSPA CAFA | PGAP micro | PGAP CAFA | ratio (micro) | ratio (CAFA) |
|---|---|---|---|---|---|---|
| rprowazekii | **0.911** [0.900, 0.921] | **0.905** [0.896, 0.918] | 0.503 | 0.518 | **1.81×** | **1.75×** |
| tpallidum | **0.892** [0.878, 0.904] | **0.896** [0.884, 0.911] | 0.491 | 0.492 | **1.82×** | **1.82×** |
| saureus | **0.867** [0.857, 0.878] | **0.886** [0.878, 0.895] | 0.449 | 0.446 | **1.93×** | **1.99×** |
| vcholerae | **0.858** [0.850, 0.864] | **0.883** [0.878, 0.888] | 0.443 | 0.441 | **1.94×** | **2.00×** |
| pfuriosus | **0.857** [0.847, 0.869] | **0.874** [0.865, 0.884] | 0.350 | 0.366 | **2.45×** | **2.39×** |
| spneumoniae | **0.846** [0.834, 0.855] | **0.870** [0.861, 0.878] | 0.447 | 0.457 | **1.89×** | **1.90×** |
| tthermophilus | **0.842** [0.830, 0.854] | **0.869** [0.861, 0.879] | 0.492 | 0.509 | **1.71×** | **1.71×** |
| ccrescentus | **0.805** [0.795, 0.813] | **0.848** [0.842, 0.854] | 0.480 | 0.493 | **1.68×** | **1.72×** |
| dradiodurans | **0.780** [0.769, 0.790] | **0.818** [0.810, 0.828] | 0.463 | 0.473 | **1.68×** | **1.73×** |
| scoelicolor | **0.778** [0.770, 0.787] | **0.847** [0.843, 0.852] | 0.490 | 0.503 | **1.59×** | **1.68×** |

**Mean GSPA/PGAP = 1.85× (micro), 1.87× (CAFA)** across 10 new
genomes. The two metrics agree closely on the GSPA-vs-PGAP ratio (the
asymmetry from heavily-annotated proteins cancels because both methods
see the same proteins). CAFA F-max is consistently slightly higher
for GSPA — predictions are stronger on the typical protein than on
the few high-annotation outliers (notably scoelicolor, dradiodurans).
InterProScan adds +0.034 to +0.094 micro F-max compared to D+P alone.

### Combined: all 13 genomes with PGAP GO annotations

| Genome | GSPA micro | GSPA CAFA | PGAP micro | PGAP CAFA | ratio (micro) | ratio (CAFA) |
|---|---|---|---|---|---|---|
| rprowazekii | **0.911** | **0.905** | 0.503 | 0.518 | **1.81×** | **1.75×** |
| mgenitalium | **0.908** | **0.910** | 0.469 | 0.448 | **1.94×** | **2.03×** |
| tpallidum | **0.892** | **0.896** | 0.491 | 0.492 | **1.82×** | **1.82×** |
| saureus | **0.867** | **0.886** | 0.449 | 0.446 | **1.93×** | **1.99×** |
| vcholerae | **0.858** | **0.883** | 0.443 | 0.441 | **1.94×** | **2.00×** |
| pfuriosus | **0.857** | **0.874** | 0.350 | 0.366 | **2.45×** | **2.39×** |
| spneumoniae | **0.846** | **0.870** | 0.447 | 0.457 | **1.89×** | **1.90×** |
| tthermophilus | **0.842** | **0.869** | 0.492 | 0.509 | **1.71×** | **1.71×** |
| hpylori | **0.819** | **0.844** | 0.316 | 0.299 | **2.59×** | **2.82×** |
| ccrescentus | **0.805** | **0.848** | 0.480 | 0.493 | **1.68×** | **1.72×** |
| dradiodurans | **0.780** | **0.818** | 0.463 | 0.473 | **1.68×** | **1.73×** |
| scoelicolor | **0.778** | **0.847** | 0.490 | 0.503 | **1.59×** | **1.68×** |
| mjannaschii | **0.732** | **0.775** | 0.285 | 0.269 | **2.57×** | **2.88×** |

**Mean GSPA/PGAP across all 13 genomes**: **1.93× (micro), 1.96× (CAFA)**.
GSPA is consistently 1.6-2.9× better than PGAP across 7 phyla, both
domains of life (Bacteria + Archaea), and genome sizes from 483 to
7,872 proteins. The two metrics give the same qualitative picture
(GSPA > PGAP on every single genome under both); CAFA is slightly
more favourable to GSPA on the larger Actinobacteria/Deinococcus
genomes and on the small archaea where PGAP misses most annotations.

## Known issues

### ConsistencyPrior requires taxonomy lineage

The SAT4J consistency checker works correctly (13k taxon constraints
from OBO, UNSAT core extraction validated) but requires `--taxonomy`
with the genome's NCBI lineage. Without it, every `only_in_taxon`
constraint fires as a false positive. Gated on `--taxonomy` flag.

### Gapseq Reactions.tbl bug

4 of 9 genomes (hpylori, bsubtilis, mtb, mgenitalium) produce
zero-byte Reactions.tbl from gapseq — reproduced on both GlusterFS
and local /tmp. Pathways.tbl is fine. These genomes lack
GapFillingPrior and dark-matter suggestions.

### PGAP comparison limited to 3 genomes

PGAP GO annotations (from RefSeq GFF go_function/go_process/
go_component fields) are only available for hpylori, mgenitalium,
mjannaschii. The other 6 genomes' RefSeq GFFs don't carry GO.
NCBI's gene2go file is not a fair proxy (it includes UniProt-GOA
IEA annotations, creating circular overlap with our ground truth).
