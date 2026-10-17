# Eukaryotic genomes

Use this workflow for fungi, protists, plants and animals. Compared with
prokaryotes, three things change: gene models must come from a eukaryotic gene
predictor, operon context does not apply, and anatomical constraints become
relevant for multicellular organisms (a plant protein cannot take part in
*neural crest cell migration*, a moss cannot develop a *camera-type eye*).

## Inputs

| Input | How to produce it |
|---|---|
| Genome FASTA + GFF3 gene models | BRAKER/AUGUSTUS, Funannotate, or an existing reference annotation |
| Protein FASTA | from the same gene models (one isoform per gene is usually enough) |
| `claims.jsonl` | predictors on the protein FASTA, converted as in [the bacterial workflow](bacterial-isolate.md#step-2-run-predictors-and-build-claimsjsonl) |
| `go.owl` | required for ancestor relaxation |
| Organism taxon | an `NCBITaxon_<id>` at species or a higher rank |

GSPA's built-in gene callers (Prodigal, Pyrodigal) are prokaryotic. Do not let
`gspa annotate` call genes on a eukaryotic genome; give it the gene models:

```bash
gspa annotate --genome genome.fna --gff3 genes.gff3 -k eukaryote \
  --go-owl go.obo -o annot/
# or: --proteins proteins.faa --gff3 genes.gff3
```

CDS rows that share an `ID` are concatenated in transcription order before
translation, so multi-exon genes translate correctly. With `--proteins`,
proteins are joined to contigs through the CDS attributes `protein_id`, `ID`
or `locus_tag`. `-k eukaryote` loads `config/eukaryote.yaml`, which disables
operon prediction.

## Integrate

```bash
gspa integrate \
  --claims claims.jsonl --out integrated.tsv \
  --provenance provenance.json \
  --go-owl go.owl \
  --taxon NCBITaxon_7227 --essential-profile eukaryote \
  --enable-priors essentiality,coherence,consistency,anatomical \
  --enforce-stage2 --enforce-consistency-mode relax \
  --enforcement-actions-out stage2_actions.tsv \
  --lite
```

(`NCBITaxon_7227` is *Drosophila melanogaster*, as used in the benchmark
tables; substitute your organism.)

Options that matter:

- `--enable-priors ...,anatomical` turns on the GO-Uberon/PO/CL constraints.
  They do two things: GO terms whose anatomical structure cannot exist in the
  lineage become inconsistent, and borderline developmental terms are boosted
  when markers of the target structure are supported in the genome. Omitting
  `anatomical` from the list also removes the anatomical constraints from the
  Stage 2 enforcement step.
- `--enable-priors` without `genomic_context`: there are no operons. Leave out
  `gap_filling` unless you supply gaps.
- `--enforce-consistency-mode` decides what happens to an inconsistent term:

  | Mode | Effect |
  |---|---|
  | `relax` (default) | replace the term with its most specific valid ancestor, recovering the unpenalised evidence |
  | `remove` | delete the term |
  | `downrank` | lower its score |
  | `flag` | keep it and record the violation |
  | `minimal-flip` | weighted MaxSAT: demote the minimum-cost set of terms that makes the protein consistent |

### Relax or remove

`results/pillar3_relax_vs_remove.tsv` compares the two modes on six eukaryotes
(CAFA-style F-max). Both modes raise F-max well above uncorrected Noisy-OR in
every species. `relax` was higher in four species (*C. elegans*,
*D. melanogaster*, *A. thaliana*, *M. musculus*), and `remove` was higher in
two (*S. cerevisiae*: 0.9604 vs 0.9592; *P. patens*: 0.9828 vs 0.9792). The
pattern is the expected one: `remove` gives the highest precision, `relax`
keeps more recall by retaining a general but correct ancestor. Use `relax` when
you want coverage, `remove` when you want a conservative set.

### What the enforcement log looks like

Real output (anatomical relaxation from the test run on the bundled data):

```text
dimension   action protein_id     term       term_label        reason
consistency relax  NC_000908.2_25 GO:0048468 cell development  anatomical-inconsistent for NCBITaxon_2: GO:0014045 (establishment of endothelial blood-brain barrier) [results in development of UBERON:0001017 (central nervous system)] violates only_in_taxon(NCBITaxon_33208) [UBERON(via UBERON:0000483)] -> relaxed to GO:0048468
```

The reason names the original term, the anatomical structure it refers to, and
the taxon constraint on that structure. Relaxed terms appear in
`integrated.tsv` with `stage2_relax:1.000` in `priors_fired`.

## Evaluate

```bash
uv run python benchmark/make_synth_gff_gaf.py --fasta proteins.faa \
  --tsv integrated.tsv --gff-out synth.gff --gaf-out synth.gaf
gspa evaluate -i synth.gff -a synth.gaf --go-owl go.owl \
  -k eukaryote --taxon NCBITaxon_7227 --lite -o quality.json
```

`make_synth_gff_gaf.py` cleans UniProt-style headers (`sp|ACC|NAME` becomes
`ACC`); use `--acc-mapper` if your protein IDs need another mapping.

## Runtime

From `results/runtime_memory_empirical_scaling.tsv`: *S. cerevisiae* (6,060
proteins) 52 s and 2.5 GB peak RSS, *D. melanogaster* (13,968) 63 s and 3.0 GB,
*M. musculus* (21,986) 51 s and 2.7 GB. Runtime is dominated by loading the
ontology and constraint tables, not by proteome size.

## Pitfalls

- **Isoforms**: claims for several isoforms of a gene are integrated as separate
  proteins. Collapse to one isoform per gene, or expect duplicated hypotheses.
- **Too coarse a taxon**: `--taxon eukaryote` cannot exclude animal-only terms
  from a plant. Use the most specific taxon in the bundled hierarchy.
- **Missing `anatomical`** in a custom `--enable-priors` list silently disables
  anatomical enforcement in Stage 2 as well.
- **No `--go-owl`**: `relax` cannot find ancestors and removes the term.
