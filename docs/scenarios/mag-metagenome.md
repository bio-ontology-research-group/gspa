# Metagenome-assembled genomes

Use this workflow for MAGs (binned genomes from metagenomes) and for any
multi-contig input where contigs may come from different replicons or
organisms. Two things differ from an isolate: completeness expectations must
allow for missing genome fractions, and quality must be computed per contig or
per bin rather than pooled over a whole assembly.

## Inputs

| Input | How to produce it |
|---|---|
| Bin FASTA (one file per MAG) or a contig set | your binning pipeline |
| CheckM2 completeness / contamination | `checkm2 predict` (also available in the Nextflow pipeline) |
| Lineage | GTDB-Tk `classify_wf`; needed for `--taxon` |
| Proteins + GFF3 | Prodigal/Pyrodigal in metagenome mode, Prokka `--metagenome`, or `gspa annotate` |
| `claims.jsonl` per bin | as in [the bacterial workflow](bacterial-isolate.md#step-2-run-predictors-and-build-claimsjsonl) |
| Operons per bin | `benchmark/make_operons.py` on the bin GFF (operons never span contigs) |

## Option A: annotate a contig set with per-contig metrics

`gspa annotate` can take a metagenome or a multi-replicon assembly plus a GFF3
and compute genome-scale metrics for each contig separately:

```bash
gspa annotate \
  --genome contigs.fna --gff3 contigs.gff3 \
  --mag -k auto \
  --go-owl go.obo --metrics-scope contig \
  --base-predictor light \
  --neural-sidecar $PWD/benchmark/neural/run_neural_predictors.py \
  --deepgo-plusplus-light-assets /path/to/dgpp_assets \
  --min-score 0.3 \
  -o mag_out/
```

- `--genome` + `--gff3`: CDS are translated against their contig and keep the
  contig assignment. Without a GFF3, genes are called.
- `--mag` switches the input type to MAG (lowered completeness expectations
  in the quality report).
- `--metrics-scope contig` (the default once `--go-owl` is given) treats each
  contig as its own genome `<genome>:<contig>`. Use `genome` for pooled
  metrics or `both`.
- `--min-score` drops low-scoring GO annotations right after prediction;
  learned predictors propagate thousands of near-zero terms per protein.

Outputs: `<genome>_annotations.tsv`, `<genome>_quality_per_contig.tsv` (one
row per contig), one `<genome>_<contig>_quality.json` per contig, the operon
files and `<genome>_timing.tsv`. Add `--infer-taxon` to write
`<genome>_taxon_inference.tsv` (candidate taxa with forbidden/support term
counts and the chosen lineage), or `--taxon` plus `--enforce-consistency` to
remove lineage-inconsistent terms.

The validated run in `GENOME_GFF3_ANNOTATION.md` (mock community of two
chromosomes, two plasmids and phage lambda) shows why per-contig scoring
matters: chromosomes score near-complete while plasmids and the phage score far
lower, a distinction that pooling would hide. Base predictor setup:
[DeepGO-PlusPlus base predictor](deepgo-plusplus.md).

## Option B: integrate each bin

Treat each MAG as one genome and give it its own lineage:

```bash
gspa integrate \
  --claims bins/bin_017.claims.jsonl --out bins/bin_017.integrated.tsv \
  --provenance bins/bin_017.provenance.json \
  --go-owl go.owl \
  --taxon NCBITaxon_<from GTDB-Tk/NCBI> --essential-profile bacteria \
  --operons bins/bin_017.operons.tsv \
  --enforce-stage2 --enforce-consistency-mode relax \
  --enforcement-actions-out bins/bin_017.stage2.tsv \
  --lite
```

For hundreds or thousands of bins, use one JVM with
`--batch-manifest` (shared ontology and constraint tables, per-row taxon,
profile and operons); see [Many genomes in batch](batch-multi-genome.md).

!!! note "Taxon inference in `integrate`"
    `gspa integrate --infer-taxon` loads the SAT checker but, in 2.1.0, did not
    assign a taxon in our tests (the log reports
    `SAT checker: ... (organismTaxon=none)`). Pass the GTDB-Tk/NCBI lineage
    with `--taxon` instead, or infer it with `gspa annotate --infer-taxon` and
    copy the chosen taxon from `<genome>_taxon_inference.tsv`.

## Completeness on incomplete genomes

The essentiality prior boosts candidates for essential functions that have no
annotation. On a 70 %-complete MAG some essential genes are absent from the
assembly, and the prior can only raise weak candidates that exist; it never
creates claims. Interpret completeness together with CheckM2 completeness,
and use `gspa evaluate --mag` when scoring:

```bash
gspa evaluate -i bin_017.synth.gff -a bin_017.synth.gaf --go-owl go.owl \
  -k bacteria --mag --lite -o bin_017.quality.json
```

## Dark matter in MAGs

MAG proteomes have a large fraction of proteins without any claim. The gap-driven
suggesters in [Metabolic gaps and dark matter](metabolic-gaps-dark-matter.md)
rank such proteins as candidates for missing pathway steps. Their per-phylum
effect on 97 MAGs is summarised in `results/pillar3_dark_matter_97mags.tsv`
(candidate fills per MAG and pathway completeness before and after).

## Pitfalls

- **Contamination**: a contaminated bin mixes lineages, and `--taxon` then
  removes correct annotations of the contaminant fraction. Filter bins by
  CheckM2 contamination first.
- **Prefixed protein IDs**: binning pipelines often rename contigs; make sure
  the IDs in claims, operons and GFF agree after renaming.
- **Pooled metrics on an assembly** (`--metrics-scope genome`) mix replicons and
  organisms; keep `contig` for metagenomes.
