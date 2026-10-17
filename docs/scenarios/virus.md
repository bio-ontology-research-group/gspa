# Viral genomes

Use this workflow for phage and other viral genomes, either isolated or as viral
contigs recovered from a metagenome. Viral proteomes are small and divergent, so
most of the value comes from per-protein evidence (structure and domain search),
while several genome-scale steps designed for cellular organisms must be adapted
or left out.

## Inputs

| Input | Notes |
|---|---|
| Viral genome FASTA | one sequence per genome or contig |
| Proteins + coordinates | Prodigal/Pyrodigal, Pharokka, or `gspa annotate -k virus` |
| `claims.jsonl` | DIAMOND/MMseqs2 against Swiss-Prot, Pfam/InterProScan, FoldSeek; viral protein families are better covered by structure and profile search than by Swiss-Prot similarity |
| Operons TSV | optional; phage genes are densely packed, so distance-based operons span long modules |

## Step 1: call genes and run predictors

```bash
gspa annotate --genome phage.fna -o annot/ -k virus --db uniprot_sprot.dmnd -t 8
```

Real output (`test-data/genomes/lambda_phage.fna`): Prodigal called 61 genes,
DIAMOND assigned 132 product annotations to 55 of them, and the operon
detector found 8 operons covering 52 genes:

```text
lambda_phage_00001  CUSTOM  Terminase small subunit             1.0000  diamond  ISS  predictor:diamond@1.000
lambda_phage_00001  CUSTOM  DNA-packaging protein NU1 homolog   0.9830  diamond  ISS  predictor:diamond@0.983
```

The `annotate` step runs Prodigal with translation table 11. That is correct for
most phages of bacteria; for viruses with other genetic codes call genes
yourself and pass `--proteins` with `--gff3`.

Build `claims.jsonl` from the predictor outputs as described in
[the bacterial workflow](bacterial-isolate.md#step-2-run-predictors-and-build-claimsjsonl).
The `annotate` operon file `annot/operons_for_integrate.tsv` can be passed to
`--operons` directly.

## Step 2: integrate

```bash
gspa integrate \
  --claims claims.jsonl --out integrated.tsv \
  --provenance provenance.json \
  --go-owl go.owl \
  --essential-profile virus \
  --operons annot/operons_for_integrate.tsv \
  --enable-priors essentiality,coherence,genomic_context \
  --lite
```

Options that matter:

- `--essential-profile virus` loads the bundled viral hallmark list (24 GO
  terms: viral process, viral life cycle, replication, virion assembly,
  packaging, entry and release). The help text lists only bacteria, archaea and
  eukaryote, but the `virus` preset exists and loads
  (`Essential functions (preset=virus): 24 GO terms`).
- `--enable-priors` without `consistency`, `anatomical` and `gap_filling`:
  phages have no metabolic network of their own and the taxon constraints are
  written for cellular organisms (see below).
- `--taxon virus` resolves to `NCBITaxon_10239`. Only add it, together with
  `--enforce-stage2`, if you want terms that require a cellular organism to be
  removed.

!!! warning "Taxon enforcement on viral proteins"
    GO taxon constraints restrict many terms to cellular organisms. In a test
    run where cellular (*M. genitalium*) claims were integrated with
    `--taxon virus --enforce-stage2`, 543 annotations were acted on, the most
    frequent being *cytosol*, *cytoplasm* and *translation*. For a phage this is
    usually what you want for location terms, but host-derived genes (auxiliary
    metabolic genes, tRNA-modifying enzymes, translation factors in giant
    viruses) lose their annotations too. Review `--enforcement-actions-out`, or
    use `--enforce-consistency-mode flag` to mark instead of remove.

## Step 3: evaluate

```bash
uv run python benchmark/make_synth_gff_gaf.py --fasta proteins.faa \
  --tsv integrated.tsv --gff-out synth.gff --gaf-out synth.gaf
gspa evaluate -i synth.gff -a synth.gaf --go-owl go.owl -k virus --lite -o quality.json
```

Completeness is scored against the viral hallmark list when `-k virus` is
given. Pathway coherence is not meaningful for most viruses.

## Prophages and viral contigs inside bacterial genomes

The Nextflow pipeline runs geNomad, CheckV and PhiSpy as a viral/prophage
region track (`--run_genomad`, `--run_checkv`, `--run_phispy`; see
[Nextflow and HPC](nextflow-hpc.md) and [Nextflow parameters](../reference/nextflow-params.md)).
Region calls are evidence about the host genome; annotate the prophage
proteins as part of the host with the bacterial profile, or extract the region
and annotate it here.

## Runtime

Bacteriophage T4 (278 proteins) integrated in 21 s with 2.3 GB peak RSS
(`results/runtime_memory_empirical_scaling.tsv`); nearly all of it is the fixed
cost of loading constraint tables.

## Pitfalls

- Using the bacterial essential profile: every phage then looks 0 % complete
  and the essentiality prior boosts cellular functions.
- Expecting Swiss-Prot similarity to cover a novel phage; most proteins will
  have no claim at all. Add FoldSeek and profile-HMM evidence.
- Combining `--taxon virus` with `--enforce-stage2` without checking the
  actions log (see the warning above).
