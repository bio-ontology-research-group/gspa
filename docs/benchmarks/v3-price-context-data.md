# Price whole-genome context inputs

This acquisition extends the reused Price development set; it does not create a
new blind test. The earlier compact artifact covered target proteins and limited
neighbors. All 132 full-proteome predictor matrices (six predictors, 22 assemblies;
23.6 GB) have now been copied from the recorded IBEX sources to
`/data/hohndor/gspa-v3/price-raw-sources/<predictor>/<assembly>/`.

Slurm job 9740 verified every destination's size and SHA-256 against the original
retrieval manifest. [The full verification report](../../benchmark/v3/results/development-3c75a38/price-source-verification.json)
retains all paths and digests. This establishes transfer integrity only.
After installing the required pandas dependency, job 9742 successfully probed one
matrix per predictor. All six are protein-by-EC DataFrames, and every protein ID
in each probed matrix matched its assembly's acquired FASTA exactly. The
[six schema probes](../../benchmark/v3/results/development-3c75a38/price-schema-probes.json)
record shapes, sample axes and overlap counts. This is not yet validation of all
132 matrices, sequence identity to the original predictor inputs, or genomic
coordinate joins; those checks remain required before context inference.

For every exact assembly accession/version, NCBI Datasets supplied protein FASTA,
GFF3, and assembly/annotation metadata. The returned accession matched the request,
and both sequence and GFF files matched the package MD5 checksums in all 22 cases.
[The acquisition log](../../benchmark/v3/results/development-3c75a38/price-genomes-acquisition.jsonl)
records URLs, archive SHA-256, retrieval time, organism, and annotation provenance.
Local packages are under `/tmp/gspa-v3-price-genomes`; a copy is under the remote
source directory's `genomes/` subdirectory.

The [NCBI Datasets documentation](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genomes/download-genome/)
describes these data packages. The equivalent CLI request for each accession is:

```bash
datasets download genome accession GCA_000006965.1 \
  --include protein,gff3 --filename GCA_000006965.1.zip
```

The recorded REST requests use the exact version, with `GENOME_GFF` and
`PROT_FASTA` included. No current-version substitution is permitted. Re-downloads
may have different archive bytes; validate identifiers and member checksums and
record the new archive provenance rather than silently replacing the frozen copy.

GFF product names, EC assignments, and other functional annotations must not
become context evidence: they can contain the functions being evaluated. Use only
coordinates, strand, contig, identifiers, and independently generated prediction
inputs. Verify sequence identity when joining predictor protein IDs to the
acquired annotations. Missing or ambiguous mappings stay explicit. Family
clustering, phylogenetic metadata, reaction mappings and the full context benchmark
remain pending; acquiring these files is not evidence of improved prediction.


## Coordinate preparation

`benchmark/v3/prepare_price_layout.py` now prepares all 22 exact-version packages:
115,542 protein sequences, of which 115,263 have one unambiguous CDS locus.
The remaining 279 retain their sequences but have unavailable local coordinates
(missing, split or multiple loci); they are not silently assigned a merged locus.
The [layout report](../../benchmark/v3/results/development-3c75a38/price-layout-report.json)
records each unavailable protein, source hashes and output digests. FASTA headers
in the combined sequence file contain only assembly and protein IDs. GFF product
names and functional annotations never enter layouts, sequences or metadata.

The output includes production-format layout TSVs, per-protein sequence digests,
assembly/contig metadata and a namespaced combined FASTA for family clustering.
Family-search completeness and phylogenetic independence are deliberately absent
until measured. These sequences describe the acquired NCBI package; the original
predictor FASTAs have not been supplied for direct sequence comparison.

`audit_price_matrices.py` checks all 132 matrix hashes again before loading them,
validates score ranges and unique EC/protein axes, and audits exact or uniquely
versionless accession joins. An explicit different protein version is never
silently joined. Its counts of nonzero and threshold-retained candidates are
resource-planning diagnostics; they are not benchmark results or a chosen pruning
policy. Slurm job 9743 runs this complete audit.
