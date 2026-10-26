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
policy. The first audit (job 9743) rejected preliminary EC identifiers such as
`1.1.1.n11`. These are valid preliminary identifiers in the
[ENZYME nomenclature](https://enzyme.expasy.org/enzyme.html); the corrected audit
counts them separately from complete numeric ECs and hierarchical classes.
Job 9744 completed all 132 matrices: no invalid score or matrix-shape failures,
89 files with fully resolved identifier joins and 43 with some unresolved IDs.
The [full audit](../../benchmark/v3/results/development-88a15e7/price-matrix-audit.json)
records 101 distinct absent identifiers (80 versioned, 21 versionless). None is a
Price target protein. They remain candidate records with unavailable local context.

The 80 exact-version proteins were fetched from NCBI. None has an identical
sequence in its newly acquired assembly protein FASTA. A further diagnostic aligns
these proteins against exact-version assembly DNA, retaining only full-length,
ungapped, identical TBLASTN alignments: 51 have a single returned exact alignment.
This establishes candidate coordinate repairs, not an exhaustive uniqueness proof;
the original layouts and family inputs have not yet been replaced with those
repairs. `recover_price_coordinates.py` reproduces the diagnostic; its
[report](../../benchmark/v3/results/development-88a15e7/price-coordinate-recovery.json)
records commands, hashes and tool version. Versionless identifiers are not assigned
an arbitrary accession version.

## Cohort families and metadata

`benchmark/v3/prepare_price_cohort.py` turns the layout directory, an MMseqs2
easy-cluster membership TSV (id 0.3, coverage 0.8 both directions, over the namespaced
combined FASTA) and the frozen NCBI taxonomy EFETCH XML into the two strict production
inputs: `cohort_orthogroups.tsv` (`genome_id, protein_id, orthogroup_id, subfamily_id,
weight`, read by `CohortOrthogroupReader.readCohortFile`) and `cohort_metadata.json`
(`CohortMetadata`). It also writes `families.tsv` and a `report.json` with the SHA-256
of every input, the MMseqs2 version and command (from `--mmseqs-log`, or
`--mmseqs-version` plus `--mmseqs-command`; one is required), and the coverage audit.

```bash
python benchmark/v3/prepare_price_cohort.py \
  --layout /tmp/gspa-v3-price-layout-v3 \
  --clusters /tmp/gspa-v3-price-families/id30_cov80_cluster.tsv \
  --taxonomy /tmp/gspa-v3-price-families/taxonomy.xml \
  --mmseqs-log /tmp/gspa-v3-price-families/id30_cov80-retry.log \
  --out <fresh directory>
```

What the outputs mean:

- **Families are sequence-similarity clusters**, not established orthologs. Paralogs
  of one genome can share a family, and the same protein accession in two genomes is
  two proteins (the combined FASTA namespace `<assembly>__<protein>` is split on the
  first `__` only). Family IDs are `fam30_` plus a hash of the representative, with no
  colon, because the production reader strips an `ns:` prefix.
- **`familySearchComplete` is true only when every acquired protein occurs exactly
  once** in the membership and no member ID is foreign (unknown genome, unacquired
  protein, stripped namespace, other protein version). Repeated or foreign IDs, and
  any missing protein, are rejected; `--allow-partial` writes an uncertified output
  with `familySearchComplete=false` and the measured coverage. Absence of a family is
  detection-limited to this cluster protocol and is not evidence that remote
  homologs are absent.
- **`clade` is the NCBI genus taxid and `strainCluster` the NCBI species taxid**, read
  from the frozen lineage. The species block is a conservative redundancy proxy, not
  an ANI cluster, and strains or species are not independent units. A missing taxon or
  rank leaves the field absent (unknown) and is listed in the report. A top-level
  `TaxId` that was not requested, including a merged ID, is an error.
- **Subfamilies** (optional 70% subclusters, `--subclusters`) must nest in one 30%
  family. A violation fails by default; `--subfamily-policy combined` uses the
  family+subcluster signature for every protein and records a warning. Without a
  subcluster file `subfamily_id` is `-`.

Only identifiers, assembly status and contig lengths enter the outputs; functional GFF
annotations never do. Tests: `benchmark/v3/tests/test_price_cohort.py`.

Status: all 115,542 acquired proteins occur exactly once in 35,298 sequence-similarity
clusters; 14,255 clusters occur in multiple genomes and 4,269 contain within-genome
paralogs. All 18 distinct requested taxids resolved. The
[cohort report](../../benchmark/v3/results/development-88a15e7/price-cohort-report.json)
binds the inputs and outputs by hash. The first clustering attempt exhausted `/tmp`;
its failure log is preserved and the successful retry used disk-backed scratch.
These are context inputs, not evidence of prediction accuracy.

The production metadata now explicitly declares the 279 acquired proteins without
one usable coordinate as `unplacedProteins`. A declared family member remains
present with unknown neighborhood context; it is never counted as absent because
its coordinate is missing. A placed unsupported paralog plus an unresolved paralog
also has unknown context. Arbitrary undeclared protein IDs still fail validation.

## Completed 7-Variant × 6-Arm Production Execution

Using the verified layouts (`115,263` placed proteins + `279` `unplacedProteins`), MMseqs2 30%-identity/80%-coverage families (`35,298` clusters, `familySearchComplete: true`), and NCBI genus/species metadata (`18` taxa), all 7 predictor variants (`clean`, `dpz`, `enzbert`, `graphec`, `mapred`, `topec`, `ensemble_pool4`) were run end-to-end through `gspa integrate` across all 6 inference arms (`stage0`, `stage1`, `stage1_rlc`, `stage1_rlc_cohort`, `learned_logistic`, `learned_gbdt`) on Unimatrix (`benchmark/v3/results/price-context-ablation/`). Full results, SHA-256 digests, and locus-decoy metrics are reported in [`docs/benchmarks/v3-price-development.md`](v3-price-development.md).
