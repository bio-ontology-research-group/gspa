# Archaeal genome

The archaeal workflow is the [bacterial isolate workflow](bacterial-isolate.md)
with three changes: the essential-function profile, the asserted taxon, and the
handling of steps whose reference data are bacterial. Read that page first; this
page lists only what differs.

## When to use it

You have an assembled archaeal genome (for example a methanogen or a
hyperthermophile) or a high-quality archaeal MAG. For incomplete MAGs also read
[Metagenome-assembled genomes](mag-metagenome.md).

## Inputs

Same as for bacteria. Gene calling with Prodigal translation table 11 (what
`gspa annotate` and Prokka `--kingdom Archaea` use) is appropriate for archaea.
`gspa annotate -k archaea` loads `config/archaea.yaml`, which differs from the
bacterial config only in the kingdom and the completeness profile.

```bash
gspa annotate --genome genome.fna -o annot/ -k archaea --db uniprot_sprot.dmnd -t 16
```

Build `claims.jsonl` and operons exactly as in
[steps 2 and 3 of the bacterial page](bacterial-isolate.md#step-2-run-predictors-and-build-claimsjsonl).
Swiss-Prot and GOA are sparser for archaea, so structure (FoldSeek) and
domain (InterProScan, Pfam) evidence carry more weight than for a model
bacterium; include them if you can.

## Integrate

```bash
gspa integrate \
  --claims claims.jsonl --out integrated.tsv \
  --provenance provenance.json \
  --go-owl go.owl \
  --taxon archaea --essential-profile archaea \
  --operons operons.tsv \
  --pathways kegg_pathways.tsv --ec2go ec2go \
  --enforce-stage2 --enforce-consistency-mode relax \
  --enforcement-actions-out stage2_actions.tsv \
  --lite
```

What changes and why:

- `--essential-profile archaea` loads the bundled archaeal core-function list.
  With the bacterial profile, `EssentialityPrior` would push bacterial-only
  processes onto archaeal proteins and completeness would be scored against the
  wrong list.
- `--taxon archaea` (resolved to `NCBITaxon_2157`) makes the SAT checker treat
  bacteria-only and eukaryote-only GO terms as inconsistent. Use a specific
  `NCBITaxon_<id>` from your lineage call if you have one; the bundled taxon
  hierarchy decides which constraints apply.
- `--enforce-consistency-mode relax` keeps the most specific valid ancestor of
  a clade-invalid term instead of deleting it. It needs `--go-owl`.

## Steps to leave out or adapt

- **Flux-certified refinement**: `--fba-gram` only offers Gram-negative and
  Gram-positive biomass templates, which are bacterial. Do not use
  `--enable-fba-refiner` on archaea unless you supply your own data directory.
- **Metabolic gaps**: gapseq's MetaCyc pathway predictions are calibrated on
  bacteria. Gap lists from it can be used, but inspect `--suggestions-out`
  before promoting anything (see [Metabolic gaps and dark matter](metabolic-gaps-dark-matter.md)).

## Evaluate

```bash
uv run python benchmark/make_synth_gff_gaf.py --fasta proteins.faa \
  --tsv integrated.tsv --gff-out synth.gff --gaf-out synth.gaf
gspa evaluate -i synth.gff -a synth.gaf --go-owl go.owl \
  -k archaea --taxon archaea --lite -o quality.json
```

`-k archaea` selects the archaeal essential-function list for completeness.

## Runtime

In `results/runtime_memory_empirical_scaling.tsv`, *Methanocaldococcus
jannaschii* (1,715 proteins) integrated in 56 s with 3.0 GB peak RSS.

## Pitfalls

- Passing `--essential-profile bacteria` by habit. Check the log line
  `Essential functions (preset=archaea)`.
- Using `-k bacteria` (the default of `gspa compare`) when scoring archaeal
  annotations.
- Expecting eggNOG/KEGG pathway coverage comparable to bacteria; low pathway
  coherence on an archaeon often reflects missing reference pathways rather
  than missing genes.
