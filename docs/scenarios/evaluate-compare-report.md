# Evaluate, compare, report

Use this workflow to measure the genome-scale quality of an annotation set,
compare multiple pipelines on the same genome, and produce publication-ready
HTML reports.

Traditional protein function assessment focuses on per-protein accuracy against
holdouts. In contrast, the Genome Annotation Evaluation Framework (GAEF) assesses
whether the annotated proteome functions as a coherent, biologically viable
organism. GAEF evaluates five orthogonal dimensions:

1. **Completeness**: coverage of indispensable, universal domain-essential functions
   (such as DNA replication, ribosome biogenesis, and core tRNA synthetases).
2. **Process coherence**: topological connectivity across Gene Ontology biological
   processes via `has_part` axioms evaluated in an ELK reasoner.
3. **Complex coherence**: structural pairing of obligate heteromeric complex subunits,
   penalizing orphan singletons.
4. **Consistency**: adherence to taxonomic and anatomical true-path rules, verifying
   zero impossible terms are assigned.
5. **Information content (mean IC)**: specificity of assigned terms, measuring whether
   annotations provide granular molecular insights rather than high-level root categories.

The toolset comprises three commands:

- `gspa evaluate`: computes GAEF metrics for a single GFF3 and GAF pair, writing a JSON report.
- `gspa compare`: benchmarks multiple GAF annotation sets against the same genome GFF3,
  producing a side-by-side comparative HTML report.
- `gspa report`: renders existing JSON evaluation reports into standalone HTML dashboards.

## What you need

| Input | Required | Description |
|---|---|---|
| Genome GFF3 (`genome.gff3`) | yes | Structural annotations with CDS coordinates and IDs |
| Functional annotations (`annotations.gaf`) | yes | Standard Gene Association File (GAF 2.0/2.2) format |
| `go.owl` | yes | Gene Ontology release file |
| Organism kingdom / domain | recommended | `bacteria`, `archaea`, `eukaryote`, or `virus` |
| Pathway TSV + `ec2go` | optional | KEGG or MetaCyc pathway definitions and EC-to-GO mapping |
| Taxonomic constraint tables | optional | Bundled NCBI taxonomy and taxon-constraint files |

## Step 1: convert GSPA outputs to GFF3 and GAF

If you integrated annotations using `gspa integrate`, the output is an
`integrated.tsv` file. Convert this file and your protein FASTA into GAEF-compatible
GFF3 and GAF files using `benchmark/make_synth_gff_gaf.py`:

```bash
uv run python benchmark/make_synth_gff_gaf.py \
  --fasta proteins.faa \
  --tsv integrated.tsv \
  --threshold 0.50 \
  --gff-out synth.gff \
  --gaf-out synth.gaf
```

- `--threshold`: posterior probability cutoff for accepting an annotation (default: `0.50`).
- `--gff-out`: writes a synthetic GFF3 where each protein is placed as a CDS feature
  matching the GAF target identifier.
- `--gaf-out`: writes standard 17-column GAF lines with GO aspect codes (F, P, C).

## Step 2: evaluate an annotation set with gspa evaluate

Run `gspa evaluate` on the generated files:

```bash
gspa evaluate \
  -i synth.gff \
  -a synth.gaf \
  --go-owl reference_data/go.owl \
  --kingdom bacteria \
  --taxon NCBITaxon_562 \
  --ec2go reference_data/ec2go.txt \
  --pathways reference_data/kegg_pathways.tsv \
  -o quality_report.json \
  --lite
```

### Options

- `-i, --input`: genome structure GFF3 file.
- `-a, --annotations`: functional annotations in GAF format.
- `--go-owl`: path to Gene Ontology OWL file.
- `-k, --kingdom`: organism domain for essential function profiles (`bacteria`, `archaea`, `eukaryote`, `virus`).
- `--taxon`: organism taxon identifier (e.g. `NCBITaxon_562`).
- `--lite`: skips ELK reasoning, accelerating runtime while bypassing process coherence.
- `-o, --output`: output JSON path (default: `quality_report.json`).

### Understanding the JSON report

The output JSON file contains comprehensive evaluation metrics:

```json
{
  "genome_id": "synth",
  "assessment_date": "2026-10-17T09:00:00Z",
  "summary": {
    "total_proteins": 4363,
    "annotated_proteins": 3812,
    "annotation_rate": 0.8737
  },
  "completeness": {
    "score": 0.9583,
    "present_functions": ["GO:0006260", "GO:0006412", "GO:0003887"],
    "missing_functions": ["GO:0006313"]
  },
  "coherence": {
    "process_coherence": 0.7420,
    "pathway_coherence": 0.8150,
    "complex_coherence": 0.8920
  },
  "consistency": {
    "consistent": true,
    "violation_count": 0,
    "violations": []
  },
  "information_content": {
    "mean_ic": 8.4210,
    "median_ic": 8.9150
  }
}
```

## Step 3: compare multiple annotation pipelines with gspa compare

To benchmark different tools (such as Prokka, Bakta, PGAP, and GSPA) on the
same genome assembly:

1. Prepare GAF files for each pipeline: `prokka.gaf`, `bakta.gaf`, `pgap.gaf`, `gspa.gaf`.
2. Execute `gspa compare`:

```bash
gspa compare \
  -i genome.gff3 \
  --go-owl reference_data/go.owl \
  -k bacteria \
  -o comparison_report.html \
  prokka.gaf bakta.gaf pgap.gaf gspa.gaf \
  --lite
```

### Comparative metrics

`gspa compare` loads the reference ontology once, evaluates each GAF file
against the identical genome structure, and generates an interactive HTML dashboard:

```text
GSPA compare: 4 annotation sets
Evaluating: prokka...
  Completeness: 0.8333 (missing 4 essential functions)
  Mean IC: 5.210
Evaluating: bakta...
  Completeness: 0.9167 (missing 2 essential functions)
  Mean IC: 7.150
Evaluating: pgap...
  Completeness: 0.8750 (missing 3 essential functions)
  Mean IC: 6.840
Evaluating: gspa...
  Completeness: 0.9583 (missing 1 essential function)
  Mean IC: 8.421
Wrote: comparison_report.html
```

The rendered HTML file displays side-by-side radar plots, tabular metric summaries,
and breakdown lists of present versus missing essential functions.

## Step 4: generate HTML reports with gspa report

If you previously generated JSON reports from `gspa evaluate` across several
genomes or benchmark runs, compile them into an HTML summary with `gspa report`:

### Single-genome report

```bash
gspa report \
  -o ecoli_report.html \
  quality_ecoli.json
```

Produces an HTML dashboard detailing the genome's quality metrics, GO aspect
breakdown, and completeness status.

### Multi-genome cohort report

```bash
gspa report \
  -o cohort_quality.html \
  runs/genome1_quality.json \
  runs/genome2_quality.json \
  runs/genome3_quality.json
```

When given multiple JSON input files, `gspa report` automatically emits a
multi-genome cohort dashboard with sorting, filtering, and cross-isolate metric
distributions.

## Pitfalls

- **GAF identifier mismatches**: column 2 of the GAF file (DB Object ID) must
  match the CDS feature `ID=` or `protein_id=` attribute in the GFF3. If IDs
  differ, `gspa evaluate` cannot link annotations to genome features.
- **Root term inflation**: annotations to root terms (`GO:0008150`, `GO:0003674`,
  `GO:0005575`) have an Information Content of zero. High annotation rates driven
  exclusively by unspecific root terms will result in low `mean_ic` scores.
- **Incomplete essential function presets**: the `--kingdom` flag sets the default
  essential functions. If annotating an atypical genome (such as an obligate
  endosymbiont with a stripped-down genome), custom essential function lists
  can be supplied via `--essential-functions <tsv>`.
