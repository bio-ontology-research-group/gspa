# Interactive browser

Use this workflow to build a self-contained, zero-dependency HTML5 genome browser
with `gspa visualize`. 

Genome annotation pipelines typically produce fragmented outputs: GFF files,
TSV tables, JSON logs, and separate predictor outputs. `gspa visualize` consolidates
the entire GSPA workspace into a single portable HTML file that can be opened
in any modern web browser. The generated file runs entirely client-side without
requiring a web server, backend database, or internet connection.

## What you need

`gspa visualize` reads the conventional workspace layout produced by
`gspa annotate` and `gspa integrate`:

```text
workdir/
  prokka_out/{genome}.{tsv,gff}                 per-CDS metadata and coordinates
  gspa_out/                                     integrator outputs
    integrated.tsv                              calibrated posterior probabilities
    provenance.json                             supporting predictors and raw claims
    operons.tsv                                 predicted operon clusters
    quality_gspa.json                           GAEF quality assessment report
  mdf_out/{genome}.mdf.tsv                      mDeepFRI sidecar outputs (optional)
  proteinfer_out/{genome}.proteinfer.tsv        ProteInfer sidecar outputs (optional)
  clean_out/{genome}.clean.tsv                  CLEAN EC predictions (optional)
  amrfinder_out/{genome}.amr.tsv                AMRFinderPlus resistance hits (optional)
  antismash_out/regions.js                      antiSMASH secondary metabolite BGCs (optional)
  ipr_out/{genome}.faa.tsv                      InterProScan domain hits (optional)
  input/{genome}_assembly.fa                    assembly FASTA for igv.js (optional)
```

## Running gspa visualize

Execute `gspa visualize` pointing to your workspace directory:

```bash
gspa visualize \
  -w workspace/ \
  -g MYISO \
  -o MYISO_browser.html \
  --fasta workspace/input/MYISO_assembly.fa \
  --go-obo reference_data/go.obo \
  --ec2go reference_data/ec2go.txt
```

### CLI options

- `-w, --workdir`: root workspace directory containing `gspa_out/`, `prokka_out/`, and sidecars.
- `-o, --out`: output path for the self-contained HTML file (default: `<workdir>/<run-dir>/<genome>_browser.html`).
- `-g, --genome-id`: genome identifier used to match sidecar filenames (e.g. `MYISO.mdf.tsv`). Defaults to the workdir folder name.
- `-r, --run-dir`: subdirectory within workdir containing integrator outputs (default: `gspa_out`).
- `--fasta`: assembly nucleotide FASTA. When supplied, the sequence is encoded as a base64 data URL so the embedded `igv.js` track functions offline.
- `--go-obo`: path to `go.obo` for human-readable GO term names and definitions.
- `--ec2go`: path to `ec2go.txt` for enzyme name lookup.
- `--python`: Python interpreter to execute the templating script (default: `python3`).

### Real output

```text
GSPA visualize
  Workdir:   workspace/
  Run dir:   gspa_out
  Genome id: MYISO
  Output:    MYISO_browser.html
... Loading prokka_out/MYISO.tsv (3840 CDS)
... Loading gspa_out/integrated.tsv (1757 annotations)
... Loading gspa_out/operons.tsv (153 operons)
... Embedding reference FASTA (4.2 MB)
... Rendering self-contained HTML bundle
Wrote HTML browser (12.4 MB) -> MYISO_browser.html
```

## Browser layout and interactive tabs

The generated HTML browser features eight tabs:

### 1. Proteins
The primary exploration table. Lists all predicted coding sequences with locus
tags, product names, gene coordinates, and integrated GO/EC annotations.
Clicking any protein opens a slide-out drawer displaying:

- **Calibrated posteriors**: final probability alongside likelihood log-odds and prior bonuses.
- **Evidence provenance**: exact raw scores from each supporting predictor (DIAMOND, FoldSeek, ESM2, InterPro).
- **Enforcement audit**: records if the term was relaxed, promoted, or pruned during Stage 2.

### 2. Functions
A Gene Ontology exploration interface. Provides interactive DAG visualization
showing the hierarchical path from assigned molecular functions up to the root
categories. Includes frequency breakdown charts across Molecular Function (MF),
Biological Process (BP), and Cellular Component (CC).

### 3. Genome browser
An embedded offline `igv.js` track viewer. Visualizes contigs and scaffolds with
interactive tracks:

- CDS features colored by functional category or posterior confidence.
- Operon blocks indicating multi-gene transcriptional units.
- Overlaid annotations for antiSMASH biosynthetic gene clusters (BGCs) and AMRFinderPlus loci.

### 4. Operons
An operon map viewer. Displays co-directional gene clusters predicted by
`OperonPredictor` with intergenic distances and joint functional coherence
scores. Users can search by operon ID or gene name to inspect cluster synteny.

### 5. Pathways
Metabolic pathway completion matrices for KEGG and MetaCyc pathways.
Displays pathway completeness percentages, highlighting present reactions,
missing steps, and candidate genes identified by the dark-matter and RLGC suggesters.

### 6. Special features
Dedicated views for specialized genomic elements:

- **AMRFinderPlus**: detected antimicrobial resistance genes, point mutations,
  and virulence factors with target drug classes.
- **antiSMASH**: biosynthetic gene clusters (BGCs) categorized by secondary
  metabolite type (non-ribosomal peptides, polyketides, terpenes).

### 7. Quality (GAEF)
Displays the complete Genome Annotation Evaluation Framework report:

- Radar chart of the five quality dimensions (completeness, process coherence,
  complex coherence, consistency, information content).
- Table of present versus missing domain-essential functions.
- Taxonomic consistency certification status.

### 8. Pipeline
A complete provenance and audit trail. Documents software versions, database
release timestamps, execution timing across each pipeline phase, and overall
per-predictor claim yield.

## Offline sharing and publication

Because the output file is entirely self-contained:

- You can email it to collaborators or host it as a static file on GitHub Pages,
  AWS S3, or an institutional web server.
- No database configuration, Node.js server, or backend CGI script is required.
- Opening the file locally (`file:///path/to/MYISO_browser.html`) works identically
  to viewing it on a web server.

## Pitfalls

- **Missing assembly FASTA**: if `--fasta` is omitted, the Genome Browser tab
  will display gene coordinate tracks but cannot display nucleotide sequence
  bases or three-frame translations.
- **Large multi-gigabase genomes**: for large eukaryotic assemblies, base64-encoding
  the entire genome into the HTML file can create files exceeding 100 MB, which
  may slow down browser rendering. For large genomes, visualize specific scaffolds
  or omit `--fasta`.
- **Locus tag mismatches**: ensure that Prokka TSVs and GSPA integrated TSVs use
  the same gene identifier column. If locus tags differ, the browser cannot link
  physical coordinates to functional calls.
