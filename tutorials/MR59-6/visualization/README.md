# MR59-6 annotation browser

Interactive HTML viewer for the GSPA annotation set on the desert *Pontibacter*
genome MR59-6.

## How to open

```bash
xdg-open MR59-6_browser.html       # Linux
open MR59-6_browser.html           # macOS
```

The file is self-contained — no server, no network. Just double-click.

## What you can do in it

- **Header**: GTDB taxonomy, environmental sample metadata, eight headline metric
  cards (genome size, CheckM2 / BUSCO completeness, GSPA coverage, total
  hypotheses, high-confidence hypotheses, GAEF Composite, process coherence).
- **Proteins tab**: 4,372-row virtualised table with per-protein roll-ups
  (best GO term, max posterior, count of high-conf annotations, P/F/C aspect
  breakdown, AMR/BGC/membrane/secreted/enzyme flags, # contributing tools).
  Filters: full-text search, posterior threshold slider, GO aspect toggles,
  predictor-of-origin, "show only" categorical (hypothetical, AMR, BGC,
  membrane, secreted, enzyme, has high-conf, no high-conf). Click a row → the
  detail panel slides in showing every (protein, GO) hypothesis grouped by
  aspect, sorted by posterior, each with its confidence band and the predictors
  that supported it as chips with multiplicity.
- **Functions tab**: confidence histogram across all 84,215 hypotheses, GO
  aspect donut for the high-conf subset, top-60 most-frequent GO terms.
- **Genome map tab**: linear plot of the single circular contig with all CDS
  drawn on +/− strands and BGC regions overlaid as transparent bands.
- **Special features tab**: AMR table (1 hit — subclass B1 metallo-β-lactamase),
  BGC table (4 NRPS / PKS clusters with the CDS in each region), and a
  subcellular signature derived from high-conf cellular-component GO terms
  (since SignalP / TMHMM weren't part of the IPR run for this tutorial).
- **Quality tab**: GAEF metric bars (coverage, completeness, three coherences,
  composite) with the missing essential GO term called out.
- **Pipeline tab**: predictors that ran (with ✓), predictors documented as
  blocked (PSORTb container, DeepEC transformers chain), and a bar chart of
  how many supporting hits each predictor contributed to the integrated set.

## How it was built

`make_viz.py` reads the cluster workspace at
`/data/hohndor/gspa-tutorial-MR59-6/`, joins integrated.tsv with
provenance.json, looks up GO term labels in `go.obo`, and emits a single
HTML file with the data embedded as JSON. ~5.6 MB on disk.
