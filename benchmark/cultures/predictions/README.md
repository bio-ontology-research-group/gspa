# Dark-matter prediction results — empty-quarter cultures

Context-based functional predictions for four user cultures from
`/data/emptyquarter/sequencing-results` on unimatrix01, generated
without re-running any wet-lab work. For each reaction that gapsmith
could not fill in the culture's metabolic network, the pipeline ranks
candidate proteins by the density of annotated pathway-neighbour
enzymes in their genomic neighbourhood, with proteins already assigned
to other reactions explicitly excluded from the candidate pool.

Two files of interest:

- `{tag}_dark_matter.tsv` — the raw per-gap predictions for culture
  `{tag}`. Top-5 candidates per unfilled reaction.
- **`validation_candidates_pident.tsv`** — the wet-lab-ready shortlist:
  predictions where the candidate protein is genuine dark matter by the
  standard < 30% sequence-identity criterion, the candidate is specific
  (predicted for ≤ 5 distinct gap reactions), and at least one prediction
  is at rank ≤ 3.
- `validation_candidates.tsv` — equivalent shortlist using an ESM2
  cosine percentile definition instead of sequence identity. Included
  for completeness; `validation_candidates_pident.tsv` is the file to
  use.

## Cultures included

| tag     | source                                                   | proteins called |
|---------|----------------------------------------------------------|-----------------|
| MR59-1  | isolate, site59 (surface water, TSB medium)              | 2,394           |
| MR60-1  | isolate, site60                                          | 5,481           |
| C-1.1   | enrichment C-1, MAG 1 (Bacillus paralicheniformis-like)  | 4,257           |
| C-1.3   | enrichment C-1, MAG 3                                    | 3,238           |

## Pipeline

For each culture genome (assembly FASTA → predictions):

1. **Gene calling** — `prodigal -p meta` over the assembly FASTA,
   producing a proteome and per-gene coordinates.
2. **Sequence annotations** — DIAMOND blastp against the reference
   Swiss-Prot DIAMOND DB (`/data/hohndor/gspa/benchmark/benchmark_data/reference_db.dmnd`)
   plus `hmmsearch --cut_ga` against Pfam-A. Results combined via
   `benchmark-py/02b_parse_predictors_to_claims.py` using the full
   UniProt-GOA as the source of GO annotations on Swiss-Prot hits.
3. **Gapsmith gap identification** —
   `gapsmith find -p all -t Bacteria -A diamond` on the proteome,
   producing a per-pathway / per-reaction status table with
   `good_blast` / `bad_blast` / `no_blast` / `no_seq_data` /
   `spontaneous` labels.
4. **Evidence integration** — `gspa integrate --lite` produces
   per-(protein, GO) integrated posteriors accounting for essential-
   function priors, pathway context, GO-ontology propagation, etc.
5. **Dark-matter prediction** — `predict_dark_matter.py`
   (see `benchmark/cultures/predict_dark_matter.py`) iterates every
   gap reaction (status `bad_blast` or `no_blast`, has EC + pathway
   context) and, for each:
   1. BFS through the SEED reaction graph from the gap reaction's
      EC-equivalents (depth 2, decay α=0.5, currency metabolites
      pruned).
   2. Maps the neighbour reactions to their GO terms via ec2go.
   3. Identifies *anchors*: proteins in the culture with integrated
      posterior ≥ 0.3 on any neighbour GO (the "A, B, C, D already
      catalyse neighbour steps" proteins).
   4. Computes the Gaussian kernel density of anchor weights
      (bandwidth 5 kb, truncated at 4σ) at every gene position on
      each contig carrying anchors.
   5. Emits the top-5 genes by density, **excluding any gene that is
      itself already assigned by gapsmith to another reaction or is
      itself an anchor for this reaction's neighbours**. The excluded
      set is anchors ∪ `good_blast` proteins (minus any gapsmith
      assignment for the *same* reaction, which we want to keep as a
      legitimate candidate).

All steps run as SLURM jobs on unimatrix01. The full pipeline for a
4 k-protein genome takes ~30 min (DIAMOND + Pfam dominate).

## Column definitions (both TSVs)

| column | meaning |
|---|---|
| culture | MR59-1 / MR60-1 / C-1.1 / C-1.3 |
| gap_pathway | MetaCyc pathway ID as reported by gapsmith |
| gap_rxn | MetaCyc / SEED reaction ID of the unfilled step |
| gap_ec | EC number for the gap reaction (4-digit) |
| rank | Candidate rank within this gap (1 = top by density) |
| candidate | Prodigal-called protein ID (contig-prefixed) |
| density | Sum of anchor-weight × Gaussian-kernel at candidate position |
| n_anchors | Number of neighbour-reaction anchors within ±20 kb of the candidate on its contig |
| n_nbr_gos | Number of neighbour-reaction GO terms used to form the anchor set |
| contig / pos | Candidate's contig + midpoint coordinate |
| near_anchors | Three nearest anchors to candidate (position + weight) |

Additional columns in `validation_candidates_pident.tsv`:

| column | meaning |
|---|---|
| max_pident | Max % identity from DIAMOND claims: candidate against any GO-annotated Swiss-Prot entry. 0 ⇒ no hit at all. < 30 is the dark-matter cutoff used for filtering. |
| max_cos_esm2 | Max cosine between the candidate's ESM-2 t30 embedding and any of 2,259 panel-wide per-EC centroids. Context signal; > 0.9 is typical base rate. |
| esm2_argmax_ec | EC whose centroid gave that max cos. If different from `gap_ec`, ESM2 alone would have predicted a *different* function. |
| n_gaps_for_candidate | # distinct gap reactions this candidate appears as top-5 for. Low = specific prediction; high ⇒ positional hub (see caveats below). |

## Dark-matter definition

We use the standard **< 30 % identity** criterion: a protein is dark if
its maximum DIAMOND pident to any Swiss-Prot entry with any GO
annotation is below 30 %. Proteins with no DIAMOND hit at all (above
the 1e-5 e-value threshold) count as `pident = 0` and are the darkest
cases.

Across the four cultures:

| culture | proteins | dark (< 30 %) | context-predicted (rank 1–5) | rank-1 predicted |
|---|---|---|---|---|
| MR59-1 | 2,394 | 831 (34.7 %) | 213 | 64 |
| MR60-1 | 5,481 | 2,173 (39.6 %) | 315 | 83 |
| C-1.1  | 4,257 | 1,170 (27.5 %) | 385 | 119 |
| C-1.3  | 3,238 | 1,521 (47.0 %) | 337 | 103 |

## How to read the shortlist

`validation_candidates_pident.tsv` contains 1,525 prediction rows over
273 unique `(culture, candidate)` pairs. Filter recipe used to build
the shortlist:

1. `max_pident < 30` (dark matter in the standard sense).
2. `n_gaps_for_candidate ≤ 5` (candidate is predicted for only a
   handful of reactions — filters out positional hubs).
3. `best rank ≤ 3` (context has a strong signal for at least one gap).

High-priority wet-lab targets have:

- `max_pident = 0` — no sequence-based annotation of any kind exists.
- `n_gaps_for_candidate = 1` — context picks them for a single,
  specific reaction.
- `rank = 1` and `n_anchors ≥ 10` — tight operon with many
  pathway-neighbour catalysts nearby.
- `esm2_argmax_ec ≠ gap_ec` — ESM2 alone would have predicted a
  different function, so a positive wet-lab result would show context
  adding information beyond what sequence-embedding methods provide.

Examples that satisfy all four:

| culture | candidate | gap_rxn | gap_ec |
|---|---|---|---|
| MR59-1 | contig_1_417 | 4-NITROPHENOL-2-MONOOXYGENASE | 1.14.13.29 |
| MR59-1 | contig_1_1640 | 2.7.1.121-RXN | 2.7.1.121 |
| MR60-1 | contig_1_4977 | PCP4MONO-RXN | 1.14.13.50 |
| MR60-1 | contig_1_4593 | GLYOXII-RXN | 3.1.2.6 |

## Caveats (read before prioritising)

1. **Positional hubs.** A small subset of candidates per culture show
   up as top-1 for 100+ different gap reactions. These are proteins
   sitting in a genuinely high-density region of the genome (big
   biosynthetic operon cluster) that will always win by density
   regardless of what reaction we ask about. The
   `n_gaps_for_candidate ≤ 5` filter removes them from the shortlist.
2. **EC centroid coverage.** The ESM-2 centroid set covers 2,259 EC
   numbers built from the 29-genome GTDB panel (excluding mg1655).
   When a candidate's gap_ec has no centroid (~24 % of gaps), the
   ESM2 columns are still informative (`max_cos_esm2` is over *all*
   2,259 ECs, `esm2_argmax_ec` is the best match regardless of target).
3. **The context signal requires anchors.** Reactions whose neighbours
   also have no annotated catalysts in the culture will produce no
   candidates at all. Roughly half the unfilled reactions across these
   cultures fall into this category and do not appear in the output —
   context cannot help where no context exists.
4. **Predictions are at operon granularity.** Within a dense operon of
   3–5 genes, density ranks the highest-density gene first; the
   "correct" dark-matter member may be at rank 2 or 3. Inspect all
   top-3 candidates in the output, not just rank-1.
5. **Not validated against orthogonal signals.** No cross-genome
   non-anchor LR scoring has been applied yet (the `nonanchor_catalog`
   infrastructure exists but isn't wired into these predictions).

## Reproducibility

The full pipeline, inputs, and outputs live on unimatrix01 under
`/data/hohndor/gspa/proteomes/cultures/` with subdirectories per
culture containing `prodigal/`, `preds/`, `claims/`, `layout/`,
`gapsmith/`, `integrated/`, and `plm/`. Scripts are mirrored in
`benchmark/cultures/` in this repository.

To reproduce for a new culture, add a row to
`benchmark/cultures/manifest.tsv` (`tag<TAB>fasta_path`) and submit
`annotate_culture.sh` as a SLURM array job, followed by
`predict_dark_matter.py`.

## Reports & validation

`benchmark/leave_reaction_out/BASELINE_REPORT.md` documents the
methodological validation of the approach on two LRO benchmarks
(*E. coli* MG1655 and *B. subtilis*, as panel-known genomes with
ground truth) — including the measurements that revealed the
operon-granularity ceiling and the anchor + gapsmith-assigned
exclusion trick that drives the rank-1 lift on mg1655 from 0.041
(density alone) to 0.196 (filtered).
