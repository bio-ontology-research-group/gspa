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
| reaction_name | Human-readable reaction name (from gapsmith Reactions.tbl, fallback SEED reactions.tsv) |
| ec_name | EC description (extracted from ec2go.txt) |

Additional columns in `validation_candidates_pident.tsv`:

| column | meaning |
|---|---|
| max_pident | Max % identity from DIAMOND claims: candidate against any GO-annotated Swiss-Prot entry. 0 ⇒ no hit at all. < 30 is the dark-matter cutoff used for filtering. |
| max_cos_esm2 | Max cosine between the candidate's ESM-2 t30 embedding and any of 2,259 panel-wide per-EC centroids. Context signal; > 0.9 is typical base rate. |
| esm2_argmax_ec | EC whose centroid gave that max cos. If different from `gap_ec`, ESM2 alone would have predicted a *different* function. |
| n_gaps_for_candidate | # distinct gap reactions this candidate appears as top-5 for. Low = specific prediction; high ⇒ positional hub (see caveats below). |
| orthogroup | Panel orthogroup assigned to the candidate via MMseqs2 easy-search against the 29-genome GTDB panel (≥50 % id, ≥80 % query coverage). `unclustered` ⇒ no panel homolog at that cutoff — this is the norm for truly dark candidates. |
| seed_rxn_lookup | Comma-separated SEED reaction IDs used to look up cross-genome evidence for this gap (gap_rxn + EC-alias equivalents that exist in the catalog). |
| n_sig_nonanc / n_sig_total | # panel genomes (out of 29) where the R-signature is present AND this orthogroup has a non-anchor member inside / total genomes with R-signature. Empty ⇒ no catalog match. |
| n_base_with / n_base_total | # panel genomes where this orthogroup has any member / panel size. Gives the baseline rate to compare against. |
| log_lr | log₁₀ of `[(n_sig_nonanc / n_sig_total + ε) / (n_base_with / n_base_total + ε)]`. Positive ⇒ orthogroup appears in R-signature contexts more often than its baseline prevalence — cross-genome evidence supports the prediction. 0 ⇒ no signal. Negative ⇒ under-represented. |
| n_catalog_lookups | # SEED equivalents of gap_rxn that had a catalog entry for this orthogroup. |

### Cross-genome evidence (29-genome panel × culture candidates)

Cross-genome support is computed by:
1. MMseqs2 easy-search of every culture protein against the
   concatenated 29-genome panel proteome, taking the best hit at
   ≥50 % seq id and ≥80 % query cov → assigns each culture protein to
   the panel orthogroup of its best hit (or `unclustered`).
2. For the gap reaction, resolving SEED equivalents via EC.
3. Looking up `(orthogroup, SEED_rxn)` in
   `nonanchor_catalog_cultures.tsv` (43 M rows, 5,500 SEED reactions ×
   29 panel genomes, restricted to SEED equivalents of the 5,238 unique
   culture gap ECs).

Of the 1,525 validation rows:

| step | count | % |
|---|---|---|
| predictions in shortlist | 1,525 | 100 |
| candidate has ≥50% id panel homolog (orthogrouped) | 301 | 19.7 |
| catalog match (orthogroup × SEED rxn) | 136 | 8.9 |
| `|log_lr| > 0.1` (meaningful cross-genome signal) | 63 | 4.1 |
| positive `log_lr > 0.1` (cross-genome *supports* prediction) | 61 | 4.0 |

Low orthogrouped fraction is expected: the shortlist is dark
candidates (pident < 30 % to Swiss-Prot), so most also lack panel
homologs at 50 % id. Cross-genome evidence is most useful exactly for
the subset that *does* have a distant panel ortholog.

### Top cross-genome-supported dark-matter predictions

Candidates with `pident = 0` (no Swiss-Prot annotation homolog at all)
AND `log_lr > 0.5` (cross-genome strongly supports):

| culture | candidate | gap reaction | EC | orthogroup | n_sig_nonanc / n_sig_tot | log_lr |
|---|---|---|---|---|---|---|
| MR59-1 | contig_1_908 | Neuraminidase | 3.2.1.18 | A0A6B0BKS9 | 20 / 30 | **1.29** |
| MR59-1 | contig_1_908 | N-acetylneuraminate synthase | 2.5.1.56 | A0A6B0BKS9 | 5 / 27 | 0.73 |
| MR60-1 | contig_1_88 | palmitoyl-ACP 9-desaturase | 1.14.19.2 | A0A3M5EC75 | 5 / 29 | 0.70 |
| MR59-1 | contig_1_1559 | UDP-NAM-pentapeptide lysyltransferase | 2.3.2.10 | A0A0D1KVH5 | 18 / 30 | 0.65 |
| MR60-1 | contig_1_617 | CMP-N-acetylneuraminate synthetase | 2.7.7.43 | A0A367M9C3 | 4 / 27 | 0.64 |
| MR60-1 | contig_1_3647 | allantoinase | 3.5.2.5 | A0A3M5E3X5 | — | 0.62 |

These are the strongest wet-lab candidates in the shortlist:
**pident = 0** (no sequence annotation) AND context puts the gene into
the right operon AND across the panel the same orthogroup
consistently appears as a non-anchor in that reaction's context
windows. A positive wet-lab result would show three independent lines
of evidence (genomic context within culture, cross-genome
non-anchor consistency, operon-level localisation) combining to
predict a function for a protein with no sequence-level annotation.

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
5. **Cross-genome evidence is now included in
   `validation_candidates_pident.tsv`**, but only for the 19.7 % of
   shortlist candidates that have any panel ortholog at 50 % id /
   80 % cov. For the remaining 80 % (true dark matter w.r.t. the
   panel), the prediction rests on within-culture context alone.
   Raw per-gap predictions (`{tag}_dark_matter.tsv`) are still
   within-culture-only.

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
