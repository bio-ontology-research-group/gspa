# v3 metabolic benchmark data (items 8, 9, 12)

`benchmark/v3/prepare_metabolic_benchmarks.py` downloads public primary-source
datasets with `urllib`, validates every payload, hashes raw and derived files, and
normalises the measured phenotypes into labels and per-case experiment manifests.
Downloads go to `/tmp/gspa-v3-metabolic-data` (nothing large is committed).

```bash
PY=/tmp/gspa-v3-testenv/bin/python
$PY benchmark/v3/prepare_metabolic_benchmarks.py all [--out DIR] [--force]   # acquire + normalize
$PY -m unittest discover -s benchmark/v3/tests                                # 23 tests, no network
```

## What was acquired (run of 2026-10-26)

| Source | Status | Role |
|---|---|---|
| Maranas-group 2026 auxotroph validation, DOI 10.1371/journal.pcbi.1014469 (PMC13345470, CC-BY-4.0), supplement via the Europe PMC supplementary-files service | **acquired** | measured labels: S8 phenotypes (61 strain-substrate rows), S9 raw OD600 curves; S3-S7 are model designs and are **not** truth |
| Same paper, full-text XML | acquired | protocol and thresholds (Methods 2.6) |
| `maranasgroup/auxotroph-design` (iML1515 SBML, `carbon_substrate.csv`, design script), pinned to commit `01d975f5` | acquired | model used to design the strains; substrate to exchange map; design medium |
| Monk et al. 2017 (10.1038/nbt.3956), supplement datasets 1 and 11 | **blocked** | PMC returns a proof-of-work HTML page. It is rejected and recorded, **not solved or bypassed**. The Nature static-content URLs tried returned nothing and no author mirror of the Keio growth profiles was found. Monk labels are therefore **not available** |
| Harcombe et al. 2014 (10.1016/j.celrep.2014.03.070) measured co-culture data | **unavailable** | not open access in Europe PMC, PMC behind the same page, no machine-readable table in the author repositories reached. Nothing is filled in |
| Harcombe-lab COMETS models (E0 old iJO1366 metB deletion, S0 Salmonella, M0 AM1) and the notebook that builds them | acquired | **models only**. S0 methionine export is a model-imposed constraint (biomass altered to force 0.5 mmol Met per gDW), never a measured phenotype |

`acquisition_manifest.json` records, per source: every URL attempted and its
outcome, URL, DOI, license where stated, strain, role, retrieval time, size and
SHA-256, and (for GitHub files) the pinned commit or why none was pinned. Derived
files are hashed in the same manifest. The Europe PMC zip is re-generated on each
request, so its byte hash is not stable; the report therefore also records the
SHA-256 of each member workbook.

## Auxotroph validation labels (items 8 and 12)

Measured: Keio single-gene deletions inoculated into M9 with 4 g/L of one carbon
source in 24-well plates, 37 C, 24 h, OD600 (glucose positive control, no-carbon
water control). Per strain-substrate: dln_sub_net and dln_glu_net (water-corrected
ln OD change) and phi = dln_sub_net / dln_glu_net. Rule (paper): **growth** if
dln_sub_net >= 0.60 and phi >= 0.55; **partial** (slight growth) if >= 0.10 and
>= 0.20; otherwise **no_growth** (validated auxotrophy).

Acquired and normalised: **61 rows, 44 loci, 24 substrates; 36 no_growth, 17
partial, 8 growth** (matches the paper). All genes and substrates map
(39 exact symbols, 22 with the per-strain index such as `tpiA5` stripped; substrate
to exchange through `carbon_substrate.csv`, cross-checked against the design
script and the model). 61 cases are ready.

Quality checks recorded per row (`quality_report.json`):

- published class vs the paper rule applied to the published metrics: 0 conflicts
  (a conflict would make the label `unknown`);
- published metrics vs recomputation from raw S9 OD: 53 reproduced, 7 not
  reproducible because S9 has **no water control** for those strains (all D-mannose
  rows but one), 1 (xylB on D-xylose) differs and is unexplained. The audit never
  changes a label; treat the 8 as a sensitivity subset;
- 10 rows are published as borderline and their class changes somewhere in the
  threshold grid (phi_low, phi_high, dln_sub); the label set across the grid is in
  `labels_across_threshold_grid`.

Truth versus prediction. Labels come only from measured OD columns. Any model
column (`wt_growth`, `rescued_biomass`, `pct_rescue`, ...) is dropped and listed in
`model_columns_dropped`. Missing or unnormalisable data (no substrate measurement,
glucose control that did not grow, ambiguous gene symbol) is `unknown` or
`excluded`, never `no_growth`. `partial` is never forced into a binary class
(`binary_growth` is empty for it).

**Limits that matter for use.**

- *Selection.* All 61 strains were chosen because the iML1515-based design
  predicted that the deletion abolishes growth on the substrate. The set has no
  model-predicted growers, so it estimates the precision of predicted auxotrophy
  (here 36 of 61 confirmed, 17 partial, 8 false) and supports no recall,
  specificity or full confusion matrix.
- *Model overlap and masking.* iML1515 GPRs designed these strains. Using them
  for a **phenotype** (knockout-growth) test is allowed. They must **not** be used
  for **independent gene-function recovery**. `masking.json` lists, per case, the
  reactions whose GPR contains the locus; a recovery run must remove the locus from
  those GPRs and from every sequence, domain, neighborhood and recurrence input.
  Nothing here is fitted or tuned on the labels.
- *Medium.* The experimental medium is M9 plus 4 g/L carbon. The model medium
  (ALWAYS_OPEN salts at the model default, O2 15, carbon 10 mmol/gDW/h, read from
  the design script) is a modeling assumption and is recorded as such in
  `model_metadata.json`. Fluxes are mmol/gDW/h; the 4 g/L concentration has no
  flux equivalent.
- *Strain.* Keio deletions; the parent (BW25113) is not restated in the paper.
- *Partial auxotrophy* (17 rows) is a measured intermediate; a static FBA call
  cannot reproduce graded growth, so score it separately.

## Exports (`<out>/derived/`)

| File | Content |
|---|---|
| `auxotroph_validation_labels.csv` | one row per strain-substrate: IDs and mappings with their mode, measured metrics, published and recomputed labels, final `label`, `binary_growth`, borderline and grid labels, raw-OD audit, `status` |
| `cases.jsonl` | one line per ready case: the row plus a `metabolic_experiments` manifest (matched wild-type `carbon_source:` control and the Keio `knockout:` under that carbon source; declared external exchanges; closed-by-default medium; explicit biomass and ATPM maintenance 6.86) |
| `model_metadata.json` | model file and hash, biomass, exchanges, medium and its provenance |
| `masking.json` | per-case masking requirements |
| `community_registry.json` | Harcombe models (hash, reaction counts, objective, notes) and the unavailable measured data |
| `quality_report.json` | all counts above |

Scoring rule for a case: predicted growth is `max_growth >= growth_threshold` of the
knockout experiment (absolute growth on the substrate; the `growth_ratio` is
relative to glucose and is not used). The wild-type control on the same substrate
shows whether the model can grow on it at all. The generated manifests need the
root fix `822ab69` for GPR reading: on the earlier reader the `knockout:` item did
not resolve (checked, the carbon-source item and boundary validation were fine).
No LP benchmark was run here; the root runs them.

## Data dependencies still open

- Monk 2017 Keio growth and knockout profiles (item 8 breadth): needs a
  retrievable copy of Supplementary Data set 11 (and 1 for the models) from a
  publisher or author repository.
- Harcombe 2014 measured co-culture growth and species ratios, with units, for the
  metB / Salmonella / Methylobacterium consortium (item 9). A static
  balanced-growth model cannot validate the spatial and temporal dynamics, so only
  comparable endpoint phenotypes should be registered, and none has been acquired.
- Measured exchanged metabolites and transporter-knockout phenotypes for item 9
  are not in any acquired source.
- Characterised secondary-metabolite production and held-out organisms with
  measured growth across media (item 12 beyond E. coli) are not covered by these
  sources.
