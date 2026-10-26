# v3 metabolic benchmark data (items 8, 9, 12)

`benchmark/v3/prepare_metabolic_benchmarks.py` downloads public primary-source
datasets with `urllib`, validates every payload, hashes raw and derived files, and
normalises the measured phenotypes into labels and per-case experiment manifests.
Downloads go to `/tmp/gspa-v3-metabolic-data` (nothing large is committed).

```bash
PY=/tmp/gspa-v3-testenv/bin/python
$PY benchmark/v3/prepare_metabolic_benchmarks.py all [--out DIR] [--force]   # acquire + normalize
$PY benchmark/v3/prepare_wintermute_coculture.py                               # measured co-culture pairs (after acquire)
$PY benchmark/v3/prepare_mibig_production.py                                   # characterised production observations
$PY -m unittest discover -s benchmark/v3/tests                                # no network
```

Sections: manuscript protocol correction, auxotroph labels, community data (item 9), secondary-metabolite
production observations (item 12), what the root can run, open dependencies.

## What was acquired (run of 2026-10-26)

| Source | Status | Role |
|---|---|---|
| Maranas-group 2026 auxotroph validation, DOI 10.1371/journal.pcbi.1014469 (PMC13345470, CC-BY-4.0), supplement via the Europe PMC supplementary-files service | **acquired** | measured labels: S8 phenotypes (61 strain-substrate rows), S9 raw OD600 curves; S3-S7 are model designs and are **not** truth |
| Same paper, full-text XML | acquired | protocol and thresholds (Methods 2.6) |
| `maranasgroup/auxotroph-design` (iML1515 SBML, `carbon_substrate.csv`, design script), pinned to commit `01d975f5` | acquired | model used to design the strains; substrate to exchange map; design medium |
| `knockout_utils.py` and `03_single_knockout.py` of the same repository, same pinned commit | acquired | the authors' knockout stage: `ATPM_MAINT = 8.39`, `O2_RATE = 20.0`, `UPTAKE_RATE = 10.0`, and the code that applies the closures |
| Same paper, S1 and S2 tables (inside the supplement zip) | acquired | per-substrate regulatory reaction closures and the growth they imply (see below) |
| Monk et al. 2017 (10.1038/nbt.3956), supplement datasets 1 and 11 | **blocked** | PMC returns a proof-of-work HTML page (not solved or bypassed); europepmc.org and static-content.springer.com answer 403; nature.com answers with a bot-check page; Europe PMC lists the manuscript as not open access; GitHub, Zenodo and figshare searches found no author mirror. Monk labels are **not available** |
| Harcombe et al. 2014 (10.1016/j.celrep.2014.03.070) measured co-culture data | **unavailable** | Europe PMC: not open access (supplement and full text both refused); the Harcombe-lab repositories reached (COMETS protocols; Dryad/Zenodo 10.5061/dryad.5jk4d) hold model and layout inputs only. Nothing is filled in |
| Wintermute & Silver 2010 (10.1038/msb.2010.66, PMC2964121), supplement `msb201066-s2.xls` | **acquired** | alternative established community dataset with measured growth and explicit medium (see below) |
| MIBiG 4.0 JSON (`dl.secondarymetabolites.org`, CC0) | **acquired** | positive characterised-production observations only (see below) |
| Harcombe-lab COMETS models (E0 old iJO1366 metB deletion, S0 Salmonella, M0 AM1) and the notebook that builds them | acquired | **models only**. S0 methionine export is a model-imposed constraint (biomass altered to force 0.5 mmol Met per gDW), never a measured phenotype |

`acquisition_manifest.json` keeps an `attempt_log` of every URL tried in every run, failed attempts included
(`acquire --only` merges into the existing manifest instead of replacing it), and records, per source: every URL attempted and its
outcome, URL, DOI, license where stated, strain, role, retrieval time, size and
SHA-256, and (for GitHub files) the pinned commit or why none was pinned. Derived
files are hashed in the same manifest. The Europe PMC zip is re-generated on each
request, so its byte hash is not stable; the report therefore also records the
SHA-256 of each member workbook.

## Manuscript protocol correction (primary source)

The frozen development run used the repository defaults (ATPM 6.86, O2 15, carbon 10 mmol/gDW/h, unregulated
iML1515). That is **not the paper configuration**. Methods 2.1-2.3 of the paper and the authors' own knockout code
(`knockout_utils.py`, same pinned commit) give **ATPM 8.39, O2 20, carbon 10**, and apply substrate-specific regulatory
reaction closures (rFBA, S2 Table) in the wild-type control and the knockout alike. The paper and the code agree with each
other; only the legacy numbers differ. The values are read from the paper text (the quoted sentence is stored) and from
the code by `ast` (never executed); none was chosen using the labels.

Three registered protocols, separate case files, frozen legacy cases untouched (`cases.jsonl` is byte-identical to the
registered baseline, SHA-256 `58504287...`; a test rebuilds it with and without the new sources):

| Protocol | ATPM | O2 | Carbon | Closures | File | Engine |
|---|---:|---:|---:|---|---|---|
| `legacy_repository_defaults` | 6.86 | 15 | 10 | no | `cases.jsonl` | supported (frozen baseline, not the paper) |
| `manuscript_unregulated` | 8.39 | 20 | 10 | no | `cases_manuscript_unregulated.jsonl` | **supported**: runnable manifests that differ from legacy only in the ATPM bound and the O2 limit |
| `manuscript_regulatory` | 8.39 | 20 | 10 | S2 | `cases_manuscript_regulatory.jsonl` | **unsupported until wired** |

`protocol_registry.json` stores, per parameter, its provenance (quote or constant name, source hash), the discrepancy
list against the legacy run, and the closure semantics. `regulatory_closures.json` holds the closures of 46 substrates
(244-298 closed reactions each; 303 unique reactions, equal to the 303 stated in the paper; the S1 counts agree with S2
for every substrate) with the hashes of the S1 and S2 members.

Regulatory records carry **typed sidecar constraints** (`regulatory_reaction_closure`: lower = upper = 0 on each listed
reaction, applying to both matched experiments, only S2 rows flagged `applied`; the `relaxed_FlagB` rows were relaxed by
the authors to keep the model viable and are listed, not closed). `scripts/metabolic_experiments.py` has no
baseline-condition closure, so these records have **no `manifest` key**: a runner that is not wired fails the case
visibly instead of running it as an unregulated one. The base manifest is stored as `manifest_without_closures`.
Closures are not encoded as gene knockouts. All 61 strain/substrate cases are fully specified (24 substrates; 21
distinct closure sets among the cases). A substrate without an S2 set, or a closure naming a reaction absent from the
model, becomes `excluded:<reason>` and is never guessed (neither occurs in the real data).

Truth independence. Case inclusion depends on gene and substrate mapping only, never on a label. Manifests and
constraints are built from `case_id`, locus and substrate exchange; labels live only under `case`. Each record has an
`inputs_sha256` over its manifest and constraints. A regression test flips every measured label and checks that the case
set, manifests, constraints and digests do not change (only the labels do).

Two caveats that sit in the primary sources and are not resolved here. (1) S2 was generated by design script 01 under
O2 15 and the model-default ATPM, while the paper text and the knockout code apply it at O2 20 / ATPM 8.39. (2) The S1
growth values (0.6919 / 0.6923 on glucose, ...) were computed under script-01 conditions; the root can use them as a
sanity check of the wiring, not as a target the manuscript conditions must reproduce exactly. They are stored as
`paper_reported_growth` and are model outputs, never labels.

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
- *Medium.* The experimental medium is M9 plus 4 g/L carbon. The legacy model medium
  (ALWAYS_OPEN salts at the model default, O2 15, carbon 10 mmol/gDW/h, read from
  the design script; the manuscript protocols above use O2 20 and ATPM 8.39) is a modeling assumption and is recorded as such in
  `model_metadata.json`. Fluxes are mmol/gDW/h; the 4 g/L concentration has no
  flux equivalent.
- *Strain.* Keio deletions; the parent (BW25113) is not restated in the paper.
- *Partial auxotrophy* (17 rows) is a measured intermediate; a static FBA call
  cannot reproduce graded growth, so score it separately.

## Community data (item 9)

Harcombe 2014 measured data could not be retrieved, so the alternative established dataset is
`prepare_wintermute_coculture.py`: Wintermute & Silver 2010, 46 conditionally lethal Keio auxotrophs, every pair
co-cultured in M9 glucose minimal medium in which no strain grows alone, total cell density measured daily for four
days in two replicates (`wintermute2010_pairs.csv`: 1,035 unordered pairs plus 46 self-pair controls;
`wintermute2010_cases.jsonl`: one `inputs` + `observed` record per pair). The supplement is a binary `.xls`, read by
the dependency-free `biff_reader.py` (tested against a generated BIFF/OLE2 fixture).

Labels use the paper's own thresholds on the mean of the four day-4 readings (fold growth = density / 1.5e7 cells/ml):
**synergy** above 50-fold, **no_synergy** below 8-fold (the bound the paper reports for self-pairs), **intermediate**
in between, never forced to binary. Result: 178 synergy, 430 no_synergy, 427 intermediate. **Data check:** the
supplement reproduces the paper's "17% of 1,035 pairs above 50-fold" (178 / 1,035 = 17.2%). All 46 strains map to
iML1515 by exact gene symbol. Three self-pairs exceed 8-fold on the mean although the paper reports none; this is
reported, not adjusted. The label is a measurement; nothing is simulated.

Limits that matter. The strains are auxotrophs, so the set measures cross-feeding synergy between deletion mutants, not
growth of prototrophs or a Harcombe-type (E. coli / Salmonella / Methylobacterium) dependence. A static balanced-growth
model cannot reproduce the daily dynamics. The paper prints the M9 composition; `inputs.medium` lists those components
mapped to exchanges with concentrations as printed (g/l and mM have no flux equivalent, so `uptake_limits` is null and no
rate is invented). **Trace metals and the oxygen regime are not in the paper** (sealed microplates): they are listed
under `assumed_not_in_paper` and `oxygen`, and the root must freeze both before any result is read. The authors modelled
the system with iAF1260, not iML1515. The Harcombe COMETS models remain model-only.

## Secondary-metabolite production observations (item 12)

`prepare_mibig_production.py` turns the MIBiG 4.0 archive into one row per (cluster, compound): organism and NCBI
taxon, locus accession, the curators' evidence method, literature references, structure and database ids
(`mibig_production_observations.csv`; 5,443 compound records, **606 ready** from 196 taxa; evidence NMR 447, mass
spectrometry 118, MS/MS 80, authentic standard 27, fewer for the rest). Eligibility is fixed: entry active, NCBI taxon,
locus accession and at least one compound evidence record carrying a reference; the others are kept with the reason
(3,688 lack referenced compound evidence, 1,149 are not active). All 606 ready rows have MIBiG entry quality
`questionable`; quality is recorded, not used as a filter.

This source gives **positive observations only**. There are no negatives (absence of a record is not evidence of
non-production), no titres, no medium, and no mapping of compounds to model metabolites (a separate, registered step,
never guessed). It can score recall of characterised compounds or route existence for registered precursors, not
precision or medium dependence, and it is curated from literature that a method's training data may contain. A source
with measured negatives or conditions was not located.

## What the root can run

- **`manuscript_unregulated`**: `python benchmark/v3/run_metabolic_benchmark.py --cases <data>/derived/cases_manuscript_unregulated.jsonl --model <data>/raw/iML1515.xml --out <new dir>`. The runner is unchanged and accepts these records (`manifest` key present, model hash registered); the manifests were validated with the engine's `validate_manifest` / `build_menu` in tests, but no LP was run here.
- **`manuscript_regulatory`**: needs a baseline-condition closure in the engine (lower = upper = 0 on the sidecar's reactions, in both matched experiments, before the knockout). Until then it must not be run or reported.
- **Community pairs**: `wintermute2010_cases.jsonl` needs a pair-growth runner with a frozen trace-metal and oxygen choice.
- Report the three protocols side by side with the legacy run labelled as repository defaults. No condition may be tuned against the 61 labels.

## Exports (`<out>/derived/`)

| File | Content |
|---|---|
| `auxotroph_validation_labels.csv` | one row per strain-substrate: IDs and mappings with their mode, measured metrics, published and recomputed labels, final `label`, `binary_growth`, borderline and grid labels, raw-OD audit, `status` |
| `cases.jsonl` | one line per ready case: the row plus a `metabolic_experiments` manifest (matched wild-type `carbon_source:` control and the Keio `knockout:` under that carbon source; declared external exchanges; closed-by-default medium; explicit biomass and ATPM maintenance 6.86) |
| `model_metadata.json` | model file and hash, biomass, exchanges, medium and its provenance |
| `masking.json` | per-case masking requirements |
| `community_registry.json` | Harcombe models (hash, reaction counts, objective, notes) and the unavailable measured data |
| `quality_report.json` | all counts above, plus the `protocols` summary |
| `protocol_registry.json`, `regulatory_closures.json` | protocol parameters with provenance, discrepancies, per-substrate closures |
| `cases_manuscript_unregulated.jsonl`, `cases_manuscript_regulatory.jsonl` | case records per manuscript protocol (see above) |
| `wintermute2010_pairs.csv`, `wintermute2010_cases.jsonl`, `wintermute2010_quality.json` | measured co-culture pairs, label-free case inputs, data check |
| `mibig_production_observations.csv`, `mibig_production_quality.json` | positive production observations and their limits |

Scoring rule for a case: predicted growth is `max_growth >= growth_threshold` of the
knockout experiment (absolute growth on the substrate; the `growth_ratio` is
relative to glucose and is not used). The wild-type control on the same substrate
shows whether the model can grow on it at all. The generated manifests need the
root fix `822ab69` for GPR reading: on the earlier reader the `knockout:` item did
not resolve (checked, the carbon-source item and boundary validation were fine).
No LP benchmark was run here; the root runs them.

## Data dependencies still open

- Monk 2017 Keio growth and knockout profiles (item 8 breadth): needs a retrievable copy of Supplementary Data set 11
  (and 1 for the models) from a publisher or author repository that serves it without a challenge. Blocked, not
  bypassed.
- Harcombe 2014 measured co-culture growth and species ratios, with units, for the metB / Salmonella / Methylobacterium
  consortium (item 9): not retrievable. Wintermute 2010 gives measured pairwise growth in a different system. Measured
  exchanged metabolites and transporter-knockout phenotypes are still in no acquired source.
- Secondary-metabolite production with measured negatives, titres or media (item 12): not located. MIBiG is
  positive-only.
- Held-out organisms with measured growth across media beyond E. coli (item 12) are not covered.
- The manuscript regulatory protocol needs engine support for baseline-condition reaction closures before it can run.
