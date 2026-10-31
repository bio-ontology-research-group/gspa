# Experimental phosphatase screen intake

Status: quarantined research intake; **no training or evaluation labels accepted**.
This investigates the independent negative/alternative evidence still missing
from acceptance items 4 and 7.

The primary source is Huang et al. (2015), *Panoramic view of a superfamily of
phosphatases through substrate profiling*, [DOI 10.1073/pnas.1423570112](https://doi.org/10.1073/pnas.1423570112),
[PMCID PMC4413258](https://pmc.ncbi.nlm.nih.gov/articles/PMC4413258/), PMID 25848029.
The full article was downloaded and read, and Europe PMC metadata independently
confirms its identity. Publisher supplement downloads returned HTTP 403; the PMC
workbook URL returned an HTML download challenge, not an Excel workbook. These
failed retrievals are retained and distinguished from missing publications.

Supplement copies and sequence files were retrieved from the public research
repository [enzyme-datasets](https://github.com/samgoldman97/enzyme-datasets/tree/627556e265e2a52d39753e684c05b550d53a9be4/data/raw/phosphatase),
pinned to commit `627556e265e2a52d39753e684c05b550d53a9be4`.
Equality to publisher-hosted supplement bytes is **unverified**. The upstream
processed data and chemical remapping code were not used to make labels.

## Raw inventory

The read-only inventory retains every assay cell from C7:N20, its sheet/cell
coordinates, exact substrate comment, lookup accession and numeric absorbance.
All values are finite. The workbook has **218 assay sheets**, **218 distinct
lookup accessions**, **36,406 substrate measurements** and **218 blank cells**.
Each assay has 167 substrate cells; substrate comments agree at every plate
position across sheets. There are no duplicated substrate texts within a sheet.

The paper reports **217 enzymes**, while metadata marks **216** entries screened.
Two assay sheets therefore need explicit reconciliation:

| EFI ID | Lookup accession | Metadata discrepancy |
|---|---|---|
| 502329 | P94592 | Present, but `Screened?` is `N` |
| 900135 | C2JBT4 | Absent from the metadata table |

All 218 sheets remain in the raw inventory. Thirteen have no values above the
paper's nominal 0.2 absorbance threshold. These descriptive counts do not choose
which entries belong in a benchmark or establish negative functions.

## Why this is not yet null or role truth

The paper measures phosphate release in a particular in-vitro assay, with
background correction and a stated detection threshold. Low activity can reflect
assay conditions or protein stability. In particular, an enzyme inactive at the
screen temperature has reported activity at higher temperature. A failed screen
therefore cannot establish universal absence of a GO function, an EC activity,
or a biologically correct candidate in a locus.

The inventory preserves raw names and values rather than interpreting the
mirror's fuzzy/manual compound mappings. One manual mapping changes a mannose
phosphate position; that mapping needs chemical review before use. Assay constructs,
organism/genome mappings, assay population, replicate/error handling and exact
substrate-to-reaction identities also require independent verification. No
inferred/propagated annotations in the third supplement are counted as measured
outcomes. These issues must be resolved before a split or label policy is frozen.

Working evidence is in `~/.cache/gspa-v3-data/had-substrate-evidence-v1/`:
retrieval receipts retain failures, successful source bytes and hashes;
`inventory_workbooks.py` records the extraction, `raw-observations.jsonl` retains
all 36,624 cells, and `inventory.json` records sample metadata and limitations.
This is a one-off intake inventory, not a validated training-data importer.

The [immutable intake archive](../../benchmark/v3/results/had-substrate-intake-v1/README.md)
preserves all retrieved artifacts, failure receipts and the full extraction.
Archive verification establishes byte preservation only.

## Primary sequence identity follow-up

A subsequent primary-database audit resolves all **218 accessions** and confirms
that every sequence in the mirror's accession/sequence table exactly matches its
primary sequence bytes. UniProtKB currently returns **192 active entries** and
**26 inactive entries**. The inactive records explicitly link UniParc identifiers;
all 26 linked archival sequences were retrieved and matched exactly. There are
zero missing sequences, sequence mismatches or retrieval errors.

The first audit assumed every returned UniProtKB record contained a sequence and
failed on an inactive entry. Its script and failure log are preserved. The
corrected audit preserves inactive status and follows only the UniParc identifier
provided by the primary record. In particular, the discrepant assay sheet C2JBT4
has an inactive primary entry and an exact linked archival sequence match.
This resolves sequence-byte identity; it does not resolve the assay-population
counts, assay constructs, organism/genome assignment or chemical mapping.
No database functional annotation was imported as an experimental label.

Receipts and source responses are in `sequence-verification-v1/` and
`sequence-resolution-v1/` under the working intake directory. The original intake
archive predates this follow-up and remains unchanged.

The [sequence identity archive](../../benchmark/v3/results/had-sequence-identity-v1/README.md)
preserves this follow-up and its failed first parser attempt.
