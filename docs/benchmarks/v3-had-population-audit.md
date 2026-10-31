# HAD assay population and identity audit

Status: **quarantined; no training or evaluation labels accepted**. This extends
the [initial intake](v3-had-substrate-intake.md) by comparing both metadata tables
with every assay sheet. All 36,624 archived cells match their original workbook
values, substrate comments and lookup accessions. The source copies remain pinned
to the same research mirror; equality to publisher-hosted bytes is unverified.

## Population reconciliation

Supplemental 2 marks 216 enzymes screened. Supplemental 3 supplies an additional
screening flag for EFI 900135, previously missing from the population metadata.
However, checking the entire tables exposes additional contradictions:

| EFI ID | Lookup accession | Supplemental 2 | Supplemental 3 | Assay sheet |
|---|---|---|---|---|
| 502329 | P94592 | N | N | Present |
| 900135 | C2JBT4 | Absent | Y | Present |
| 900190 | Q5ZS98 | N | Y | Absent |
| 508031 | D2QPE6 | Y | N | Present |

There are 218 distinct enzymes marked Y in at least one metadata table. Of these,
217 have assay sheets, matching the article's reported count numerically. This
intersection is a possible explanation, not confirmation of the intended study
population. EFI 502329 still has an assay despite two N flags; EFI 900190 has a Y
flag but no assay. The audit keeps all 218 measured sheets and never invents a
missing assay. Supplemental 3's propagated functional annotations are not used as
measured truth. These discrepancies require source clarification before choosing
a benchmark population or making a split.

## Assay identifiers versus lookup sequences

The previous sequence audit established that the mirror's sequences match their
lookup accessions. This follow-up checks a different link: the identifier written
in each assay sheet's D2 cell. All 194 numeric GI identifiers resolve in retained
NCBI protein records. Of these, **190 match the lookup sequence exactly and four
do not**. The other 24 sheets directly name the same UniProt accession as the
lookup table; those sequences also pass the previous primary-sequence check.

| EFI ID | Lookup accession | Assay GI → NCBI version | Assay / lookup length |
|---|---|---|---:|
| 501280 | Q8Z202 | 29144292 → NP_807634.1 | 252 / 252 |
| 508563 | Q5AM80 | 68466982 → XP_722406.1 | 271 / 265 |
| 900153 | P42962 | 16077472 → NP_388286.1 | 249 / 249 |
| 900165 | Q8NU30 | 19551370 → NP_599372.1 | 201 / 209 |

Equal lengths in the first and third rows conceal amino-acid substitutions.
These are discrepancies between two database identities linked by the workbook;
the audit cannot establish which sequence was expressed. No sequence is replaced
automatically and no assay label is transferred across a mismatch. The earlier
218/218 accession check remains correct, but was insufficient to verify the
assayed proteins. Exact matches also leave tags, engineered mutations and sample
handling unverified. All 218 entries remain in the report, including the four
unresolved cases.

## Measurement uncertainty

The [primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC4413258/) describes
duplicate endpoint phosphate-release assays, background correction, a nominal
0.2 absorbance threshold with approximately 0.03 variability, and loss of linear
response above 0.7. Some inactive proteins may need different temperatures or
have been unstable. These results measure activity under the tested conditions;
they do not establish physiological substrate preference or universal absence
of a function.

The complete workbook contains one numeric value per enzyme/substrate pair.
No separate replicate readings or per-pair error estimates are identifiable.
Its `Averages` sheet uses formulas across all 218 enzyme sheets, so it cannot
serve as an independent replicate. The only formula on an assay sheet counts
threshold hits; there are no extra numeric assay values outside C7:N20 from row 7.

| Descriptive signal stratum | Measurements |
|---|---:|
| Below 0.17 | 30,113 |
| Between 0.17 and 0.23, inclusive | 1,485 |
| Above 0.23 | 4,808 |
| Above the nominal 0.2 threshold | 5,440 |
| Above the 0.7 linear-response limit | 64 |

The first three rows partition all 36,406 substrate measurements; the final two
overlap them. The 218 blank cells are excluded. The 0.17–0.23 interval is a
descriptive sensitivity band based on the paper-wide variability, **not** a
confidence interval for individual observations. Saturated values cannot support
linear comparisons of catalytic strength. None of these strata becomes a
negative GO/EC label or a label that no valid gene exists in a locus.

## Consequences for the benchmark

This screen remains potentially useful for testing condition-specific substrate
activity after the population, assay construct and chemical identities are
resolved. It cannot currently supply the missing null-locus training truth.
A future preregistered evaluation should distinguish assay detection from
physiological function, retain uncertain measurements, report sensitivity to the
detection threshold, and separate related sequences across training and test data.
The complete population must be retained in the audit even if eligibility rules
later exclude entries. No model performance or improvement is claimed here.

The [reproducible audit](../../benchmark/v3/results/had-population-audit-v1/README.md)
retains all population checks, cell comparisons, source hashes and per-enzyme
signal counts.
