# Experimental negative annotation audit

The June 2026 GOA snapshot contains useful negative assertions, but it does not
provide a complete null/role benchmark for the existing development cohort.
The audit and candidate intersection are complete; source review is selective
and explicitly recorded. No negative assertion is automatically a training label.

## Source and complete scan

`benchmark/v3/audit_gaf_negatives.py` streamed the entire local
`benchmark/neural/cafa6_recon/goa_uniprot_all.gaf.gz` file, checked gzip integrity,
and verified compressed-source, code, evidence-policy and identity-mapping hashes
before/after processing. The source declares GAF 2.2, generation on
2026-06-17 14:49, and GO release 2026-06-15. It contains **667,198,446 lines**
and decompresses to **126,032,427,777 bytes**. Its compressed size is
11,664,243,116 bytes and SHA-256 is
`f315375b07946a0649142b2f4de2e15e282316989677a04e7a561203186dd2ff`.

The scan found **10,526 negative assertions**, including **4,117** with an
experimental GO evidence code. These are annotation rows, not independent
experiments or proteins. All 17 source fields, relation, references, evidence,
isoform, annotation extension, interacting taxon and original row identity are
preserved. The reader validates negative rows; it is not a validator of every
positive annotation in the source.

[GAF 2.2](https://geneontology.org/docs/go-annotation-file-gaf-format-2.2/)
represents negation with `NOT|relation` and keeps contextual extensions and
gene-product forms in separate columns. The
[GO annotation guide](https://geneontology.org/docs/go-annotations/)
distinguishes `enables`, `contributes_to`, process relations and localization
relations. Missing annotations are unknown; negative annotations propagate
downwards rather than to ancestors. This inventory performs **no propagation**.
The experimental-code set follows the
[GO evidence guide](https://geneontology.org/docs/guide-go-evidence-codes/);
curator inference, similarity and automatic evidence are not experimental truth.

## Intersection with production candidates

Thirteen negative rows join accessions in the existing sequence-matching table.
One has similarity evidence (`ISS`) and is excluded. The other twelve require
reference and sequence-version review before use. An accession join alone does
not prove identity between annotation and protein snapshots.

`benchmark/v3/assess_gaf_candidate_coverage.py` verified the hashes/counts of all
18 production feature exports, containing **484,202 candidate rows**. Eight of
the twelve experimental annotations have no exact GO/protein candidate. Four
have exact candidates:

| Protein | Negated term | Reference | Exported competitors |
|---|---|---|---:|
| Q97VW8 / AAK42622.1, ThrRS-cat | GO:0002161, aminoacyl-tRNA deacylase activity | PMID:15240874 | 10 |
| Q03513 / ABB57453.1, CcmM | GO:0004089, carbonate dehydratase activity | DOI:10.1139/b05-057 | 5 |
| Q980D1 / AAK40712.1, ThrRS-ed | GO:0004829, threonine-tRNA ligase activity | PMID:15240874 | 2 |
| Q980D1 / AAK40712.1, ThrRS-ed | GO:0006435, threonyl-tRNA aminoacylation | PMID:15240874 | 2 |

Competitors are all exported proteins for the same GO term, role and genome,
including the listed protein. Counts do not establish truth. This intersection
uses exact identifiers without alternate-ID normalization or descendant closure;
it does not establish that the source contains no other potentially useful labels.

## Source review and exclusions

The full text of Korencic et al.,
[A freestanding proofreading domain is required for protein synthesis quality control in Archaea](https://pmc.ncbi.nlm.nih.gov/articles/PMC478561/),
identifies SSO2486 (ThrRS-cat) and SSO0384 (ThrRS-ed), cloned from ATCC 35092D.
The threonylation experiments distinguish the catalytic protein from the editing
protein. The editing assay tests Ser-tRNAThr; it does not exhaust all possible
misacylated tRNA substrates. The GO:0002161 competition also contains ten
proteins, so this pair alone cannot certify its complete truth. The process
annotation GO:0006435 is kept separate: absence of aminoacylation activity does
not establish absence of every form of involvement in the broader process.

The Q03513 reference has not been verified against its full text and remains
unaccepted. Its DOI resolves to *Cyanobacterial carbonic anhydrases*, a review;
the reported direct-assay evidence must be traced before use.

A separate source discrepancy illustrates why code `IDA` is insufficient by
itself. Q6LZY8 has `NOT|enables GO:0004662` citing PMID:24291792. The pinned GO
ontology defines that term as **CAAX-protein geranylgeranyltransferase activity**.
The retrieved full text of
[Mechanism of farnesylated CAAX protein processing by the intramembrane protease Rce1](https://doi.org/10.1038/nature12754)
reports that the archaeal protease did not cleave geranylgeranylated peptides.
That is a different reaction from transferring the lipid group. This review
does not establish a direct negative transferase assay; the row remains
unaccepted. It has no exact candidate in the present export anyway.

The full text of
[Crenarchaeal arginine decarboxylase evolved from an S-adenosylmethionine decarboxylase enzyme](https://pmc.ncbi.nlm.nih.gov/articles/PMC2533785/)
supports the enzyme-specific negative for Q9UWU1, but there is no exact
GO:0004014 candidate in this export. Its separate negative process annotation
is not transferred to a generic process-association label. Other unmatched
assertions remain a source-review queue; they were not silently discarded or
converted into training examples.

## First reviewed complete training competition

`benchmark/v3/prepare_reviewed_thrs_locus.py` creates an augmented copy of the
existing cohort. Both versioned protein accessions and full sequences match the
acquired GCA_000007005.1 proteome. UniProt reports sequence version 1, last
updated 2001-10-01, for both proteins, predating the GAF snapshot. Sequence
hashes, source record hashes and the archived full-text hash are retained.

For **GO:0004829 / UNSPECIFIED**, Q97VW8 is positive and Q980D1 is explicitly
negative, both supported by PMID:15240874. These are the entire two-member
exported competition. The standard builder verifies this before/after label
joining and family masking; the production trainer's certificate validator
accepts the resulting membership and provenance. This establishes completeness
of the exported competition, not experimental characterization of every protein
in the genome or absence of activity under all conceivable assay conditions.

All 4,262 existing candidate identities, feature values, binary labels and splits
remain unchanged. One previously unconfirmed zero now has explicit experimental
negative support. The two competition members retain their original **training**
assignment. There is **one complete training group**, zero calibration/test
groups, and no all-negative group. No model was fitted or promoted from this
panel; it cannot establish null calibration or independent performance.

```bash
"$HOME/tmp/gspa-v3-testenv/bin/python" benchmark/v3/prepare_reviewed_thrs_locus.py \
  --original-manifest "$HOME/.cache/gspa-v3-data/locus-truth-v2-final/table/training_table_manifest.json" \
  --provenance benchmark/v3/results/training-cohort/provenance_manifest.json \
  --inventory "$HOME/.cache/gspa-v3-data/gaf-negative-audit-v1/inventory.json" \
  --references "$HOME/.cache/gspa-v3-data/gaf-reference-review-v1" \
  --out "$HOME/tmp/reviewed-thrs-locus"
```

The original cohort/model artifacts remain historical evidence. The new table
is in `~/.cache/gspa-v3-data/reviewed-thrs-locus-v1/`; its compact provenance,
labels, policy, locus declaration and two certified rows are archived in
`benchmark/v3/results/reviewed-thrs-locus-v1/`. Two additional tests reject
sequence, accession, taxon, version-date and truncated-record mismatches.
Actual null outcomes, further independent competitions and role-distribution
coverage are still needed for items 4 and 7.

## Reproduction and evidence

Use new output directories:

```bash
export TMPDIR="$HOME/tmp"
"$HOME/tmp/gspa-v3-testenv/bin/python" benchmark/v3/audit_gaf_negatives.py \
  --source benchmark/neural/cafa6_recon/goa_uniprot_all.gaf.gz \
  --matches benchmark/v3/results/training-cohort/swissprot_matches.tsv \
  --out "$HOME/tmp/gaf-negative-audit"
"$HOME/tmp/gspa-v3-testenv/bin/python" benchmark/v3/assess_gaf_candidate_coverage.py \
  --inventory "$HOME/tmp/gaf-negative-audit/inventory.json" \
  --features "$HOME/.cache/gspa-v3-data/training-cohort/features" \
  --out "$HOME/tmp/gaf-negative-audit/coverage"
```

The full inventory remains in
`~/.cache/gspa-v3-data/gaf-negative-audit-v1/`; source text and retrieval records
are in `~/.cache/gspa-v3-data/gaf-reference-review-v1/`. Compact inventory,
cohort assertions, candidate-intersection results, source dispositions and
verification logs are archived under
`benchmark/v3/results/gaf-negative-audit-v1/`.

Ten focused inventory/coverage tests cover polarity, malformed rows, contextual
exclusions, CRC errors, changed files, incomplete inventories and exact matching.
The final benchmark suite passed **504 tests**, with three optional skips,
including eight live Nextflow tests and the two additional sequence-identity
tests. These checks validate the audit machinery and preparation; they do not
establish accuracy of a learned model.
