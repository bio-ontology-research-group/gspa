# Captured antiSMASH 8.0.4 records

`records.json` contains four unchanged records selected from an existing real
antiSMASH run on GCA_003598085.1: all three records containing regions and the
shortest record with no region. Top-level metadata is retained; only the records
list and JSON whitespace change. `provenance.json` pins the full source, selected
indices and resulting bytes. The source has 145 records and 4,915 CDS features.

The positive regions are T3PKS, betalactone and terpene-precursor. Their coordinate
intervals overlap 35, 22 and 13 CDS features respectively. These are predicted
region memberships, not 70 experimentally established biosynthetic enzymes.
The zero-region record establishes a parser/coverage case, not true biological
absence. The test genome deliberately renames the CDS proteins with `GSPA_` to
ensure antiSMASH locus tags are never used as protein identity joins.

Circular regions, altered identities, missing modules and boundary CDS tests
are explicitly constructed mechanism cases. The production invocation uses
FASTA and does not convey circular topology to antiSMASH. Source-version 8.0.4's
`Feature.start` and `Feature.end` explicitly use the first/last location parts,
and its serializer uses those properties for region areas; origin-crossing areas
therefore differ from Biopython's minimum/maximum location coordinates.

The external process in CLI verification is a declared captured-output replay.
No fresh antiSMASH inference or experimental BGC benchmark is claimed here.
