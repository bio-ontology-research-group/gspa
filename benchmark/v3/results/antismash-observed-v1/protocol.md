# antiSMASH region ingestion and protein mapping

Status: repair passes remote JVM gates and four CLI comparisons. Job 10439
provides the tested JAR; the completed CLI verification is job 10459.
This addresses input integrity for secondary-metabolism context. It does not
establish function recovery, BGC accuracy or completion of recommendations 3,
5, 7 or 12.

The old adapter read `cds_members` or `genes` inside region objects. Captured
antiSMASH 8.0.4 output instead supplies region coordinates and product classes
in `records[].areas[]`. The actual 145-record source has three regions, but the
old parser returns an empty map. Its hand-written fixture contained the invented
member lists and therefore failed to detect this loss. The adapter also returned
antiSMASH's re-called CDS identifiers, which need not match GSPA proteins, and
assigned every supposed member a fixed score of 0.8. Separately, enabling
`antismash` in configuration did not construct an adapter.

The repair reads native schema-4 regions, preserves product lists and source
identity, and projects their coordinates onto the original GSPA CDS using the
existing genomic-region helper. It never joins re-called locus tags or transfers
per-gene enzyme functions. Partial overlap remains distinct from containment.
Original regions survive independently of protein overlaps and are exported in
GFF3. All memberships are explicitly unscored: raw evidence is retained, while
calibrated probability and functional inference are unavailable. The GFF score
is `.`; an actual numeric zero from another source remains numeric.

This distinction follows the [upstream region interpretation](https://docs.antismash.secondarymetabolites.org/understanding_output/regions/):
region boundaries include rule-defined neighborhoods and are not experimentally
determined BGC borders. Proximity alone does not show that a protein catalyzes
secondary-product synthesis. No EC assignment or model product exchange is
inferred from a class name such as T3PKS.

## Configuration and coverage

```yaml
predictors:
  enable: [antismash]
  antismash:
    executablePath: /path/to/antismash
    database: /path/to/antismash/databases
    strictness: relaxed
    threads: 4
```

The CLI thread override reaches `--cpus`. The existing minimal, Prodigal-based
antiSMASH invocation is retained; output regions are mapped by coordinates to
GSPA's existing proteins. FASTA input does not convey circular topology to the
tool. The parser can represent origin-crossing regions in supported native JSON,
but that mechanism test is not a circular-genome execution benchmark.

An output record must match an exported contig and its nucleotide sequence.
Missing records, missing detection modules, absent input sequence and invalid
CDS coordinates remain explicitly partial. Wrong identities, malformed regions,
unsupported schemas and contradictory detection records fail. A valid scanned
record with no regions remains distinct from an unassessed record. A full
native-results document must be unambiguous; auxiliary metadata JSON does not
replace it. Region and record-coverage metadata are cleared on a failed rerun.

## Evidence and limitations

The versioned fixture retains three positive-region records and one genuine
zero-region record, unchanged from the source identified in its provenance.
The three regions overlap **70** source CDS features (35, 22 and 13). The tests
rename every protein to verify that coordinate projection preserves GSPA IDs.
Additional tests cover adjacent versus partially overlapping CDS, raw claim
round trips, unscored inference, malformed results and unavailable coverage.

The first remote job, 10437, reproduces the original loss and rejects the first
repair because private static helpers do not dispatch through the replay
subclasses. That failure is retained. Job 10439 passes **124 focused cases** and
the full **1,502-case JVM suite** (six optional skips, no failures or errors).
Its initial CLI checker incorrectly reads the appended GFF FASTA section as
feature rows and exits nonzero; that checker failure is also retained. Job 10459
uses the same tested JAR and the corrected checker and completes all four CLI
calls with exit zero.
The CLI experiment supplies all **145 records and 4,915 proteins** from the
captured output and uses an explicitly declared process replay. It does not run
antiSMASH inference again. The older packaged CLI is compared with repaired
annotation, annotation plus integration, and direct integration.

Both repaired annotation routes retain all 70 raw memberships, all three original
regions and the complete protein denominator. Their raw claims agree exactly;
annotation/direct integrated prediction files are byte-identical and contain no
functional calls from these unscored observations. The older CLI emits no
antiSMASH claims and registers no antiSMASH status. A separate checker uses native
`region` SeqFeature locations and qualifiers, rather than `areas`, to verify all
70 output memberships and their containment/partial-overlap labels independently.

The [evidence archive](../../benchmark/v3/results/antismash-observed-v1/README.md)
retains both failed attempts, the successful CLI run, source snapshots, reviews,
test XML and the compressed full native output. The validated development JAR is
`/data/hohndor/gspa-v3/antismash-observed-v1/jobs/10439/gspa-3.0.0.jar`, SHA-256
`b8b28f522732f8873ff2f899d2f6f4119c4ced2ca187f6c79f4d1a33b7fd0d7d`.

The Nextflow antiSMASH process remains a separate unresolved route: it publishes
a JSON artifact but does not yet feed this region evidence into shared
integration. Learned BGC contributions, exact product mapping, experimental
role/function truth and full benchmark acceptance remain open.
