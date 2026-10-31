# Review dispositions

The first Sonnet review identified the real missing member-list/identity/score/registration problems. Its suggestion that native records lack areas is contradicted by the captured 8.0.4 output and pinned serializer; areas are used directly. No need to re-call GSPA genes or change antiSMASH gene calling to establish coordinate membership, and no per-gene function annotations are transferred.

The implementation review found no concrete parser/configuration bug. It explicitly did not execute tests. The first remote gate then exposed a Groovy private-static helper dispatch failure through replay subclasses; helper visibility is corrected and the failed run is retained.

Circular representation is verified against pinned upstream feature.py: start/end select the first/last location parts, not Biopython min/max, and serialiser.py uses those properties. The mechanism test is synthetic and FASTA invocation does not convey circular topology. No circular-genome execution claim is made.

Missing detection identity/structure, blank original IDs, invalid products, empty record lists and boundary overlaps received tests. Malformed extra JSON remains a fail-closed contract choice. Large-document streaming is outside this bounded repair; the 7.5 MB full source is used for CLI validation. Circular segment copies retain common region identity and stay unscored; no unique functional evidence count is derived from them.

Nextflow integration is explicitly outstanding; no route-completion claim is made. A region with no CDS still provides a valid genomic observation and is retained in metadata/GFF, so zero CDS overlap alone is not classified as a tool failure.
