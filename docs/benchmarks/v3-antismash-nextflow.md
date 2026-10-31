# antiSMASH evidence in Nextflow

The main Nextflow workflow previously ran antiSMASH but never passed its results
to integration. This follow-up to the [native parser repair](v3-antismash-observed.md)
adds `antismash-claims`, using the same native schema-4 parser, sequence validation
and region projection as the annotation adapter. Nextflow imports and merges the
resulting unscored protein memberships. A receipt retains coverage and every
region, including regions without overlapping CDS.

The importer binds the native results to the supplied nucleotide sequences and
uses explicit GFF protein identities or Prodigal FASTA header IDs to join CDS.
For Prodigal headers, contig, start, end and strand must also agree. It rejects
ambiguous identities, duplicate FASTA IDs, unmapped CDS/proteins and invalid CDS
coordinates. It never pairs proteins by file order. It accepts a separate genome
FASTA and GFF; embedded FASTA is rejected. Protein translation accuracy is not
established by this import. Missing antiSMASH record coverage remains partial.

The antiSMASH process now fails on missing or multiple result JSON files, instead
of suppressing a failed rename. Enabling antiSMASH without its database fails.
Detection strictness is explicit (`relaxed` by default, matching the JVM adapter).
The raw claims remain `SECONDARY_METABOLITE` region observations with
`score_available=false` and no calibrated probability. This wiring does not train
a BGC contribution, infer an EC number or verify a biosynthetic product.

## Validation design

The controlled replay uses the complete captured antiSMASH 8.0.4 result already
identified in the native-parser report: 145 records, 4,915 proteins, three regions
and 70 memberships. The independent overlap oracle reads native `region`
SeqFeatures rather than the `areas` used by the parser. A second sample contains
the genuine zero-region record QZKU01000103.1 and its six CDS.

The test deliberately assigns native-shaped Prodigal GFF identifiers that differ
from FASTA protein IDs. Actual GSPA annotation, explicit import, and the complete
Nextflow workflow must preserve identical raw claims. Direct and Nextflow
integration must emit identical, header-only predictions. A native process that
exits zero without a result must fail the Nextflow workflow. The external
antiSMASH and gene-calling executables are declared captured-output replays,
**not fresh native tool runs or a biological accuracy benchmark**.

The Spock checks cover reordered proteins, exact region endpoints, partial
overlap, intergenic regions, missing records, genuine zero regions, and rejection
of identity/coordinate/sequence inconsistencies. A Sonnet 5.5 source review found
no correctness defects and identified repeated linear contig lookups in the new
importer; these now use an index. The shared region projector's scaling on very
large fragmented communities remains unevaluated.

## Validation status

Job 10480 retained a test-compilation failure: the new CLI test imported a JSON
library unavailable in that module. The production code compiled and all 27
existing adapter cases passed. The test now uses the module's existing Jackson
library. Job 10481 passes all **46 focused cases** and the full **1,521-case JVM
suite** (six optional skips, zero failures or errors). Its annotation/import raw
claims agree, but the first Nextflow launch fails because the container's default
home is read-only. Job 10482 fixes the cache binding, then exposes a missing
existing claims-parser script in the isolated test package. Both failed attempts
remain archived; neither is counted as a successful workflow validation.

Job **10483** includes that script with its checksum, reuses the same tested JAR,
and completes with exit zero. All 70 memberships agree with the independent
native-region oracle. Annotation/import/Nextflow raw claims are byte-identical;
direct/Nextflow predictions are byte-identical and contain no functional calls.
The zero-region sample remains empty and the intentionally missing native result
fails the workflow. No production code changed after the full JVM gate.

The [archive](../../benchmark/v3/results/antismash-nextflow-v1/README.md) retains
source and input hashes, scoped source snapshots, the native JSON, test XML,
review, failed attempts and the successful replay. The validated JAR is
`/data/hohndor/gspa-v3/antismash-nextflow-v1/jobs/10481/gspa-3.0.0.jar`, SHA-256
`069263f7b876222815928a637d02106b9c3b31a756edd1cf692432acc2e52303`.
These results establish ingestion and route parity; learned BGC evidence,
biological function recovery and all whole-item acceptance remain open.

## Follow-up correction: operon identities

The earlier synthetic input used GFF IDs `1_1` and `1_2`, while its protein IDs
were `chr_1` and `chr_2`. That demonstrates an incompatible-input failure but
**does not reproduce the output of the pinned Pyrodigal 3.7.1 caller**. Subsequent
real single/meta runs (job 10484) show all 4,359/4,291 GFF CDS IDs matching their
protein FASTA IDs, with all 3,867/3,789 operon member IDs intact. Numeric IDs appear
inside FASTA header metadata; the GFF uses full protein IDs. The earlier claim of
an unresolved production ID defect is withdrawn for this pinned version.

The historical archive preserves the original synthetic probe, which must not be
cited as evidence of a Pyrodigal 3.7.1 workflow defect. No normalization patch is
needed on that evidence. The [native Nextflow follow-up](v3-pyrodigal-operon-route.md)
now verifies both caller modes and context routing with explicit synthetic function
controls. Nextflow's actual operon algorithm is same-strand distance
grouping, not the three-predictor ensemble previously claimed in comments; those
comments are corrected. The BGC replay itself disabled functional priors and did
not establish operon scoring or biological validity.
