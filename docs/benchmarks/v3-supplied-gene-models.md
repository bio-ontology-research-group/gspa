# Supplied GFF3/protein identity and neighborhood preservation

The actual `gspa annotate --genome ... --gff3 ... --proteins ...` route lost all
positions for native Prodigal output. Its GFF CDS IDs are numeric (`1_1`), whereas
FASTA protein IDs include the contig name; the connecting numeric ID occurs in
header metadata. The previous join discarded that metadata. This prevented
neighborhood context even when complete coordinates were provided.

## Observed failure, not an accuracy benchmark

Job **10523** used the validated 10520 JAR
(`1ae842e5905eecfff4ecbcb07bfa21d320ff9de0797273e4aeafd2abaccbdb5d`)
and the complete retained native outputs from gene-caller jobs 10507/10515.
All four annotation commands exited zero. Both tools used the same 145-contig
assembly; these are genuine caller outputs, not experimental gene-function truth.

| Supplied input | Proteins | Correct positions before repair | Neighborhood groups / members |
|---|---:|---:|---:|
| Pyrodigal 3.7.1 single | 4,359 | 4,359 | 819 / 3,867 |
| Pyrodigal 3.7.1 meta | 4,291 | 4,291 | 812 / 3,789 |
| Prodigal 2.6.3 single | 4,359 | 0 | 0 / 0 |
| Prodigal 2.6.3 meta | 4,286 | 0 | 0 / 0 |

The assembly SHA256 is
`76a0087d8d90585696bd21d4c01b0e1f6699c72cce2a9a4682898c0255154ed5`.
All native protein IDs were retained, so a protein-count-only check would have
missed the defect. The earlier suspected Pyrodigal Nextflow identity problem was
not reproduced; this is a distinct, verified problem in the supplied-input route.

## Repair and boundaries

`AnnotationPipeline` now retains full protein headers for identity matching and
shares the native alias/contig/coordinate/strand resolver with `GeneModelReader`.
Native coordinate evidence constrains all aliases before choosing a unique match;
otherwise numeric contig names can collide with native numeric IDs. The strict
complete-model reader used by gene callers and antiSMASH imports remains strict.
Supplied annotation can still be partial, use GFF-embedded sequences, or omit
nucleotide FASTA when protein sequences are supplied.

Ordinary declared multi-segment proteins retain the representative-CDS convention:
leftmost start, then shortest end for a deterministic tie. Exact duplicate rows,
ambiguous cross-contig/strand joins and competing proteins on any segment are
unplaced. This does not add transcript assembly or full multi-exon spans. Shared
protein IDs can still falsely group tandem loci if the input itself assigns them
one identity; annotation identities are not biological proof. Invalid input such
as duplicate FASTA identifiers is rejected; arbitrary malformed numeric headers
are not promised graceful recovery. Unmatched sequences remain available for
sequence-based prediction. Synthetic contig names avoid existing assembly names.

## Validation and experiment records

- Local working root: `/home/leechuck/tmp/gspa-supplied-gff-v1`.
- Remote root: `/data/hohndor/gspa-v3/supplied-gff-v1`.
- 10523: complete observed baseline, exit zero.
- 10524: intermediate v1; 13 new regressions reproduce 11 failures before repair,
  63 focused cases and the full intermediate gate pass after repair.
- 10525: v2 cancelled while verified pending, before execution; retained receipt.
- 10526: final v3 gate passes **67 focused cases and 1,569 full cases**
  (1,563 executed, six optional skips, no failures/errors). All 17 new cases run;
  the exact old source fails 14 of them. Exit zero.
- 10527: all **ten actual CLI cases pass**, exit zero. Prodigal single/meta
  now place all **4,359 / 4,286** proteins and yield **819 / 813** groups
  containing **3,867 / 3,790** proteins. Pyrodigal exports remain byte-identical
  across all seven checked output files in both modes. Numeric cases place all
  294 / 291 proteins. The single conflicting-header control retains its sequence
  unplaced; the other 4,358 proteins keep correct positions.

Final source v3 contains 845 hash-bound files. The validated JAR is
`/data/hohndor/gspa-v3/supplied-gff-v1/jobs/10526/gspa-3.0.0.jar`, SHA256
`51adf1c96de8466e4373de88fa4d2f295f482fa4cc0de3bdf9c2e7a7f17124f1`.
It supersedes 10520 for new development; queued CLEAN/context and Price jobs keep
their frozen versions.

The final CLI checker verifies every native sequence and coordinate, compares
complete neighborhood memberships to an independent adjacency calculation, and
requires byte-identical Pyrodigal protein/annotation/GFF and neighborhood exports.
Additional cases cover native numeric contigs, gzip, no nucleotide FASTA, reversed
protein order, and one deliberately conflicting header coordinate. That injected
case must retain the protein but withhold genomic context. These are identity and
route tests, not measured operon accuracy or EC recovery.

Two read-only Sonnet 5.5 reviews were used. The repair adopts deterministic segment
selection and group-wide collision checks. It does not adopt alias-first native
matching or coordinate-only search across contigs, which could choose a wrong
identity. The claimed leading-whitespace parser disagreement was rejected after
checking that `FastaReader` trims each line before identifying headers. Complete
native CLI validation supplies the coverage absent from synthetic unit fixtures.

Archive: `benchmark/v3/results/supplied-gene-models-v1`, **1416 files plus
`archive.sha256`**. The portable independent verifier passes all native coordinate,
sequence, neighborhood, baseline parity and test-result checks. It verifies 1,412
retained/decompressed bindings, all 845 frozen sources and 58 runtime input bindings.
Two JARs and Java are omitted from the archive and explicitly identified by hash;
the build/runtime packages remain at the recorded remote paths. See
`local-verification.json` and run `python3 verify_archive.py` inside the archive.

All twelve whole acceptance requirements remain open. This fixes an upstream
loss of context; whether preserved context improves function predictions must
still be demonstrated by the registered biological comparisons.
