# antiSMASH ingestion evidence

The repair preserves native region membership without inventing enzyme function
or confidence. The native 8.0.4/schema-4 source has 145 records, 4,915 CDS features,
three regions and 70 coordinate memberships. No fresh antiSMASH inference or
experimental function-accuracy result is claimed. CLI execution uses a declared
captured-output replay, validated against the full input sequences.

- Job 10437 reproduces old-parser loss and records the initial private-static
  helper dispatch failure. Only its actually executed predictor XML is retained.
- Job 10439 passes 124 focused cases and all 1,502 JVM cases (six optional skips).
  Its CLI checker fails after parsing an appended FASTA section as feature rows;
  that failure is preserved. The tested JAR remains unchanged.
- Job 10459 corrects only that checker, verifies the prior test gate/JAR hashes,
  and passes four real GSPA CLI calls plus an independent native-feature audit.
  Both repaired annotation routes retain all 70 memberships and three GFF regions;
  raw evidence agrees, and direct/annotation integration outputs are identical.
  The unscored observations produce no functional calls.

`build-source-v1.tar.gz` and `build-source-v2.tar.gz` contain exact source snapshots
(801 and 802 files), checked against their manifests. The final CLI checker is
stored separately because it changed after the tested Java build. Upstream
serializer/location source and its AGPL licence support the schema interpretation.
The four-record observed fixture is independently checked against the full source.

The full native output is losslessly compressed here. Generated FASTA/GFF files
from CLI runs remain in the remote job, with their hashes retained in
`job-10459/cli/verification.json`; they can be regenerated from the pinned source
and checker. Small actual CLI artifacts, raw observations and reports are retained.
The portable verifier checks retained files and independently derives membership
from native region SeqFeatures rather than the adapter's areas objects. It does
not rerun Java or claim to reverify omitted remote output bytes.

Remote root: `/data/hohndor/gspa-v3/antismash-observed-v1`.
Tested JAR: `jobs/10439/gspa-3.0.0.jar` (hash in job-10439/jar.sha256).
The expired scheduler entry for 10437 is represented by its retained execution,
logs and exit-code receipt; the scheduler snapshot retains 10439 and 10459.

Run `python3 verify_archive.py` in this directory. It uses a private directory
under `~/tmp` for the decompressed source and removes it afterwards.

All twelve whole acceptance requirements remain open. The Nextflow BGC route,
learned contribution, exact product mapping and independent function recovery
remain outstanding. No release or commit was made.
