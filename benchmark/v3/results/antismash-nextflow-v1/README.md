# antiSMASH import and Nextflow replay evidence

See the [report](../../../../docs/benchmarks/v3-antismash-nextflow.md) for scope and
limitations. This is captured-output route verification, not fresh antiSMASH
inference, an EC accuracy benchmark or evidence that a biological function was
recovered.

- Job 10480: production code and 27 adapter cases pass; new CLI test compilation
  fails on an unavailable test JSON library. The failing test source is retained.
- Job 10481: corrected tests; 46 focused and 1,521 full JVM cases pass (six optional
  skips). Annotation/import parity passes, Nextflow launcher fails on a read-only
  default home directory. The successful JVM gate and failed replay are separate.
- Job 10482: explicit Nextflow cache binding; workflow fails because the isolated
  package omitted the existing claims-conversion script. Failure retained.
- Job 10483: same tested JAR, explicit cache and hash-bound parser overlay. All
  six top-level command checks pass their expected outcomes: the old command is
  unavailable, repaired import and annotation agree, complete Nextflow and direct
  integration agree, and missing native output fails. Positive and zero-region
  samples preserve 70 and zero memberships respectively, without functional calls.

Run `python3 verify_archive.py` here. It checks all archived hashes, the nine
scoped source snapshots against the full 828-file JVM build manifest, both replay
overlays, the 1,521-case test gate, raw byte parity and all 70 memberships against
native region SeqFeatures (independent of the parser's `areas`). It also checks
complete zero-region coverage and the expected success/failure outcomes.

`source/` is the source-v2 scoped snapshot, including the original replay checker.
The root-level checker is the corrected replay overlay, accompanied by its hash
and execution scripts. The full source packages, build caches and JAR remain on
unimatrix; source manifests and package hashes are retained here. Omitting those
large artifacts means this archive verifies retained evidence, not an independent
rebuild of the JAR. The replay's large generated FASTA/GFF/native copies and
Nextflow work/cache directories are omitted; the scripts regenerate them from the
compressed native JSON, and remote receipts bind them to the actual run.

The native JSON and provenance are the same as in
[the prior observed-output archive](../antismash-observed-v1/README.md): antiSMASH
8.0.4 on assembly GCA_003598085.1, native SHA-256
`423d21fab26bfcfacb0dafd4a21f10a7172c659ee0dc3bd36c2b4cfcbc2ec9a6`.
The tested JAR SHA-256 is
`069263f7b876222815928a637d02106b9c3b31a756edd1cf692432acc2e52303`.

`operon-id-finding.json` records a separate synthetic native-shaped ID mismatch
in the existing operon script. It remains unresolved. This replay disables
functional priors and does not certify the operon route.
