# Sonnet 5.5 review adjudication

Review was read-only, based on source and fixtures; it did not execute dbCAN.

- Structured metadata: confirmed. Keep existing in-memory lists/maps and add explicit `rawFieldsJson` and `supportingMethodsJson` copies. The existing shared claims writer stringifies nested values; the new JSON copies are parseable after export without changing other predictors. Tests recover subfamily identifiers and domain intervals independently of the expected writer.
- Protein-level tool count: confirmed potential ambiguity. Rename this newly introduced field `proteinToolCount`; do not treat it as per-family consensus or probability.
- Version: implement the pinned version-5 `run_dbcan version` command. Legacy version remains explicitly unknown. Interface selection is not verification of the installed release. Native import receipts bind input/output bytes and explicitly do not attest execution/version/database identity. Actual database and executable identity remain a limitation of this controlled route experiment.
- Obsolete class documentation: corrected for both interfaces. The executable in upstream 5.2.8 is `run_dbcan`, as established by its console-script entry point.
- Integration evidence: add an assertion that standalone integration loads all four observations; also retain raw claims and compare annotation/import/Nextflow bytes. This does not establish that the web interface presents unscored observations usefully; that remains open.
- Suspected CPU oversubscription: not substantiated for the pinned v4 source. Its core DIAMOND and HMMER calls are sequential; the `split_uniInput` implementation calls one hmmsearch using `hmm_cpu` (despite stale parallelization prose). Do not divide the configured limit based on speculation. This does not certify resource use of all releases or optional SignalP.
- Implicit database directory: require an explicit directory for both supported interfaces.
- Unknown protein IDs: validate every overview row against the provided FASTA, including zero-hit rows, in native execution and import.
- Dual FASTA readers: assert their identifier sets agree; retain the separate duplicate-header check because the shared FASTA reader overwrites duplicates.
- Partial import publication: stage claims and receipt in the destination directory; publish only after both writes succeed, and remove our published claims on receipt-publication failure. This is ordinary exception cleanup, not a transaction across process termination or filesystem failure.
- Real token/header coverage: still not certified. Controlled fixtures and upstream source inspection are not a genuine full dbCAN search. Pinned v4 code normalizes GT2 profile variants before overview generation; do not weaken token validation based on unverified recollection about other releases.
- Concatenation: all GSPA writers used here terminate JSONL records with newlines. The route experiment must verify final merged bytes; broader third-party sidecar contracts are outside this dbCAN repair.
