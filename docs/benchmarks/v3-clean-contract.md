# CLEAN score semantics and integration audit

Status: controlled repair, full native execution, imports and independent archive audit verified.
No EC improvement, biological recovery, or whole-item acceptance is established.

## Defect and scope

The old neural sidecar interpreted every native CLEAN value as a confidence and
applied a lower confidence cutoff. CLEAN's native writer emits distances when
`gmm=None`, and GMM confidence values when a GMM is supplied. The inspected
installed driver and pinned upstream driver both supply a GMM. Thus the old
adapter's distance inversion is a real supported-mode defect, but is **not**
evidence that those particular GMM runs inverted their scores. The historical
full-matrix Price benchmark uses a separate ingestion route; its alleged gains
must be assessed using its own retained outputs.

Controlled before-repair examples show a distance of 0.01 discarded while 8.0627
is retained, malformed/empty output accepted as no calls, NaN accepted, distinct
small confidence values rounded together, and duplicate native rows accepted.
A shared `query` basename and original protein IDs also allowed stale or
concurrent native files/embedding caches to contaminate execution.

Primary source: CLEAN commit
[`f2bf2a4f497fa2cc87dac2a1bb314fee587c0a15`](https://github.com/tttianhao/CLEAN/tree/f2bf2a4f497fa2cc87dac2a1bb314fee587c0a15),
particularly `app/CLEAN_infer_fasta.py`, `app/src/CLEAN/infer.py`,
`app/src/CLEAN/utils.py`, and `app/src/CLEAN/evaluate.py`.
The upstream writer rounds to four decimal places. Preserving input precision
in GSPA cannot recover precision already lost upstream. CLEAN uses ESM-1b and
truncates sequences beyond 1,022 residues.

## Repair

CLEAN now requires an explicit `gmm_confidence` or `distance` mode checked against
the driver's literal GMM argument. Each execution uses an isolated filename and
opaque protein IDs, then restores original IDs. Every native row is validated;
unknown, duplicate, missing, malformed and nonfinite outputs fail. Provisional EC
identifiers such as `EC:1.1.1.n11` in the native training vocabulary are retained.
An explicitly empty input produces an empty receipt without invoking the model.

The typed TSV preserves native values, text, direction and availability. Its
receipt binds the original FASTA, full ID map, TSV hash, recorded native hash and coverage
counts. Confidence filtering defaults to zero. Distances are retained as
**unscored raw observations**, not transformed into invented probabilities. The
existing learned-candidate preparation excludes these unscored observations;
learning from native distance features remains future work. GMM confidence is
not asserted to be calibrated in GSPA.

`gspa clean-claims` and annotation share the same parser. Nextflow imports typed
CLEAN separately and keeps distance observations out of confidence ensembles and
evaluation. The generic TSV claim converter rejects CLEAN instead of dropping
its metadata. Direct ensemble loading skips explicitly unavailable scores.
The visualizer obtains CLEAN support from integrated results/provenance instead
of adding unchecked CLEAN TSV values as tool-only posterior scores.

This changes the CLEAN input contract: configurations must declare `scoreKind`
(or Nextflow `clean_score_kind`), and old four-column outputs need regeneration.
See [configuration](../reference/configuration.md) and
[Nextflow parameters](../reference/nextflow-params.md). Receipts validate internal
consistency and input binding; they are not proof of model execution or accuracy.
The importer does not read the native CSV and therefore cannot independently
verify its recorded hash. The separate native execution audit retains that CSV.
The native application directory must be writable for its staged inputs, outputs
and embeddings; UUID names prevent concurrent reuse, but an abrupt process kill
can leave owned temporary files behind.

## Validation record

Python regressions cover both native modes, precision, provisional ECs, complete
protein coverage, malformed output, native-file isolation, empty input, ensemble
exclusion and rejection by the generic claim converter. JVM tests cover
annotation/raw-claim preservation, configuration forwarding, exact input binding,
CLI import parity, replacement/tampering, and refusal to overwrite claims.
The controlled route experiment runs the actual Python sidecar, JVM annotation,
CLI import and Nextflow with deliberately synthetic native processes. It tests
both modes and propagation of missing/malformed/unknown-ID output failures.
These controlled processes do **not** run CLEAN's neural models.

Full native execution is separately registered on every protein in the authors'
`new` dataset (392) and `price` dataset (149): 541 records, no outcome filtering.
Reference EC labels are retained separately and never supplied to inference.
The scope is native execution, identity/format validation and reference agreement,
not independent context improvement or complete EC candidate coverage: maxsep
returns only selected candidates.

Initial native job 10195 stopped before inference because node002's loaded NVIDIA
kernel driver (580.159.03) differed from its userspace library (580.178.04).
A retained diagnostic reproduced the failure with and without `LD_LIBRARY_PATH`.
Node004 passed CUDA initialization, but attempt 10199 stopped before inference
because its local home cache lacked the ESM checkpoint. Neither attempt provides
prediction results. Subsequent execution uses a separately retained shared cache.
Job 10201 completes both full datasets successfully. The independent auditor
matches all **646 native selections** to typed outputs, including 24 zero GMM
confidence values and one provisional EC. Every registered native code/model
hash is unchanged after execution. This checks the recorded files; it is not a
claim that every transitive dependency was hashed before execution.

| Dataset | Proteins | Selected EC calls | Reference EC labels | Matches | Reference precision | Reference recall | Micro-F1 | Exact reference sets |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Authors' new | 392 | 465 | 503 | 241 | 0.5183 | 0.4791 | 0.4979 | 181 |
| Authors' Price | 149 | 181 | 152 | 71 | 0.3923 | 0.4671 | 0.4264 | 50 |

These are **fixed maxsep-selection reference agreement** results with no further
confidence cutoff, not F-max or a GSPA-versus-CLEAN comparison. An unmatched
prediction is absent from the supplied reference, not an experimentally verified
negative. Training/pretraining overlap and equivalence of the installed model
weights to a particular published benchmark were not established. The native
checkpoint hashes and complete predictions are retained; these numbers must not
be substituted for the historical dense-matrix Price results.

Validation passes 14 Python contract tests, 436 predictor and 88 CLI cases
(four optional skips), and 84 visualization checks (20 optional skips). The
unchanged core passes all 909 cases (two optional skips) in job 10200. The final
21 controlled invocations in job 10208 verify CLEAN and dbCAN annotation/import
parity, both CLEAN score modes, full Nextflow success/failure, and no promotion
of unscored distances. Native import job 10210 succeeds for all 541 proteins and 646 calls. The independent
[768-file archive auditor](../../benchmark/v3/results/clean-contract-v1/verify_archive.py)
verifies every file hash, exact primary-reference sequences, native/typed call
agreement, imported raw scores and complete metadata, all recorded gate counts,
and all 21 controlled invocations. [Full results and limitations](../../benchmark/v3/results/clean-contract-v1/report.json)
and [native reference agreement](../../benchmark/v3/results/clean-contract-v1/native-audit.json)
are retained with the failed attempts and read-only Sonnet reviews.

One import-audit assertion initially assumed `calibrated_prob` would be null for
GMM confidence. That was incorrect: the existing shared `ClaimExtractor` applies
its default heuristic transform. The repaired adapter preserves native GMM values
in `raw_score` and metadata; the derived `calibrated_prob` is **not learned or
validated calibration**. Distance observations remain unscored. The corrected
audit compares every raw score and all metadata and separately records the
finite derived probability range. The failed assertion and corrected run are
retained, not represented as a production prediction failure.

The controlled annotation failure case also confirms a broader CLI limitation:
`annotate` returns exit zero on predictor failure under its current best-effort
policy. Its predictor status manifest correctly records `FAILED`, and the raw
claims are empty. The route check distinguishes this from `SUCCESS_ZERO` rather
than interpreting process exit as predictor success. Nextflow native process
failure must still fail the workflow. A strict CLI exit policy remains outside
this adapter repair.
