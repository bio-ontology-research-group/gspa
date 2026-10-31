# Sonnet 5.5 read-only review dispositions

- Frozen source: the reviewer saw the local editable worktree, but execution uses
  a separate remote source snapshot, checked from source.sha256. No local edits
  propagate to it. Failed v1 and final v2 sources remain separate. No requirement
  to commit unrelated shared work or add another copying layer.
- Accepted package pin check. Preparation now requires the same eight package
  versions as successful native job 10201; execution checks them again. Native
  code and weights are content-bound, but this is not whole-system attestation.
- Genome identity: production integration already scopes claim files by genome;
  it must not flatten files across organisms. Added genome_id in claim metadata
  and the native audit for an explicit additional record.
- Scratch: the code already used parents=True and the wrapper passes a unique
  /scratch/$SLURM_JOB_ID. No shared scratch pathname or source edit is needed.
- Accepted alphabet preflight using the installed ESM-1b alphabet. Reject unknown
  residues rather than silently changing sequence bytes. No target is filtered.
- Accepted explicit plan-to-receipt driver hash and nonempty-input checks.
  Added deterministic opaque-ID/input-order validation to detect swapped maps.
- Accepted per-genome sequence-overlap counts. This is reference-membership
  evidence, not a complete homolog/pretraining audit. No reference EC annotations
  are extracted for case selection, scores, or training labels.
- Added native-ID swap, wrong driver, and existing-claims regressions. Eight
  new observer tests pass locally. The earlier five-case version passed remotely.
- Initial remote attempt 10324 failed before environment creation or inference:
  the test package omitted build_ensemble_preds.py and sidecar_to_claims.py.
  The corrected seven-file package preserves the required repository layout.

The review ran no models and changed no files. Live native prediction and all
eighteen genome outputs still require their execution receipts and audits.
