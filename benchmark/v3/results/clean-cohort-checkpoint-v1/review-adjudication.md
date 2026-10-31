# Sonnet 5.5 collection-audit review dispositions

Accepted: exact typed metadata keys, explicit boolean/number treatment, canonical
claim comparison, and rejection of extra/missing TSV columns. Added an attack
case that changes both TSV metadata and claims and reseals the receipt/inventory.

Accepted: independent reconstruction of all four predictor-contract digests,
parent-cohort population equality, staged-file identity, and AST verification of
the native driver's explicit GMM call. The auditor does not import clean_contract:
importing its validator, as suggested, would weaken the independent check.
Added valid/damaged contract, changed parent population and driver-mode tests.

Accepted: exact native command and manifest checks against this registered Slurm
wrapper; the scratch path is /scratch/<job>/clean-app. Added a resealed command
mutation test. This is an acquisition-specific auditor, not a generic portable
native-job validator. Serialize the report before opening its exclusive output.

Retained conservative treatment of unfinished attempts. A missing marker does
not prove that a process has stopped. No attempt is currently killed or retried.
If a future OOM/timeout leaves an unfinished directory, require explicit terminal
scheduler evidence in a versioned finalization record before extending the
policy; do not silently accept incomplete work as a failure. Multiple successful
attempts also remain an error, preventing result-dependent selection. Two tests
for these cases were added while the review was running.

Qualified: the final audit runs after successful completion of the full array;
the acquisition snapshot and attempts are not edited in place. Repeated hashes
do not provide adversarial filesystem attestation. The report explicitly treats
reference-overlap counts as the hash-bound acquisition report, rather than an
independent homology audit. Failed-attempt markers remain reported separately
and cannot substitute for any successful genome. Neither limitation justifies
changing frozen acquisition jobs or claiming extra biological validation.

The review was read-only. Final local validation: 17 collection-audit tests and
8 acquisition tests pass (25 total). Full remote collection verification remains
pending; these tests are integrity fixtures, not biological benchmark results.
