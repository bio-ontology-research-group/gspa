# Sonnet 5.5 review disposition

The read-only reviewer found no blocking semantic defect. It did not independently
check the real source counts. The archive verifier reconstructs all 8,472 entries
and all 1,309 legacy dispositions independently; the counts are not inferred from
reviewer approval or from six invented examples.

Accepted: require expected source and legacy SHA-256 values, check them before and
after parsing, bind the generator source, reject duplicate legacy JSON keys and
non-string or malformed target IDs, and test pinned end-to-end construction and
no-overwrite behavior. Eight unit tests pass. Real source assertions, including
counts and R/S stereochemistry, live in the mandatory portable archive verifier;
they are not optional tests skipped when a source download is absent.

Qualified suggestions: absolute paths are provenance, while portable verification
uses role-specific archived inputs and their hashes. No canonical-scoring consumer
is added yet: that consumer must use the complete registry, not only legacy alias
rows, and explicitly retain unknown/ambiguous cases and full coverage denominators.
Unknown predictor IDs are not evidence of absence, and must not silently disappear.
The preliminary flag is retained on each successful resolution. A failed output
write leaves an incomplete immutable attempt, not a successful report to reuse.

The audit does not alter raw exports or historical projection statistics. It does
not establish whether any observed prediction or benchmark score changes. That
requires the complete scoring comparison after the registered inference arms finish.
