# Sonnet 5.5 review disposition

Reviewer found no blocking defect and confirmed that verified_batch output_hashes
and reusable are unchanged by extraction. The review was read-only and ran no tests.

Accepted: preserve frozen Price execution code; use only read-only verification
for the old dense run; document old-receipt invalidation and the large-ledger
runtime dependency; add parser code/library/backend cache invalidation checks;
reject an old verified receipt missing the helper binding; test that malformed
trailing metadata becomes a recorded generic-run failure.

The local/remote gate has 106 cases (105 pass, one optional JVM skip). The shared
parser tests measure Python retention and acceptance parity. The small generic
runner test checks dispatch and unchanged sidecar hashing, not process memory.
The real 2.68GB replay measures total process RSS through both entry points.

Additional actual-backend probes are now regression tests: invalid UTF-8 in an
ignored string is rejected by both paths; an escaped lone surrogate and 1e999 in
ignored values are accepted by both. They do not validate those values as usable
prediction scores; score validation is separate.

Limits retained explicitly: bundled native-library binding does not attest to
arbitrary system-linked libyajl or the whole Python runtime. Parser/source hashes
are not a full machine-image guarantee. Parsing and hashing are separate reads
and require stable finished artifacts. No claims of concurrent-mutation safety
or arbitrary-JSON fixed memory bounds were added.

Temporary directories in these test invocations use explicitly set TMPDIR
(local ~/tmp, remote /scratch/job). No frozen production output was rerun or
replaced. Packaging failure 10241 is retained; 10244 uses a fresh complete bundle.
