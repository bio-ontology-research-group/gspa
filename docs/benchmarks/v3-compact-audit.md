# Compact storage for complete prediction-change ledgers

The audit previously retained a seven-entry map for every changed candidate and
built a union set of all before/after keys. The largest complete CLEAN matrix has
45,290,880 candidates. Keeping its initial change records as maps adds substantial
memory use even though the ledger only needs identities and before/after scores.

`AuditChanges` now captures shared key references, two primitive score arrays and
presence bits. It creates each record when a reader or Jackson requests it. The
serialized fields, order, tolerance and scores remain the compatibility boundary;
no candidates are removed and no probabilities or evidence rules are changed.
The before-key order followed by after-only keys preserves the previous union
order. Captured values survive later changes to the source maps. Snapshot lists
and records now reject mutation; the consumer review found no code that mutates
them. A malformed changed key still fails immediately with a null-pointer
exception, but its message is now explicit rather than Groovy's incidental text.

The source-level improvement is smaller retained audit storage, not a claim about
enzyme prediction. Source score maps, annotations and provenance still consume
memory. Building a snapshot scans the maps twice, and reading records still
allocates temporary maps and parsed identity strings. Whole-application resource
measurements are required before claiming an overall memory or runtime gain.

## Observed baseline

Job 10230 runs the frozen production JAR on the largest CLEAN matrix, without
priors. Its JVM exits successfully after 35:39.63 with 180,227,072 KiB peak RSS.
The [subsequent complete verification](v3-price-local-execution.md) preserves all
45,290,880 candidate identities, with 5,087.89 seconds total diagnostic time.
This diagnostic lacks the default-prior arm and cannot complete a registered
two-arm collection cell.

A runtime recording and thread dump locate an extended pause in
`FileDescriptor.close0`, called while closing the integrated TSV. The worker was
also observed waiting on page writeback. The job reached its 192 GiB cgroup limit
with approximately 172 GiB of anonymous memory and 20 GiB of file cache; the
cgroup reported no OOM kill. The output eventually closed successfully. These
observations support memory and filesystem pressure, not a JVM crash or a proof
that the audit alone caused the delay. The sampled allocation profile highlights
maps, arrays and boxed numbers, but allocation weights include temporary garbage
and do not measure retained object sizes.

The separate original-JAR DPZ job 10235 completes both arms with all 9,790,188
candidates retained. Its default-prior TSV is byte-identical to the earlier
indexed implementation (`64701c3dd55cd6dc753139df7b93eb10e49f1efad361ea943f19f42f7aa6136f`).
This establishes a complete baseline for the compact-audit comparison.

## Verification and remaining scale check

Job 10273 retains the first compilation failure, caused by byte-array compound
assignment under Groovy static compilation. Explicit result casts repair it.
Job 10275 passes the 11 focused tests: seven snapshot tests and four existing
audit integration tests. These compare the complete pretty-printed JSON with the
previous implementation, including a seeded 10,000-candidate panel, additions,
removals, reversed key order, mutation of source maps, nulls, NaN in either
position, infinities, signed zero, the tolerance boundary and protein identifiers
containing pipes. A read-only Sonnet review and its adjudication are retained.

Job 10275 finishes successfully: the complete JVM gate runs 1,443 cases with no
failures or errors and six optional skips, and builds the JAR. Both TopEC pilot
arms preserve all 3,481 candidates, every TSV byte, and the complete ledger except
stage elapsed times. The portable verifier independently compares the retained
old/new small outputs rather than relying only on the comparison report.

The controlled one-million-record comparison calls the actual old/new private
audit helper with the same inputs in separate JVMs. The measured additional heap
after explicit GC falls from 582,822,880 to 22,292,544 bytes, approximately 96%.
All serialized records have the same SHA-256:
`9678ec5b302e7e1dabe7aa9c4d2baafdf0fb5405fdb5759605b48cfc48aa2625`.
This measures audit-snapshot retention, not whole-application peak RSS. It uses a
3 GiB heap; reference widths may differ in the much larger production JVM.
This should not be substituted for the whole-application measurements below.

Job 10279 completes the full 9,790,188-candidate DPZ case through both
production arms on node002 with the same 64 GiB heap and 80 GiB allocation as
the baseline. Every TSV byte and every ledger line except exactly two
stage-duration fields match; complete candidate and input checks pass. It
explicitly places Java temporary files in job scratch. The final Slurm state is
`COMPLETED`, exit `0:0`. This success releases the dependent production waves.

| Arm | Baseline peak RSS, KiB | Compact peak RSS, KiB | Observed reduction | Baseline JVM time | Compact JVM time |
| --- | ---: | ---: | ---: | ---: | ---: |
| No priors | 60,425,520 | 51,117,972 | 15.4% | 7:28.16 | 6:44.92 |
| Default priors | 66,600,956 | 63,801,160 | 4.2% | 37:56.95 | 40:30.11 |

These are single paired observations, not repeated performance experiments.
Default-prior time increases, and total inference/verification time rises from
3,144.10 to 3,213.60 seconds. There is no general speedup claim. The
[dense receipt archive](../../benchmark/v3/results/compact-audit-dense-v1/)
retains both runs, comparison, resources and scheduler status; its verifier checks
the receipt chains, not the large files that remain remote.

The [261-file checkpoint](../../benchmark/v3/results/compact-audit-v1/) includes
both source snapshots, the final 771-file source manifest, failed-build evidence,
review, tests, memory measurements, baseline receipts and small old/new outputs.
Run `python3 benchmark/v3/results/compact-audit-v1/verify_archive.py` to check it.
Large profiles and production artifacts remain on the remote system at their
recorded paths; the portable archive does not re-read those remote bytes.

All twelve whole acceptance items remain open. This work does not establish
biological improvement, a full benchmark result or readiness for v3.0.0 release.
