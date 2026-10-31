# dbCAN route validation

See the [change report](../../../../docs/benchmarks/v3-dbcan-route.md).
This archive verifies software integration with controlled native executables.
It is not a complete dbCAN database search or a biological benchmark.

The final focused gate passes 27 cases. Eight of 12 route cases fail against the
before sources. Job 10187 passes all 421 predictor and 81 CLI cases, with four
optional skips. Job 10185 separately passes 909 core cases, with two optional
skips; `unchanged-core-proof.json` confirms all 308 packaged core source/resource
files are unchanged. These three module gates cover 1,411 cases in total.

Job 10189 exits zero after ten command invocations: both native-interface
annotation/import/integration routes, standalone builtin-parser default/allowed
empty controls, and complete Nextflow success/missing-output runs. Four unscored
family observations survive byte-for-byte through each successful route. The
missing-output trace fails specifically at DBCAN and does not complete downstream
claim conversion or integration. No observation is promoted to a confidence score.

Run this independent retained-output check from any directory:

```sh
python3 verify_archive.py
```

It checks all retained XML counts, input/output receipt hashes, parseable raw
metadata, integration input counts, merged claims and Nextflow process statuses.
`SHA256SUMS` binds the entire archive. The 45 MB tested JAR is retained remotely
at `/data/hohndor/gspa-v3/dbcan-route-v1/jobs/10187/gspa-3.0.0.jar`; its hash is
recorded in the accepted replay report and `job-10189/jar.sha256`. The replay
uses that exact JAR, with only Python/Nextflow changes afterward.

`source-after/` retains the final scoped sources; `source-before/` retains the
prior sources. The final complete source package's 737 input hashes are in
`source-package-manifest.json`. Runtime/package hashes and exact Slurm scripts
are retained. Source packages v2–v5 and the runtime packages remain in the remote
paths named by those scripts. Application execution stayed on node002 under
Slurm. The source archives deliberately exclude concurrent unrelated edits.

The upstream 5.2.8 source receipts establish executable/argument/output contracts.
The version-5 controlled overview adapts the version-4 formatter fixture; it is
not a genuine full version-5 search result. `sonnet-review.txt` and
`review-adjudication.md` retain the read-only Sonnet 5.5 review and dispositions.

Failed attempts 10183–10188 are retained and explained in `validation.json`.
Some failure archives contain stale earlier-phase XML copied by the exit trap;
only the explicitly named accepted gates above support the final passing count.
In particular, 10187 is the retained reproduction of the real Nextflow
empty-builtin-stream failure. Standalone parser failure remains the default;
Nextflow explicitly opts into empty component streams.

Real tool/database identity, broader native format coverage, biological accuracy,
calibrated family evidence and useful web presentation of unscored observations
remain open. No whole audit item is accepted by this archive.
