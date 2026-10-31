# dbCAN annotation and Nextflow integration

The dbCAN adapter was not registered by the annotation pipeline, and Nextflow's
native overview output did not enter the integration claims stream. This follow-up
to the [parser repair](v3-dbcan-contract.md) connects those routes. It does not
establish improved enzyme prediction or complete any whole audit item.

## Changes

Explicitly enabling `dbcan` now registers the adapter. Configuration supplies the
database directory, executable, interface version and CPU limit; CLI `--threads`
is forwarded. Explicit disable takes precedence. Both interfaces require an
explicit database directory rather than relying on a default in a temporary
working directory. Version 4 uses the positional interface; version 5 uses
`run_dbcan CAZyme_annotation`. The version-5 probe uses `run_dbcan version`;
the legacy tool version remains unknown rather than being invented.

The command is checked against the pinned
[dbCAN 5.2.8 entry point](https://github.com/bcb-unl/run_dbcan/blob/3c2b5403ee9df87932fb7881b3e6b51101b5e7f0/pyproject.toml)
and [CLI source](https://github.com/bcb-unl/run_dbcan/blob/3c2b5403ee9df87932fb7881b3e6b51101b5e7f0/dbcan/main.py).
The executable in that release is `run_dbcan`. Selecting an interface does not
certify which release or database is installed.

`gspa dbcan-claims` imports an existing native overview through the same parser.
It rejects malformed tables, duplicated FASTA identifiers and unknown output
identifiers, including zero-hit rows. Its receipt records input/output SHA-256
hashes and counts, and explicitly says it imports existing output rather than
attesting successful tool execution. Claims and receipt are staged before
publication, with cleanup on ordinary publication errors.

Nextflow converts native dbCAN results through this command and merges them with
builtin and sidecar claims. A missing overview now fails the native process even
if the executable exited zero. The builtin intermediate claims file has its own
published filename, preventing it from overwriting the merged claims file.
The complete replay also exposed an empty-component failure: the builtin parser
aborted a dbCAN-only run because it produced no builtin claims. Its new explicit
`--allow-empty` option lets Nextflow compose empty builtin streams with other
predictors. Standalone benchmark parsing still fails on empty output by default.

Family observations remain unscored. `proteinToolCount` describes methods with
any hit for the protein; per-family supporting methods are separate. Explicit
JSON metadata copies retain original fields, subfamily identifiers and domain
coordinates despite the shared writer's stringification of nested objects.
No EC assignment or calibrated confidence is inferred from the overview count.

## Configuration example

```yaml
predictors:
  enable: [dbcan]
  dbcan:
    database: /data/dbcan
    executablePath: /tools/run_dbcan
    interfaceVersion: '5'
    threads: 4
```

With Nextflow, use `--run_dbcan true --dbcan_db /data/dbcan` and enable integration
with its usual reference inputs. Imported families are retained in raw claims;
they are not automatically promoted to scored integrated predictions.

## Verification and limits

The [validation archive](../../benchmark/v3/results/dbcan-route-v1/README.md)
retains the source snapshots, tests, Sonnet 5.5 review and adjudication, failed
attempts and route logs. The route runner uses three synthetic proteins and
controlled native executables, with real GSPA CLI and full Nextflow execution.
Its version-5 fixture adapts the validated version-4 formatter output's column
name; it is not output from a complete version-5 search.

The final focused run passes all 27 cases. In the before comparison, eight of
12 route cases fail, demonstrating missing support rather than a biological
accuracy result. The complete final predictor/CLI gates pass 502 cases with
four optional skips. The unchanged core separately passes 909 cases with two
optional skips; all 308 packaged core source/resource files are hash-identical
between those runs. This is 1,411 cases across the three modules, not a single
new aggregate test execution.

The first Nextflow attempt failed before execution because an old launcher
rejected Java 21. A matching
[pinned Nextflow 25.10.4 launcher](https://github.com/nextflow-io/nextflow/blob/v25.10.4/nextflow)
is included in the isolated replay package. The corrected launcher exposed the
empty builtin component failure described above. After repair, job 10189 exits
zero: all ten replay invocations have their expected status. Both CLI interfaces
preserve all four observations byte-for-byte through annotation/import; the full
Nextflow run preserves the same version-5 claims. The fault run fails specifically
at DBCAN when its executable exits zero without an overview. Independent checks
of retained receipts, JSON metadata, input counts and process traces agree.
No unscored observation becomes a scored prediction. The exact tested JVM JAR
is reused for the final Python/Nextflow replay.

Actual dbCAN database searches, broader native output coverage, displayed
unscored observations, calibrated family evidence and independent biological
accuracy still require validation. All twelve whole-item acceptance gates remain open.
