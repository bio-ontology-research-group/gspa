# Gene-caller execution and coordinate-join repair

Checkpoint: 2026-10-28. **The scoped repair and native CLI validation pass.**
All twelve whole acceptance requirements remain open. This repair
addresses execution and identity integrity, not improved gene or function accuracy.

The [pre-repair audit](v3-gene-caller-status.md) independently reproduced
successful-zero misreporting, acceptance of incomplete output, order-dependent
wrong protein assignments and rejection of valid numeric contig names.

## Implementation

`AbstractProdigalCaller` now supplies the shared Pyrodigal/Prodigal contract:

- Redirect process output to a file before waiting, enforce a real deadline,
  terminate descendants on timeout, retain a bounded diagnostic tail and preserve
  interruption. Availability checks require the caller's recognizable version
  output and successful exit.
- Report fresh typed SUCCESS, SUCCESS_ZERO, FAILED, PARTIAL or UNAVAILABLE records,
  with native exit codes where available. An empty input is a failure to perform
  the requested scan; missing sequence on some contigs is explicit partial coverage.
- Require native GFF coverage records for all supplied contigs, including those
  with zero calls, and validate the complete CDS/protein correspondence before
  attaching anything. Existing CDS/proteins are preserved and refused as inputs
  for a new call; unrelated non-CDS features are retained.
- Preserve valid native output's public protein IDs and coordinate ordering.
  Reordered protein FASTA no longer controls the protein/CDS join.

`GeneModelReader` extracts the existing strict antiSMASH single-row CDS reader into
core. It resolves the union of native ID aliases using contig, start, end and
strand **before** checking uniqueness. This handles native metadata collisions on
numeric contigs, including identical coordinates on different contigs. The CLI
`loadCoordinates` entry point delegates to this reader. It is not a spliced
transcript assembler.

`AnnotationPipeline` preserves fresh caller statuses, reports configured caller
unavailability, catches caller exceptions and continues independent nucleotide
predictors. Legacy callers without a fresh typed status retain the conservative
zero-output failure fallback. No runtime fallback or confidence change was added.

## Native Prodigal acquisition

Job **10515** ran the official [Prodigal 2.6.3 Linux release](https://github.com/hyattpd/Prodigal/releases/tag/v2.6.3)
on the same inputs as native Pyrodigal observation 10507. The downloaded ELF is
13,599,280 bytes, SHA256
`80a3d7acb08ba9f63dde6b6ef4e95ec356128ac355f36960ea722f381150f246`.
The release API receipt, download receipt, version output and commands are retained.
All seven cases complete; the job exits zero after 110 seconds.

| Input/mode | Native exit | CDS/proteins |
|---|---:|---:|
| Full assembly, single | 0 | 4,359 |
| Full assembly, meta | 0 | 4,286 |
| Nine bases, meta | 0 | 0 |
| Nine bases, single | 10 | 0 |
| 10,000 Ns, meta | 0 | 1 |
| 10,000 Ns, single | 10 | 0 |
| Empty input, meta | 18 | 0 |

The version probe `-v` exits zero and prints to stderr. Whole-assembly single-line
FASTA and native contig coverage comments work. Each caller must be compared with
its own native reference: Pyrodigal meta mode returned 4,291 proteins, not 4,286.
Neither count is a biological gene-accuracy score. Successful zero output does
not establish biological gene absence.

## Validation history and current jobs

| Job | Result and scope |
|---|---|
| 10513 | Failed predictor test compilation due to a Spock multiple-variable declaration. Old CLI numeric-contig regression separately failed as intended. A missing test-result directory also caused the wrapper's exit 23. Retained; not a passing gate. |
| 10514 | Original production files: 41 cases, 22 expected failures. Repaired focused predictor run: 31/32 pass; version-name matching fails and is corrected. CLI XML in its final-results directory still belongs to the earlier failing-before phase because Gradle stopped before rerunning CLI tests. |
| 10517 | Source v3: 42 original-production cases with 23 failures; **53 focused cases pass**; full **1,548 JVM cases, six optional skips, no failures/errors**. Exit zero, 404 seconds. |
| 10519 | Cancelled only while verified pending, before any execution, to replace its planned v3 CLI validation with final v4 validation. |
| 10520 | Final v4 source: **57 focused cases pass; 1,552 full JVM cases, six optional skips, no failures/errors**. Exit zero. All 843 frozen source files still match the worktree. |
| 10521 | All **23 actual native/CLI invocations pass**, including before/after annotation and antiSMASH imports; exit zero after 432 seconds. |

The interim 10517 JAR has SHA256
`912151490e0ef1b697ac79ec24df88a31fb14c30f0adc6d88daf8dfb0c490300`.
It is not substituted into the frozen Price collection or source/context job 10406.

The final 10520 JAR has SHA256
`1ae842e5905eecfff4ecbcb07bfa21d320ff9de0797273e4aeafd2abaccbdb5d`,
at `/data/hohndor/gspa-v3/gene-caller-status-v1/jobs/10520/gspa-3.0.0.jar`.
Only the new native route validation uses it; existing frozen benchmarks stay unchanged.

Two read-only Sonnet 5.5 repair reviews and adjudication are retained. Applied
findings include identical-coordinate contig discrimination, second-contig
atomicity, wrong-end rejection, exact version arguments, a native Prodigal zero
fixture, interrupt preservation, more tolerant timeout-test timing and complete
CLI timeout records. Speculative changes to ID ordering, runtime fallback and
unverified translation rules were not adopted.

## Final experiment and independent verification

`benchmark/v3/check_gene_caller_contract.py` invokes the **actual annotation CLI**
before and after the repair. It runs both native callers in single/meta modes on
the full assembly, compares all sequences and coordinates against each caller's
independent native GFF/FASTA, and requires byte-identical valid annotation exports.
It also tests native short inputs, numeric contigs, explicitly reordered/missing
captured output, unavailable execution, and byte-identical unscored antiSMASH
imports through the extracted reader. Source, JAR, native files, wrappers and
runtime identities are bound and rechecked. All checks pass. Native failures keep
their actual exit codes; successful zero calls have SUCCESS_ZERO; unavailable
execution is UNAVAILABLE. Reordered output maps correctly, and incomplete output
is FAILED with no attached CDS or proteins. Both antiSMASH imports produce the
same **65 unscored memberships** for the 4,359-protein Pyrodigal gene set. This is
a different input gene set from the earlier 70 memberships on 4,915 supplied CDS.

The [portable repair archive](../../benchmark/v3/results/gene-caller-repair-v1/README.md)
retains 826 source/output bindings, failed attempts, complete passing gates and
the actual native outputs. Its independent verifier passes the full checksum
inventory, all four scoped source snapshots, native sequence/coordinate and export
parity, typed status controls, atomicity, numeric joins and import parity. It
checks 156 retained runtime input bindings and explicitly enumerates seven omitted
runtime binaries/wheels. The archive contains 830 files plus its checksum inventory.

```bash
python3 benchmark/v3/results/gene-caller-repair-v1/verify_archive.py
```

Timeout/process tests use controlled processes; neither
those nor native output parity certify every possible wrapper process tree,
caller version, genome class or biological prediction.

Remote root: `/data/hohndor/gspa-v3/gene-caller-status-v1`.
Local preparation: `/home/leechuck/tmp/gspa-gene-caller-status-v1`.
Frozen final source: `repair-source-v4.tar.gz` / `repair-source-v4.json` (843 files).
Launchers: `repair-v4.sbatch`, `cli-v2.sbatch`.
