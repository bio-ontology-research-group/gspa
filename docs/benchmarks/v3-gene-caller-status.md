# Native gene-caller outcomes and production mapping audit

Follow-up: [repair implementation and current validation](v3-gene-caller-repair.md).
This document and its portable archive preserve the original pre-repair evidence.

Checkpoint: 2026-10-28. This is **pre-repair evidence**. Four remote jobs exited
zero and independently retained checks reproduce three production defects. No
production code was changed; no gene-prediction or functional-accuracy improvement
is claimed. All twelve whole acceptance requirements remain open.

## Native observations

Job **10507** ran pinned Pyrodigal 3.7.1 with archspec 0.2.6 on node004, using
CPython 3.12.3, two CPUs and 4 GiB. The full 145-contig GCA_003598085.1 assembly is
the same nucleotide input as the earlier [operon route audit](v3-pyrodigal-operon-route.md).
Additional short, all-N and empty inputs are explicit contract probes, not
biological truth. The seven cases completed in a 36-second job.

| Input | Mode | Native exit | CDS/proteins | Output interpretation |
|---|---|---:|---:|---|
| Full assembly | single | 0 | 4,359 | Complete observed output |
| Full assembly | meta | 0 | 4,291 | Complete observed output on the same assembly |
| Nine bases, `ATGAAATAA` | meta | 0 | 0 | Three-comment-line GFF and existing empty FAA: genuine successful zero-call run |
| Same nine bases | single | 1 | 0 | Training error: fewer than 20,000 bases; empty GFF and FAA |
| 10,000 Ns | meta | 0 | 1 | Partial CDS and an all-X protein; **not** a zero-call fixture |
| 10,000 Ns | single | 1 | 0 | Training-length error; empty outputs |
| Empty FASTA | meta | 0 | 0 | Empty GFF and FAA; no input was assessed |

A successful zero-call run does not establish that an organism is gene-free.
Likewise a tool's success on an empty input is not evidence of a completed genome
scan. The all-N output is retained without treating its score as functional
evidence. No native Prodigal executable was found on this worker's PATH.

## Actual production behavior

Job **10509** invoked the unchanged production `AnnotationPipeline` and
`PyrodigalCaller` through a Groovy harness using the validated 10481 JAR:

```text
/data/hohndor/gspa-v3/antismash-nextflow-v1/jobs/10481/gspa-3.0.0.jar
SHA256 069263f7b876222815928a637d02106b9c3b31a756edd1cf692432acc2e52303
```

This is execution of the shared production classes, not a new CLI route or a
mock inference engine. Six calls completed in a 69-second job. Four used actual
native Pyrodigal: full assembly single/meta and nine-base input single/meta.
Two used explicitly corrupted captured native output to test the adapter contract.
The pipeline had integration disabled and no protein predictors registered, so
these checks make no assertion about downstream function prediction.

1. **Successful zero and native failure collapse into the same status.** Both
   short-input modes produce `FAILED`, identical “failure vs genuinely gene-free
   input ... unverified” detail, and no exit code. Native meta mode actually exited
   zero; single mode exited one. The caller interface supplies no typed status.
2. **Protein/CDS identity depends on record order.** Reversing the complete native
   FAA, with its headers and sequences unchanged, yields 4,286 wrong sequence/
   coordinate assignments among 4,359 proteins. The pipeline still reports SUCCESS.
   The existing strict antiSMASH coordinate reader on the same reordered input
   has zero mismatches. This corruption is injected; actual native full-assembly
   output in both modes matches every independently parsed sequence/coordinate.
3. **Incomplete output attaches orphan features and shifts assignments.** Removing
   only the first protein record leaves 4,359 CDS features and 4,358 proteins,
   including 190 wrong sequence/coordinate assignments, yet reports SUCCESS.

The independent verifier derives coordinates directly from native FASTA headers
and compares full protein sequences. It does not infer correctness from counts
alone or reuse the caller's order-based join.

Source inspection also confirms that both JVM callers read `proc.inputStream.text`
before `waitFor(30, MINUTES)`, then ignore the wait result and call `exitValue()`.
This does not enforce the intended process timeout. It remains a code finding;
these jobs did not run a 30-minute hang experiment.

## Numeric-contig ambiguity in the proposed shared reader

Sonnet 5.5 review identified a separate potential defect in
`AntiSmashClaimsCommand.loadCoordinates`. Jobs **10510** and **10511** reproduce
it using native outputs, so it must be repaired before extracting that reader.

Job 10510 preserves the first two assembly nucleotide sequences exactly and only
renames their FASTA headers to `2` and `1`. Actual Pyrodigal single/meta calls with
one thread produce 294/291 proteins. For the first protein, the native FASTA
identifier is `2_1`, header metadata says `ID=1_1`, and GFF ID is correctly `2_1`.
The reader unions aliases from both identities before checking coordinates, so
an unrelated CDS on contig `1` introduces ambiguity. Job 10511 invokes the
unchanged reader and reproduces rejection in both modes:

```text
missing or ambiguous CDS identity for 2_1
```

Independent reconstruction finds 208/202 such alias collisions. Filtering the
candidate set using the native header's contig, coordinates and strand yields
exactly one match for every protein. This is a strict-reader ambiguity, **not**
the previously withdrawn claim that pinned Pyrodigal's GFF IDs need normalization.
The native GFF/FASTA identifiers still match.

## Next repair and experiments

Implement a shared gene-caller process/status path that actually bounds waiting,
cleans up descendants, preserves diagnostic tails and native exit codes, and
distinguishes SUCCESS_ZERO, FAILED, PARTIAL and unavailable execution. Preserve
independent nucleotide predictors after caller failure and the conservative
legacy-caller fallback. Validate all output before attaching any model to the
input genome. Preserve established IDs and ordering on valid full output.

Extract a shared single-row CDS reader after fixing the observed numeric alias
case. Resolve candidates using header identity **and** contig/coordinate/strand
agreement, with exactly one consistent match; reject missing/duplicate/ambiguous
matches. Do not add an unrelated runtime fallback, change ID ordering, or impose
unverified translation-length rules. Genuine zero must retain the successful
native status; missing output or missing CDS/protein partners must fail without
partial mutation. Missing input sequences need explicit coverage semantics.

Required proof: focused failing-before/passing-after regressions; actual native
short/normal/numeric cases; missing/reordered/duplicate output injections; bounded
timeout and descendant cleanup; stale-status and independent nucleotide-tool
continuation; full remote JVM gate; affected antiSMASH import/CLI route checks;
and a complete before/after sequence/coordinate/ID comparison on the 145-contig
assembly. The four jobs here supply the pre-repair fixtures and observations.
They do not replace the eventual biological evaluation.

## Evidence and reproduction

- Remote root: `/data/hohndor/gspa-v3/gene-caller-status-v1`.
- Local preparation: `/home/leechuck/tmp/gspa-gene-caller-status-v1`.
- [Portable archive](../../benchmark/v3/results/gene-caller-status-before-v1/README.md):
  native inputs/outputs, commands, exit receipts, six pipeline outputs, numeric
  reader failures, ten source files matching the frozen build, Sonnet review and
  adjudication. Large files are gzip-compressed with decompressed hash bindings.
- The unchanged production JAR, runtime and wheels are hash-bound but not duplicated.
  An inherited unused base-runtime JAR copied by the 10509 finish hook is explicitly
  excluded; it was not the JAR executed by the harness.
- No new JVM suite was compiled or run. The earlier 1,521-case gate belongs to
  10481, not to these diagnostic jobs.

```bash
python3 benchmark/v3/results/gene-caller-status-before-v1/verify_archive.py
```

The verifier checks the archive inventory, all 121 retained/decompressed bindings,
all four exit receipts, native file hashes/counts, every normal sequence/coordinate
pair, all corrupted-output mismatches, the successful strict reordered join and
the independently reconstructed numeric alias collisions. All checks pass.
