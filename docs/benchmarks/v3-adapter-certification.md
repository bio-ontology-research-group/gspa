# Row 3 adapter certification: genuine tool outputs, dispatch paths, injected failures

Bar (`tasks/v3-acceptance.md`, row 3): "Versioned genuine tool outputs and complete
parser/dispatch paths, including injected failures and coverage loss."

Frozen plan: `tasks/v3-row3-adapter-certification-2026-11-01.md` (sha256 `4e51cf11…1148`,
recorded before any check ran in `benchmark/v3/results/adapter-certification-v1/PLAN.sha256`;
all later changes are under its "Dated deviations" heading). This is adapter evidence
(format, dispatch, attribution, typed outcomes), not a biological benchmark.

## What was run

One generic harness, `benchmark/v3/check_adapter_route.py`, driven by a per-adapter
`spec.json`, with oracles in `benchmark/v3/adapter_oracles.py` that read the native tool
output without GSPA code. For every adapter it runs the real CLI JAR
(`gspa-cli/build/libs/gspa-3.0.0.jar`, rebuilt after every JVM change, `./gradlew build` green:
538 gspa-predictors tests, 0 failures) through:

1. `annotate` with the genuine tool on PATH → raw claims (`--claims-out`) must equal the oracle
   (protein, type, function, score, hit); `annotate --integration-config` and direct `integrate`
   on the raw claims must give byte-identical integration tables.
2. Injected cases through a stub executable that replays variants of the genuine output:
   `header_only` (documented empty result), `malformed` (non-numeric statistic),
   `truncated_row`, `dropped_protein`, `unknown_protein` (a row naming a protein never
   submitted), and `unavailable` (executable absent). Every expected outcome is an exact
   `PredictorOutcome` read from `<genome>_predictor_status.json`.
3. Genome-level tools: an independent Python overlap of region coordinates with GFF CDS must
   equal the proteins carrying the region claim; a region overlapping no CDS must still be in
   the GFF3 output; a region outside the contig or on an unknown contig must make the run
   `PARTIAL` naming it.

Every input, database file, model file and the JAR are hashed before and re-hashed after
execution (`verification.json`, schema `gspa.adapter-route-replay/1`). Each adapter folder
under `benchmark/v3/results/adapter-certification-v1/<adapter>/` holds `spec.json`,
`observed/` (genuine output, tool version, command, wall time), `before/` (engine as found)
and `after/` (engine after repairs), each with `runs.json`, per-run logs and `checks.json`.

Inputs: Swiss-Prot P18843 (NadE), P0ADB7 (EcnB) and the synthetic 80-alanine control;
a 285-entry Swiss-Prot reference (UniProt REST, hashed) with GO/EC in the FASTA headers;
VFDB `VFDB_setA_pro.fas` (2026-11-01, sha256 `5b86445f…`); AntiFam 8.0; dbCAN
`db_v5-2_9-13-2025`; ProstT5 `prostt5-f16.gguf` (sha256 `b072042c…`); and for nucleotide
tools *E. coli* BW25113, RefSeq GCF_000750555.1 (FASTA, GFF3, GenBank). All work ran under
`/mnt/data2/leechuck/gspa-cert` on ws (pinned downloads hashed in `downloads/SHA256SUMS`);
geNomad, CheckV and the gapseq nucleotide route on Shaheen compute nodes (`benchmark/v3/shaheen/
adapter-cert-{genomic,genomic-harness,gapseq,gapseq-harness}.sbatch`, specs from
`make_genomic_specs.py`), and the gspa-nf containers through Slurm + Singularity
(`nf-shaheen.config`, `nf-shaheen-run.sh`; head process on a login node, every task on a
compute node). Shaheen outputs are copied to the result folders; every file over 5 MB is
listed in `LARGE_FILES_REMOVED.tsv` and kept under `/mnt/data2/leechuck/gspa-cert/`.

## Results per adapter

| Adapter | Tool, version | Genuine run | Checks (after) | Before (engine as found) | Verdict |
|---|---|---|---|---|---|
| diamond | DIAMOND 2.1.11 | 15 hits → 85 GO/EC claims = oracle | 3125, 0 failed | 58 failed: truncated rows skipped silently (SUCCESS, 78 claims); a row naming an unsubmitted protein accepted (SUCCESS, 85 claims) | PASS after repair |
| mmseqs2 | MMseqs2 17-b804f | 19 hits → 94 claims = oracle | 139, 0 failed | UNAVAILABLE with the tool on PATH (no `--version` flag); script ignored `executablePath`; truncated/unknown rows accepted | PASS after repair |
| foldseek | FoldSeek 10-941cd33, ProstT5 (CPU) | 276 hits → 1937 claims = oracle | 1974, 0 failed | same availability defect as MMseqs2 (not separately re-run; repaired in the shared probe) | PASS after repair |
| vfdb | DIAMOND 2.1.11 vs VFDB setA | 5 hits (FimE, EcpA) = oracle | 3045, 0 failed | ABSENT from the manifest: the adapter was never registered by the pipeline | PASS after repair (now `predictors.enable: [vfdb]`) |
| antifam | HMMER 3.4, AntiFam 8.0 | seed member Q8EZL8/38-87 → ANF00001; the three real proteins clean | 57, 0 failed | no engine defect; the earlier "synthetic model only" fixture replaced by genuine models | PASS |
| metapredict | metapredict 3.1.01, torch 2.5.1 | 2 disorder regions = oracle | 79, 0 failed | UNAVAILABLE (probe imported torch, >10 s); with the tool reachable the driver crashed on numpy arrays (`… or []`); single-protein failures invisible | PASS after repair (`PARTIAL` names unscored proteins) |
| dbcan | run_dbcan 5.2.8 | GH13_19/GH13 (amylase), GT35 (glycogen phosphorylase) = oracle | 72, 0 failed | no engine defect | PASS |
| crispr | MinCED 0.4.2 on BW25113 | 2 arrays (2,871,060–2,871,822; 2,897,373–2,897,766), intergenic, retained in GFF3 | 43, 0 failed (CDS-translation route); 23, 0 failed with the RefSeq FAA, where the run is `PARTIAL` naming the 6 multi-copy WP_ accessions that cannot be placed | ABSENT (never registered); contig-keyed results would have been dropped at protein assignment; after the first repair the placeholder contig was reported as "contig without sequence" | PASS after repair |
| phispy | PhiSpy 4.2.21 (biopython 1.81) via genomic sidecar | 5 prophages (DLP12, e14, Rac, Qin, CP4-57); same calls from NCBI GenBank and from GSPA's GenBank, but the e14 end differs between runs (PhiSpy is not deterministic), so the oracle is captured from the same run | 61, 0 failed; attribution 130 of 130 oracle proteins of the same run after the CdsTranslator repair (128 of 130 before it) | with biopython 1.88 PhiSpy crashes; the adapter reported `FAILED` with the missing-output diagnostic (sidecar exits 0) | PASS |
| genomad | geNomad 1.13.0 (sif), DB v1.9 | genuine sidecar runs on Shaheen (nid00155, 276 s and nid00349, 264 s; identical regions, so geNomad is deterministic between runs): 5 proviruses (559,304–579,886; 1,192,323–1,208,563; 1,406,257–1,428,481; 1,621,774–1,650,406; 2,456,491–2,469,989; scores 0.887–0.963), the same five prophage loci PhiSpy reports; 139 CDS attributed = oracle | 150, 0 failed (job 15745896) | no engine defect; two harness-caused failed runs kept: job 15745762 (CDS-translated ids not mapped back to GFF ids, earlier JAR) and job 15745876 (harness re-uploaded during the run, `unchanged:` check) | PASS |
| checkv | CheckV 1.0.3 (sif), DB v1.5 | genuine sidecar run on Shaheen (nid00349, 68 s): one `viral_contig` region = the whole 4,631,469 bp chromosome (High-quality, completeness 100.0, 92 viral genes), so all 4,217 CDS carry the claim = oracle; see the tool-domain observation below | 150, 0 failed (job 15745876; the earlier JAR also 150, 0 in job 15745762) | no engine defect (first direct attempt: `checkv_reps.dmnd` missing from the DB tarball, built in-job) | PASS |
| signalp | SignalP 6 | not runnable: DTU licence download required | – | – | NOT CERTIFIABLE |
| deeptmhmm | DeepTMHMM (pybiolib) | not runnable offline: remote DTU execution, account required | – | – | NOT CERTIFIABLE |
| gapseq (nucleotide route) | gapseq 2.1.0 (Shaheen conda env; sequence DB 1.5; blast 2.17.0), `find -p all -t Bacteria -M nucl` | genuine run on BW25113 (449 s, nid00635): 1,929 pathways, 155,682 reference sequences, 39,913 Reactions.tbl data rows (good_blast 8,201; bad_blast 26,717; no_blast 3,626; no_seq_data 1,024; spontaneous 345). Through the JAR (`-K 32`, about 5 min per run) the documented outcome: `FAILED` "nucleotide pathway hits require coordinate-to-protein attribution; use target=proteome for protein claims", zero claims, the same in both genuine runs. gapseq 2.1.0's `-M nucl` calls its own ORFs (4,314) and names hits by them, which is exactly what the contract refuses to attribute | 608, 0 failed (job 15745997): header-only table `SUCCESS_ZERO`, truncated row `FAILED` ("malformed reaction row 5"), absent executable `UNAVAILABLE`; parity on the empty claim set byte-identical | first harness attempt crashed on the `# gapseq version` comment lines of Reactions.tbl (harness fixed); the second (job 15745883) passed `--threads 1`, which overrides `pathway.threads`, so the JAR ran gapseq with `-K 1` and the harness's 3600 s limit expired before any check (harness fixed: spec keys `cli_threads`, `timeout_seconds`); no engine defect | PASS (the typed failure is the adapter's contract) |
| gspa-nf `main.nf` containers | Nextflow 25.10.4 + Singularity 4.3.2 on Shaheen (Slurm, k10018) | fresh run on BW25113, `nextflow rc=0`, 9 min 25 s: all six processes `COMPLETED` exit 0 in their images (quay.io biocontainers pyrodigal 3.7.1, diamond 2.1.9, hmmer 3.4, barrnap 0.8, minced 0.3.0; python:3.12 for the merge); outputs: 4,314 proteins, 19 DIAMOND rows, 14,072 Pfam domain rows, 22 rRNA rows, 37 MinCED rows, merged table 14,188 rows over 4,138 proteins, annotated GFF3 with 4,314 CDS (4,116 with Pfam) | `trace.tsv` with status, exit and container per task; independent recount: every merged row names a called protein, Pfam rows = domtbl rows (14,072), GO/EC rows (94) = regex over the DIAMOND `stitle` column, rRNA rows = barrnap rows (22) | first attempts: `$HOME` quota broke SIF builds (caches moved to /scratch); attempt 3: DIAMOND database not yet built (setup); attempt 4: five processes completed but MERGE_ANNOTATIONS exit 1 before its script ran — `python:3.12-slim` has no `ps`, which Nextflow's task metrics need under `-with-trace`/`-with-report` (pipeline defect, repaired), and Nextflow had kept the stale attempt-3 trace (`trace.overwrite` off) | PASS after repair |
| antismash, interproscan, hmmer, eggnog-mapper, amrfinder, gapsmith, prodigal, pyrodigal | previously certified | carried forward (antiSMASH remains "captured, not fresh") | | | carried |
| operons; neural/region/site sidecars | JVM-only; model weights | out of the bar's "tool outputs"; CLEAN certified separately | | | out of scope |

Every typed outcome in the after runs is exact: `SUCCESS` / `SUCCESS_ZERO` / `PARTIAL` /
`FAILED` / `UNAVAILABLE` as the spec demands, with zero claims in every failure case and the
offending protein or contig named in the detail.

## Finding repaired in a second pass (gspa-core)

- `gspa.io.CdsTranslator` grouped every GFF CDS sharing a `protein_id` as segments of one
  protein (the split-CDS rule), so identical multi-copy proteins (IS transposases:
  WP_000019403.1 ×9, WP_085947771.1 ×5 in BW25113) collapsed to their first locus and copies
  inside a prophage were never attributed (PhiSpy attribution 128 of 130 before the repair).
  CDS rows are now grouped per locus (protein key qualified by the gene `Parent`, falling back
  to the feature `ID`), so a programmed frameshift still translates as one protein while each
  copy of a multi-copy gene is its own protein with its own coordinates (BW25113: 4,192 →
  4,217 translated proteins; the 25 extra are the additional loci of multi-copy genes). Regression: `CdsTranslatorSpec`
  ("identical proteins at several loci stay separate", "segments of one split CDS still
  translate as one protein"). After the repair PhiSpy attribution is 130 of 130 with no excused
  gap, and CRISPR (43 checks) and PhiSpy (61 checks) re-verified. The `--proteins` route is
  unchanged: a RefSeq FAA still cannot place its 6 shared accessions and stays `PARTIAL`
  naming them.
- PhiSpy is not deterministic between two runs on the same GenBank (the e14 prophage ended at
  1,212,452 in one run and 1,212,091 in the next, one protein's membership changed). The
  harness therefore captures the sidecar TSV of the run that produced the claims
  (`capture_native`), checks attribution against that run, integrates the claims of the
  integrated run itself for route parity (`nondeterministic: true`), and records the
  between-run difference as an observation.

## Tool-domain observation (adapter correct, tool misapplied)

CheckV 1.0.3 run through the genomic sidecar on the complete BW25113 chromosome reports the
whole 4,631,469 bp contig as a `High-quality`, completeness 100.0 viral contig (92 "viral
genes"), so every CDS receives a `VIRAL_CONTIG` claim with score 1.0. The adapter route is
verified (150 checks, the attribution oracle agrees because the region is the whole contig),
but CheckV is a quality estimator for viral contigs, not a classifier; on a bacterial
chromosome its output is wrong and GSPA passes it through unchanged. Recorded here so nobody
reads the CheckV claims on complete bacterial genomes as evidence; gating CheckV on prior
viral classification (geNomad) is a workflow decision outside this round.

## Defects found and repaired (defaults unchanged)

- `AbstractToolPredictor.isAvailable` required `--version`; MMseqs2 and FoldSeek only have a
  `version` subcommand → installed tools were reported `UNAVAILABLE`. Adapters now declare their
  `versionArguments()`; `--help` remains the fallback.
- `MMseqs2Predictor` hard-coded `mmseqs` in its shell script; it now uses the configured
  executable and refuses an unset database.
- DIAMOND, MMseqs2, FoldSeek and VFDB parsers skipped short rows (`if (fields.length < n) return`)
  and did not catch `NumberFormatException`; both are now `ToolOutputContractException`
  (zero-byte output remains the documented no-hit result).
- Output naming a protein absent from the submitted input was silently dropped at assignment.
  `AbstractToolPredictor.predictBatch` now fails the run (`requireKnownProteins`).
- `VfdbPredictor` and `CrisprPredictor` were never instantiated by `AnnotationPipeline`; both are
  now registered through `predictors.enable: [vfdb|crispr]` with `VfdbConfig` / `CrisprConfig`.
- `CrisprPredictor` returned contig-keyed results that `applyAnnotations` would have dropped; it
  is now a `GenomePredictor` that keeps arrays in `genome.metadata["genomic_regions:crispr"]`
  (exported to GFF3) and projects them onto overlapping CDS through the shared
  `projectRegions`, unscored (`score_available=false`), with `PARTIAL` for invalid rows.
- `MetapredictPredictor`: availability probe uses `importlib.util.find_spec` (no torch import);
  the driver no longer tests numpy arrays for truth and writes `coverage.tsv` (one `ok|failed`
  row per protein); a protein the driver failed on or never reached makes the run `PARTIAL`
  naming it; malformed rows are contract failures instead of being skipped.
- Genome-level predictors described the pipeline's placeholder contig for proteins without a
  CDS placement as "a contig without nucleotide sequence". They now report the unplaceable
  proteins by id ("no CDS placement") — still `PARTIAL`, now naming the right thing.
- gspa-nf `MERGE_ANNOTATIONS` ran in `python:3.12-slim`, which has no `ps`; Nextflow's task-metric
  collection (`-with-trace`, `-with-report`) exits the task with status 1 before the script runs.
  `slurm_singularity.config` had already overridden this for unimatrix; the shipped default
  (`modules/quality.nf`, `nextflow.config`) now uses `python:3.12` (same interpreter). `MDF` and
  `FLATTEN_PROPHAGE_CDS` still use the slim image and were not run in this round.

Regression tests: `gspa-predictors/src/test/groovy/gspa/predictor/AdapterCertificationSpec.groovy`
(each case failed against the engine before the repair); existing specs updated where the
old behaviour was asserted (`AdapterRunContractSpec`, `MetapredictPredictorSpec`,
`AllPredictorParsingSpec`, `AmrFinderObservedSpec`).

## Stored observed outputs

`gspa-predictors/src/test/resources/observed-output/{diamond,mmseqs2,foldseek,vfdb,crispr,antifam,metapredict,dbcan/v5,phispy,genomad,checkv,gapseq}/`
each with `provenance.json` (tool version, command, wall time, input and output hashes). The
gapseq fixture keeps the full Pathways.tbl and a deterministic subset of Reactions.tbl (comment
lines, header, first 500 of 39,913 data rows) with the full file's sha256 and row counts.

## Honest limits

- The two Swiss-Prot proteins and one synthetic control are a format/dispatch fixture; nothing
  here measures annotation quality.
- The before-state runs for antifam used a stub whose banner the HMMER probe rejected; its six
  "before" failures are harness-caused, not engine defects, and are reported as such.
- SignalP 6 and DeepTMHMM cannot be obtained offline; their parsers are checked only against the
  documented formats (existing Spock fixtures).
- Shaheen compute nodes resolve `python3` to a 3.6 pyenv shim; the sidecar needs ≥3.7, so the
  job pins `/project/k10018/sw/miniforge3/bin/python3`.
- The gapseq nucleotide route is certified as a typed failure, not as a producer of claims: with
  gapseq 2.1.0 `-M nucl` is a translate-then-blastp with gapseq's own ORF calls, so the route
  cannot yield protein claims; `target=proteome` (gapsmith or gapseq prot) is the working route.
- gspa-nf `MERGE_ANNOTATIONS` reports "With product name: 0" here because the certification
  reference carries GO/EC in plain headers rather than `sp|…|… OS=` titles; it also does not read
  the MinCED GFF (CRISPR arrays stay in `minced/` only). Both are observations, not repairs.
- The Nextflow certification is a fresh run (not `-resume`) so that the trace records every
  process as `COMPLETED`; the previous attempts' logs are kept under `gspa-nf/logs/`.

## Proposed row-3 text for `tasks/v3-acceptance.md` (not applied; the supervisor edits that file)

> Row 3 — Versioned genuine tool outputs and complete parser/dispatch paths, including injected
> failures and coverage loss: **Met with stated limits.** Eleven tool adapters (DIAMOND, MMseqs2,
> FoldSeek/ProstT5, VFDB, AntiFam, Metapredict, dbCAN 5, MinCED/CRISPR, PhiSpy, geNomad, CheckV),
> the gapseq nucleotide route and the six gspa-nf default predictor containers were run genuinely
> on public inputs with pinned, hashed versions (ws and Shaheen); every raw claim equals a GSPA-free
> oracle, the three CLI routes agree byte-for-byte, each injected failure (missing tool, malformed
> row, truncated row, dropped protein, unknown protein, region outside/unknown contig) gives its
> exact typed outcome with zero claims, and the Nextflow trace records every process COMPLETED in
> its container with mutually consistent outputs. Eight engine defects were repaired (availability
> probes, strict parsers, unknown-protein rejection, VFDB/CRISPR registration, CRISPR attribution,
> Metapredict driver and coverage, placeholder-contig wording, CdsTranslator multi-copy loci), each
> with a regression spec, plus one pipeline defect (the merge container lacked `ps`, failing every
> traced run). Limits: SignalP 6 and DeepTMHMM cannot be obtained offline (parser fixtures only);
> PhiSpy is non-deterministic (oracle taken from the same run); CheckV misreports a complete
> bacterial chromosome as a viral contig (tool domain, passed through unchanged); the gapseq
> nucleotide route is certified as its documented typed failure, not as a claim producer; antiSMASH
> remains captured-not-fresh; neural/region/site sidecars and the JVM-only operon predictor are out
> of the bar's "tool outputs". Evidence: `docs/benchmarks/v3-adapter-certification.md`,
> `benchmark/v3/results/adapter-certification-v1/`.
