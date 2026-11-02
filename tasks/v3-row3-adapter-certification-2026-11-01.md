# Row 3 adapter certification plan (frozen 2026-11-01)

Bar (quoted from `tasks/v3-acceptance.md`, row 3): "Versioned genuine tool outputs and
complete parser/dispatch paths, including injected failures and coverage loss."

This plan is frozen before any check runs. Its SHA-256 is recorded in
`benchmark/v3/results/adapter-certification-v1/PLAN.sha256`. Later changes are appended
under "Dated deviations" and never edit the frozen text above that heading.

## 1. Adapter inventory (enumerated from the code, not from the audit)

Source of truth: `gspa-predictors/src/main/groovy/gspa/predictor/**` and the registration
code in `AnnotationPipeline.createAllPredictors()` (read 2026-11-01).

| # | Adapter name | Class | Kind | Reachable from `gspa annotate`? | External tool it runs |
|---|---|---|---|---|---|
| 1 | diamond | similarity.DiamondPredictor | tool | yes (`similarity.tool: diamond`) | diamond blastp |
| 2 | mmseqs2 | similarity.MMseqs2Predictor | tool | yes (`similarity.tool: mmseqs2`) | mmseqs (createdb/search/convertalis via shell script) |
| 3 | foldseek | structure.FoldSeekPredictor | tool | yes (`structure.enabled`) | foldseek easy-search (ProstT5 or structure input) |
| 4 | interproscan | domain.InterProScanPredictor | tool | yes (`domains.tool: interproscan`) | interproscan.sh |
| 5 | hmmer | domain.HmmerPredictor | tool | yes (`domains.tool: hmmer`) | hmmsearch (Pfam) |
| 6 | antifam | domain.AntiFamPredictor | tool | yes (`antifam.enabled`) | hmmsearch (AntiFam) |
| 7 | eggnog-mapper | function.EggNogMapperPredictor | tool | yes (`enable: [eggnog-mapper]`) | emapper.py |
| 8 | signalp | localization.SignalPPredictor | tool | yes (`localization.enabled` + `signalP`) | signalp (licensed) |
| 9 | deeptmhmm | localization.DeepTmhmmPredictor | tool | yes (`localization.enabled` + `tmhmm`) | deeptmhmm (biolib, licensed/remote) |
| 10 | metapredict | disorder.MetapredictPredictor | tool (python driver) | yes (`disorder.enabled`) | python3 + metapredict module |
| 11 | gapseq / gapsmith | pathway.GapseqPredictor | tool, genome-level | yes (`pathway.enabled`, `tool`, `target`) | gapsmith find / gapseq find (prot or nucl) |
| 12 | amrfinder | specialized.AmrFinderPredictor | tool | yes (`enable: [amrfinder]`) | amrfinder --plus |
| 13 | dbcan | specialized.DbCanPredictor | tool | yes (`enable: [dbcan]`) | run_dbcan (v4 positional / v5 CAZyme_annotation) |
| 14 | antismash | specialized.AntiSmashPredictor | tool, genome-level | yes (`enable: [antismash]`) | antismash |
| 15 | crispr | specialized.CrisprPredictor | tool, genome-level | **no** (never instantiated by the pipeline) | minced |
| 16 | vfdb | specialized.VfdbPredictor | tool | **no** (never instantiated by the pipeline) | diamond blastp vs VFDB |
| 17 | genomad | viral.GenomadPredictor | genomic sidecar | yes (`viral.genomad.enabled`) | run_genomic_predictors.py → genomad (native or .sif) |
| 18 | checkv | viral.CheckVPredictor | genomic sidecar | yes | run_genomic_predictors.py → checkv |
| 19 | phispy | viral.PhiSpyPredictor | genomic sidecar | yes | run_genomic_predictors.py → PhiSpy.py |
| 20 | prodigal | genecalling.ProdigalCaller | gene caller | yes | prodigal |
| 21 | pyrodigal | genecalling.PyrodigalCaller | gene caller | yes | pyrodigal (python) |
| 22 | operons | context.OperonPredictor | JVM only | yes | none |
| 23–33 | deepgoplus-esm2, proteinfer, clean, esm2-centroid, deepgo-plusplus, deepgo-plusplus-light, psortb, deepfri, mdf, deepec, deeparg | neural.* / localization.PSORTbPredictor / structure.* / specialized.DeepArgPredictor | neural sidecar | yes (per `neural.*` / `foss.*` flags) | run_neural_predictors.py / run_term_predictors.py + model weights |
| 34–36 | deepsig, tmbed, tppred3 | localization.* | region sidecar | yes (`localization.deepSig` etc.) | run_region_predictors.py + model weights |
| 37–38 | musitedeep, scannet | sites.* | site sidecar | yes (`foss.musiteDeep` etc.) | run_site_predictors.py + model weights |
| – | CheckM2Runner, GtdbTkRunner | taxonomy.* | not Predictor implementations | n/a | checkm2, gtdbtk |

The audit's count of "22 adapters" corresponds to rows 1–22. Rows 23–38 are sidecar
adapters whose "tool" is a GSPA Python sidecar plus third-party model weights; they are
listed so that nothing is silently left out, and their disposition is stated in §4.

## 2. Recipe (identical for every adapter certified in this round)

Follows `benchmark/v3/check_amrfinder_route.py` and `docs/benchmarks/v3-amrfinder-observed.md`.

1. **Pinned install.** A public release pinned by tag/version and verified by SHA-256 of the
   downloaded archive, unpacked into its own empty directory under
   `/mnt/data2/leechuck/gspa-cert/` (ws; `/tmp` and `$HOME` are not used: `/` is full).
   Container images are not available on ws (docker storage on the full root disk) and
   Shaheen's ControlMaster is down at plan time (see §5), so Shaheen/unimatrix bundles are
   used only if the master comes back.
2. **Genuine run** on public inputs: the two Swiss-Prot proteins already used (P18843 NadE,
   P0ADB7 EcnB) plus the synthetic 80-alanine control from
   `gspa-predictors/src/test/resources/observed-output/interproscan/input.faa`; for
   nucleotide tools the Price genome *E. coli* BW25113, RefSeq GCF_000750555.1 (FASTA + GFF
   from NCBI Datasets). Reference search databases are built from a 285-entry Swiss-Prot
   set (UniProt REST, `reviewed:true`, NadE/EcnA/EcnB families; TSV hashed) with GO/EC
   labels in the FASTA headers, exactly the header convention the adapters parse.
   Recorded: tool version string, command line, SHA-256 of every input, database file and
   output file, wall time. The tool output is stored under
   `gspa-predictors/src/test/resources/observed-output/<tool>/` with a `provenance.json`.
3. **Independent oracle.** The expected claims are derived from the stored native output
   by a Python reader that does not import GSPA code (identity/coverage/score fields
   recomputed from the native columns).
4. **Route checks** through the real CLI JAR (`gspa-cli/build/libs/gspa-3.0.0.jar`, built on
   ws after any JVM change, `./gradlew build` gate must pass):
   - `annotate` (raw claims via `--claims-out`), `annotate --integration-config` and direct
     `integrate` on the raw claims: raw claims equal the oracle, annotation-integrated and
     direct-integrated predictions byte-identical.
   - **missing executable** (configured `executablePath` absent, or `PATH` without the tool
     for adapters that resolve by name) → status `UNAVAILABLE`, zero claims.
   - **malformed output** (stub executable replays the genuine output with one row
     truncated / a non-numeric statistic / an unknown protein id) → `FAILED` with a
     diagnostic naming the problem, zero claims.
   - **header-only / empty output** (stub writes the documented empty result) → `SUCCESS_ZERO`
     for tools whose empty output is a genuine "no hit"; `FAILED` for tools that write one
     record per input protein (SignalP, DeepTMHMM, Metapredict driver).
   - **dropped-protein coverage loss**: stub output omits one input protein. Pass criterion
     depends on the tool's output model: tools with per-protein records must report
     `PARTIAL` listing the dropped id; search tools (DIAMOND, MMseqs2, FoldSeek, VFDB) have
     no per-protein record, so the check instead (a) verifies the status is `SUCCESS` with
     the documented caveat that absence of a hit is not evidence, and (b) verifies that
     output naming a protein absent from the input is rejected (`FAILED`), never assigned.
     Genome-level tools: a contig without sequence → `PARTIAL` naming the contig.
   - Every status read from `<stem>_predictor_status.json`; every expected outcome is an
     exact enum value. The CLI exit code stays 0 for recorded per-predictor failures (as
     documented for AMRFinder); the manifest is the contract.
5. **Nucleotide attribution** (genome-level tools): region claims are projected onto CDS by
   `AbstractGenomicRegionSidecarPredictor.projectRegions`. Check: for a real genome, an
   independent Python overlap of region coordinates against GFF CDS coordinates must equal
   the set of proteins carrying the region claim; a region overlapping no CDS must still
   appear in `genome.metadata["genomic_regions:<name>"]` (not be silently dropped); a region
   with coordinates outside the contig must make the run `PARTIAL`.
6. **Pass/fail per adapter.** PASS = all of: version recorded; genuine output stored and
   hashed; raw claims = oracle; three-route parity; the four injected cases give the exact
   typed outcome; coverage-loss criterion met; JVM gate green. FAIL = any check fails and is
   not repaired within this round (repairs are allowed, documented before/after, defaults
   unchanged). NOT CERTIFIABLE = tool cannot be obtained legally or offline (reason given).
   BLOCKED = infrastructure prevented the genuine run in this round (reason given).

## 3. Harness

`benchmark/v3/check_adapter_route.py` — one generic checker driven by a per-adapter JSON
spec (config block, fixture dir, oracle, stub variants, expected outcomes). It binds and
re-hashes every input after execution, writes `runs.json`, per-run logs and
`verification.json` (schema `gspa.adapter-route-replay/1`). Each adapter's spec and
results go to `benchmark/v3/results/adapter-certification-v1/<adapter>/`.

## 4. Per-adapter disposition for this round

| Adapter | Planned action | Pass criteria beyond §2 |
|---|---|---|
| diamond | Re-run under the generic harness with the 285-entry reference DB (cheap; uniform failure cases) | GO/EC claims equal the oracle for every hit row |
| mmseqs2 | Pinned MMseqs2 17-b804f; same reference DB; genuine search | As diamond; the shell-script route must honour `executablePath` or the plan records the defect and the repair |
| foldseek | Pinned FoldSeek 10-941cd33, ProstT5 weights via `foldseek databases ProstT5`, target DB from the reference FASTA (ProstT5 mode, CPU) | GO/EC claims equal the oracle; prostMode column layout verified |
| vfdb | VFDB `VFDB_setA_pro.fas.gz` (download date and hash recorded; VFDB publishes no version string), diamond makedb, genuine blastp | Adapter must become reachable (`predictors.enable: [vfdb]`, `predictors.vfdb.database`) or the plan records it as unreachable; claims equal the oracle |
| crispr | Pinned MinCED 0.4.2 on BW25113 | Adapter must become reachable (`predictors.enable: [crispr]`); arrays must not be silently dropped by `applyAnnotations` (contig-keyed results): expected repair = project arrays through `projectRegions` and keep whole regions in genome metadata; attribution oracle per §2.5 |
| metapredict | Pinned metapredict 3.1.1 (+ torch 2.5.1 CPU) in a uv venv; python driver as shipped | Per-protein driver exceptions must surface as `PARTIAL` with the protein id (currently "UNVERIFIED" in the source); if they do not, record defect and repair |
| dbcan | Pinned run_dbcan 5.2.8 + its database (if download completes within the round); interface version 5 | overview.tsv claims equal the oracle |
| antifam | HMMER 3.4 built from source + genuine AntiFam models (EBI release, version recorded) | Replaces "synthetic model only"; outcomes per §2 |
| phispy | Pinned PhiSpy 4.2.21 natively (no container) on BW25113 via the genomic sidecar | Region TSV, attribution per §2.5 |
| signalp | NOT CERTIFIABLE: SignalP 6 requires a DTU academic licence download; no binary on ws, unimatrix or Shaheen | Record the documented CLI (`signalp6 --fastafile/--organism/--format`) against the adapter's `--fasta/--org/--format`; parser checked only on documented format |
| deeptmhmm | NOT CERTIFIABLE: DeepTMHMM runs through pybiolib (remote DTU execution, account required); no offline binary | As signalp |
| gapseq (nucleotide route) | BLOCKED unless gapseq becomes runnable: not installed on ws (no conda), docker storage full, Shaheen master down. The protein route (gapsmith) is already certified | Record; dispatch-level typed failure for `target=genome` is covered by existing JVM tests, not by a genuine run |
| genomad, checkv | BLOCKED in this round: need singularity images or native installs with 1–2 GB databases; Shaheen down, ws cannot pull images | Record |
| antismash, interproscan, hmmer, eggnog-mapper, amrfinder, gapsmith, prodigal, pyrodigal | Previously certified (linked reports); not re-run | Carry forward; antiSMASH remains "captured, not fresh" |
| operons | JVM-only, no external tool | Out of scope of the bar's "tool outputs"; covered by Spock tests |
| neural/region/site sidecars (rows 23–38) | Out of scope of this round: model weights and GPUs are not part of the row-3 bar's "tool outputs"; CLEAN is certified separately (`v3-clean-contract.md`); others keep their controlled-fixture tests | Stated explicitly in the report |
| Nextflow `main.nf` predictor containers | BLOCKED: requires Singularity on Shaheen (master down) or docker on ws (storage full). `nextflow` parse check is run if a Nextflow binary can be placed on `/mnt/data2` | Record |

## 5. Constraints observed

- Shaheen: `ssh -O check shaheen` reports a master pid but every session is refused
  (`Session open refused by peer`; keepalive log "master down" since 07:17). Per the rules,
  Shaheen work stops here and the final report says an OTP is needed.
- ws root filesystem (`/`, incl. `$HOME` and `/var/lib/docker`) is at 100 %; all work goes to
  `/mnt/data2/leechuck/gspa-cert` (732 GB free) with `TMPDIR` and `UV_CACHE_DIR` set there.
- No edits to `tasks/v3-acceptance.md` or `README.md`; defaults unchanged; JVM gate after any
  JVM change; downloads treated as untrusted (own directories, no interpreter run from them).

## Dated deviations

(none at freeze time)

- 2026-11-01 11:12 — VFDB: the planned URL `/VFs/download/VFDB_setA_pro.fas.gz` returns 404; the
  file is served from `http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz` (sha256 recorded in
  `gspa-cert/downloads/SHA256SUMS`). Same file, different path.
- 2026-11-01 11:17 — dbCAN: `run_dbcan database` from bcb.unl.edu returned HTML pages instead of the
  database files (invalid tar archive); the `--aws_s3` mirror (`db_v5-2_9-13-2025`) was used and hashed.
- 2026-11-01 11:40 — The Shaheen ControlMaster came back during the round (§5 said "only if the master
  comes back"): geNomad (sif on /project, DB v1.9 from Zenodo) and CheckV (quay biocontainer 1.0.3,
  DB v1.5) were attempted on Shaheen through the genomic sidecar; PhiSpy 4.2.21 ran natively on ws
  after pinning biopython 1.81 (PhiSpy uses `SeqFeature.strand`, removed in biopython 1.82).
- 2026-11-02 — Pass criterion refinement for the dropped-protein case: when every genuine hit row
  belongs to the dropped protein (FoldSeek: all 276 hits are P18843; AntiFam: the seed member is the
  only hit) the remaining table is the documented empty result, so the exact expected outcome is
  `SUCCESS_ZERO`, not `SUCCESS`. Recorded in the specs with a note.
- 2026-11-02 — CRISPR attribution: the RefSeq `protein.faa` carries 6 WP_ accessions shared by several
  identical CDS copies; the pipeline cannot place them and parks them on a placeholder contig. The
  certification runs the CDS-translation route (no `--proteins`) so every protein has coordinates,
  and a second run with the RefSeq FAA expects `PARTIAL` naming those 6 proteins. The engine's wording
  for the placeholder contig ("contig has no nucleotide sequence") was repaired to name the
  unplaceable proteins (§2.6 repairs allowed; before/after runs kept).
- 2026-11-02 — Harness defects found and fixed during the round (reruns kept, failed runs renamed
  `*.failed-harness-bug`): the python-driver malformed/truncated variants altered the header row
  instead of a data row; the hmmsearch truncation only cut the free-text description; a mis-indented
  loop turned the region branch into a `for/else`; CDS-translated ids `<contig>:<id>[.n]` were not
  mapped back to GFF protein ids.
- 2026-11-02 (second pass, on the supervisor's instruction) — The open items were continued under the
  same plan: (a) geNomad and CheckV route harness rerun on Shaheen with the fixed harness (jobs
  15745762 → 15745876); (b) `gspa.io.CdsTranslator` repaired in gspa-core (CDS grouped per locus, not
  per protein_id; `CdsTranslatorSpec` regressions; JVM gate green) and the PhiSpy/CRISPR attribution
  checks rerun; (c) the gapseq nucleotide route certified on Shaheen with the native gapseq 2.1.0 env
  (`adapter-cert-gapseq.sbatch`, job 15745873) — its documented outcome with hits is a typed `FAILED`
  ("nucleotide pathway hits require coordinate-to-protein attribution"), zero claims; (d) gspa-nf
  `main.nf` default predictor containers (pyrodigal, diamond, hmmer, barrnap, minced, merge) run on
  Shaheen through Slurm + Singularity (`nf-shaheen.config`, `nf-shaheen-run.sh`) for BW25113.
- 2026-11-02 — Observation recorded, not a repair: CheckV reports the whole BW25113 chromosome as a
  complete high-quality viral contig; the adapter passes it through as designed.
- 2026-11-02 (third pass) — gapseq nucleotide route, harness job 15745883 (rc=1): no check failed; the
  harness crashed before its first check. It passes `annotate --threads 1` for every adapter, and
  `--threads` overrides `predictors.pathway.threads` (AnnotateCommand), so the JAR ran
  `gapseq find … -K 1`; the genome-wide search against 155,682 reference sequences was still running
  when the harness's hard 3600 s `subprocess` limit expired (the direct run with -K 32 took 449 s).
  Harness defect, engine as documented: the spec now carries `cli_threads` (32 for gapseq) and
  `timeout_seconds` (7200), the failed run is kept as `after.failed-harness-timeout`, and the
  harness-only job was resubmitted (15745997, `--exclude=nid00022,nid00497,nid01033`).
- 2026-11-02 — gspa-nf on Shaheen: the 03:29 run (nf-shaheen.attempt4.log) completed PYRODIGAL,
  DIAMOND_BLASTP, HMMSEARCH, BARRNAP and MINCED, but MERGE_ANNOTATIONS failed with exit 1 before its
  script ran: `python:3.12-slim` ships no `ps`, which Nextflow's task-metric collection requires
  whenever `-with-trace`/`-with-report` is on. `gspa-nf/slurm_singularity.config` (unimatrix) already
  overrode that process to `python:3.12` for this reason; the shipped default (`modules/quality.nf`,
  `nextflow.config`) did not. Repaired in both (`python:3.12`: same interpreter, full Debian image);
  `MDF` and `FLATTEN_PROPHAGE_CDS` keep the slim image and are outside this round. Second finding:
  Nextflow refused to overwrite the previous attempt's `trace.tsv`/`report.html` (`trace.overwrite`
  off), so that run's trace was never written and the stale attempt-3 trace stayed (its DIAMOND
  failure was the DIAMOND database not yet built — a setup error of the certification, fixed before
  the 03:29 run, not a pipeline defect). `nf-shaheen.config` now sets `trace.overwrite`,
  `report.overwrite` and adds `container` to the trace fields. The rerun is a fresh run instead of
  `-resume`: a resumed run records cached tasks as CACHED, and the evidence needs every process
  COMPLETED with its container in one trace (previous outputs kept as `nf-out.attempt4`).
- 2026-11-02 — Results kept in the repository: every file over 5 MB (the per-run BW25113 GFF3 exports,
  the gapseq Reactions.tbl replays) is moved to `/mnt/data2/leechuck/gspa-cert/large-files/` and
  listed in `LARGE_FILES_REMOVED.tsv`; the gapseq observed-output fixture stores a deterministic
  subset (first 500 data rows) plus the full file's sha256.
