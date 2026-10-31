# InterProScan observed-output and parser repair evidence

The old parser interpreted HAMAP/ProSite profile scores as E-values, could emit
negative annotation scores, and converted malformed statistics into confident
annotations. The repaired parser distinguishes statistic kinds, preserves raw
values/method provenance, bounds the legacy transforms and rejects malformed
values. The pipeline honors threads and can disable remote precalculated lookup.
The remaining 0.9 match-only fallback and E-value transform are explicitly legacy
heuristics; neither this work nor its metadata establishes calibration.

`observed/` is actual InterProScan 5.72-103.0 output from Slurm 10095: four analyses,
two public proteins and a synthetic control. It retains the exact command, input,
stdout/stderr, TSV and pre/post identities for 3,951 files. The initial controller
query returned COMPLETED/exit 0. `capture-job-state.json` transcribes that response
and records that the controller entry expired before archive retrieval; it is not
presented as a raw controller receipt. The successful generator receipt and every
captured artifact remain available. `launch-capture/` freezes the submitted files.

`tests/interpro-observed-before.json` and `interpro-observed-after.json` are real
old/new JAR parser results on that observed TSV. Nine annotation identities are
preserved; all three negative profile-derived scores are fixed and all six
E-value-derived scores remain unchanged. Nineteen score/observed-output tests
and the pipeline regression pass. The initial focused run failed 17 cases.
Full JVM validation totals 1,046 tests, zero failures/errors, two skips. The first
full run skipped six additional CLI Python checks because of stale test-env paths;
the corrected CLI defaults and final CLI rerun execute those checks successfully.
The benchmark suite passes 528 tests, three optional skips, including eight live
Nextflow cases against the repaired JAR.

`corpus/` binds old/new parser output for seven historical full-proteome TSVs.
All 379,740 annotation identities survive; 33,758 old scores are negative, while
all repaired scores are in [0,1]. There are 40,273 changed scores, including
reported-zero E-values. Every ordinary E-value in (0,1] preserves its score.
The 167 MB full JSONL output stays in the recorded local cache; all its hashes and
every input hash were independently rechecked when archiving. Historical TSV
execution provenance is not certified here, and annotation counts are not
biological accuracy or independent protein counts.

`launch-route/` freezes actual annotation/direct-integration job 10096, on JAR
`c0def1656fc17170869b1d8f25b80e1bf433ecbc583794f63c27aff23b1bc821`.
The job completed with exit 0 (`slurm-10096.txt`). Both actual annotation arms
return SUCCESS with nine identical raw claims; annotation-plus-integration and
direct integration predictions are byte-identical. All dependencies were rehashed
after execution. `route/` and `route-job.log` preserve the receipt and outputs;
all 24 recorded output hashes were independently verified after archive transfer.
The large JAR is stored under the content-addressed local engine cache and the
recorded remote launch directory, not duplicated here.

`upstream/` preserves the primary-source TSV writer and revision/hash. `code/`
freezes the changed sources, regressions and reproduction scripts; `tests/` has
before/after failures and passing logs. `sha256.json` binds every other archive
file and is verified after writing. See `docs/benchmarks/v3-interpro-observed.md`
and `tasks/v3-acceptance.md` for remaining gates. No whole audit item is accepted.
