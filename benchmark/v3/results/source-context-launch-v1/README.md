# Full source/context launch checkpoint

Job 10406 is queued behind successful independent CLEAN audit 10368 and corrected
real-genome probe 10403. This archive records the frozen source and launcher,
submission and scheduler confirmation, seven passing prerequisite-sealing tests,
and Sonnet's complete read-only review. All 37 probe-shared source files are
identical within the 60-file source bundle. No full-comparison result is included.

Before inference, `seal.py` checks the frozen inputs and both jobs' actual exit
files, then writes prerequisite hashes without overwriting an existing seal.
The existing full runner rechecks actual probe outputs, shared sources, native
CLEAN inputs and original training/runtime data before fitting. It requires all
20 exports, 40 frozen models, 720 genome replays and 61 EC analysis arms before
writing its final receipt. No model promotion is requested.

`static.sha256` references remote inputs, including the large parent archive,
parent supplement and tested JAR. Those large files are not included here. The
45.9 MB parent supplement is described in `parent-supplement.json`; its original
registered inputs passed runtime preflight 10400. `packages.sha256` and
`sealed-prerequisites.json` do not exist yet: job 10406 creates them only after
the completed prerequisite artifacts can be hashed, before any model fitting.
The full runner still checks their semantic content rather than trusting a hash
or scheduler exit alone.

To verify the files stored here, run `sha256sum -c archive.sha256` from this
directory. This does not verify external inputs or certify completed execution,
biological accuracy, or any of the twelve whole acceptance requirements.
