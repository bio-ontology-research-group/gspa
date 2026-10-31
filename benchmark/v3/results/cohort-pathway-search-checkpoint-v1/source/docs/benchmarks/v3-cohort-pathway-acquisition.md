# Observed pathway searches for the full context cohort

This acquisition extends the verified single-genome search to all **18 previously
resolved cohort assemblies**, comprising **49,463 proteins**, eight bacterial
genomes and ten archaeal genomes. The existing cohort membership and sequence
hashes are fixed before searching; no function labels are consulted to choose
genomes or reactions. This is a prerequisite for the missing context comparisons,
not their completion.

The [checkpoint evidence](../../benchmark/v3/results/cohort-pathway-search-checkpoint-v1/README.md)
retains frozen inputs and code, all first-array attempts, raw completed outputs,
independent audit results, regression logs and scheduler snapshots.

Status: **in progress**. Slurm array **10112** has finished: all eight bacterial
searches pass the independent audit, retaining **147,818 reaction rows** and
**15,384 pathway rows**; all ten archaeal attempts failed reference preflight.
The archaeal reference acquisition job is **10120**, with retry array **10131**
waiting on its successful completion. Search arrays run at most four tasks at
once, with four CPUs and 24 GiB per search. Current scheduler/receipt
state takes precedence over these checkpoint counts; do not restart live jobs
because a local observation call times out.

## Frozen execution

`observe_gapsmith.py` now accepts an explicit proteome, taxonomy and thread budget.
`observe_gapsmith_cohort.py` binds each array index to a unique genome, checks the
manifest, source, tool and proteome hashes before and after execution, and refuses
to overwrite an existing attempt. Every failure gets a failed execution receipt.
Each search copies its reference inputs to node-local Slurm scratch, checks them
before and after use, and publishes full output tables, resource usage and logs.

The original manifest SHA-256 is
`b9491693f7b34b29abe99c7f077c44f2b70b3003db258816d487d3badf2e8fa9`.
It records the exact acquired sequence hashes and the original selection,
provenance and split-file hashes. No split labels enter the search command.
The search remains gapsmith **0.1.0**, SHA-256
`8bb3cddd832092937533af3672d918e985ec1c0176ebfa4ad511a267a9e4bd5b`,
with DIAMOND **2.0.15**, SHA-256
`edf54a3ed21e0c2348b916d3c794f794097176c30e4fb5c32bc3ec3cfde6e6a2`.
The initial submission was rejected for using the SSH alias `unimatrix01` in
Slurm's exclusion list; its actual Slurm node name is `unimatrix`. No job was
created by that rejected submission. Only the node-exclusion spelling changed.

## Archaeal reference failure and repair

Genuine execution exposed two issues hidden by the first synthetic fixture:

1. The installed `seq/Archaea` directory contains only two user FASTA files.
   It lacks the `rev`, `unrev` and `rxn` reference categories. Its presence is not
   evidence that the reference database is available.
2. The installed gapsmith sequence-root resolver requires a `Bacteria` directory
   marker even when `-t Archaea` is requested. Copying only the archaeal tree
   therefore fails before searching.

The recorder now rejects missing reference categories and creates the required
empty marker for archaeal staging. No bacterial sequences are substituted into
an archaeal search. The corrected synthetic fixture reproduces the installed
resolver requirement; the failing-before and passing-after logs are retained.
These changes apply to a new frozen retry, leaving the active bacterial searches
and their source files unchanged.

`acquire_gapsmith_archaea.py` acquires the binary's pinned reference release,
[Zenodo record 16908828](https://zenodo.org/records/16908828), into a separate
directory. The record supplies the 16.1 MB archaeal archive and checksums.
The script checks the published outer MD5 and all three nested archive MD5s,
rejects archive paths/links that escape extraction, preserves the existing
curation tables and two user references, and records SHA-256 for every file.
It rechecks the published copy. Existing shared references are never modified.
Successful download alone does not certify completed acquisition.

The archive has downloaded and expanded to **46,477 reference files**. Publishing
those small files to shared storage is slow; the existing acquisition job's time
limit was extended from 30 minutes to three hours. It was not restarted. Retry
array 10131 has an `afterok:10120` dependency and cannot start before successful
publication and verification. Its manifest SHA-256 is
`00f67553b10f89fc288e28eb9f113e3e569163813d433728cf2b4297b21916f5`.

## Independent audit and tests

`audit_cohort_pathways.py` reads the frozen manifests and each execution receipt,
checks every output hash and recomputes table counts from the raw files. It
requires the declared taxonomy and thread budget, all reference categories,
identical reference bytes within each taxonomy and identical curation tables
across genomes. Missing and failed attempts remain visible. It cannot report
`verified` until every member of the original cohort has a completed search.

Eighteen new benchmark cases cover arbitrary genome names, both taxonomies,
allocation bounds, altered inputs and code, immutable attempts, failed execution,
the actual directory-marker failure, missing reference categories, audit tampering
and archive extraction boundaries. The full benchmark gate passes **591 tests**,
with twelve optional skips. Production JVM code and the verified engine remain
unchanged by this acquisition work.

To audit downloaded attempt directories:

```bash
python benchmark/v3/audit_cohort_pathways.py \
  --root "$HOME/.cache/gspa-v3-data/cohort-pathway-search-v1" \
  --out "$TMPDIR/cohort-pathway-audit.json" --require-complete
```

Add a second `--root` for a verified retry directory. Partial inspection may omit
`--require-complete`; its report still says `incomplete` and lists every missing
genome. An incomplete report must not certify the final ablation cohort.

## Remaining work and limits

Finish and verify the searches, replay their complete tables through production
gap construction, then run the matched pathway-only, universal-network and
organism-specific comparisons with the fixed cohort splits. Preserve the raw
candidate universe and evaluate suggestions separately from functional TSVs.
Existing models do not establish calibrated behavior for these new contexts.

The acquired reference databases have not yet been audited for overlap with
target proteins or held-out sequence families. Searches use no benchmark labels,
but that does not prove absence of reference or pretrained-predictor leakage.
Gap hypotheses and tool-predicted pathways are not biological truth. Independent
enzyme-locus/role/null labels, calibration and the other acceptance requirements
remain open; no whole item is accepted.
