# Extracellular product experiments

This archive retains the original failures, source snapshots, Sonnet reviews,
numerical comparisons and actual CLI outputs for the
[secretion experiment repair](../../../../docs/benchmarks/v3-extracellular-products.md).

- `history/10416`: initial Python run; COBRA imports failed because its cache
  defaulted to a read-only home directory. The original omissions also reproduce.
- `jobs/10419`: corrected scratch-cache launcher. Nine new tests fail against the
  original engine; 168 experiment and 70 benchmark tests pass after repair.
  All 8 standard Sco-GEM product comparisons, 36 masked endpoints and 18 endpoint
  classifications agree with independent COBRApy/GLPK calculations.
- `jobs/10422`: the Java builder regression fails before repair; 55 focused checks
  and the 1,474-case full JVM suite pass (six optional skips). Four actual CLI runs
  verify the added evidence-guided condition; all 42 masked endpoint comparisons
  and nine product classifications agree independently.
- `source/`: separately frozen Python v1/v2 and Java-build sources, original
  builder, launchers and hash manifests. The source archives verify 260, 263 and
  798 files respectively. Their parent runtime archives remain on Unimatrix at
  the locations bound in the package manifests.
- `review/`: three completed read-only Sonnet reviews and explicit adjudication.

The new JAR remains at
`/data/hohndor/gspa-v3/extracellular-builder-v1/jobs/10422/gspa-3.0.0.jar`, with
SHA-256 `e08b6e180c36281c2d296b3cd67b03f4b3bd78cecf97d188537306439c9ea2ec`.
It is not duplicated here; its job-computed hash is retained. This is a validated
development build, not a release or promotion of a learned model.

The model and minimal medium come from the authors' Sco-GEM v1.2.1 release,
commit `63003f8b01a692dc2b73f2c9e0bbde871a3fe6f0`. See Sulheim et al. (2020),
*Enzyme-Constrained Models and Omics Analysis of Streptomyces coelicolor Reveal
Metabolic Changes that Enhance Heterologous Production*,
[DOI 10.1016/j.isci.2020.101525](https://doi.org/10.1016/j.isci.2020.101525).
The model license is retained in `source/Sco-GEM-LICENSE.md`; the full prior
acquisition provenance is in `../secondary-source-intake-v1/`.

Real-model runs use literal model product IDs and controlled gene assignments.
They establish numerical execution, not experimentally measured production,
unknown-function recovery or agent superiority. No observed product measurement
is passed as a flux constraint. Standard before/after Java menus, experiment
records and annotation bytes match exactly; masked hypotheses remain identical.
The associated CLI fixtures and callback-agent tests are explicitly synthetic.

`local-verification.json` records checks on the downloaded evidence;
`archive.sha256` covers every retained file other than itself. No whole acceptance
item or v3.0.0 release is signed off by this archive.
