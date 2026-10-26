# Development phenotype benchmark: 61 published auxotroph cases

This is a baseline evaluation of the production metabolic experiment engine at
`d4e9895`, **not a v3 release result or evidence of improved function prediction**.
The [full report](../../benchmark/v3/results/development-d4e9895/auxotroph61-report.json)
contains every case, model/input/code hashes, execution status and sensitivity
subsets. Slurm job 9727 ran on node002; the preceding job 9725 failed before
execution on node007 (signal 53). Results were written to a separate retry
directory, not substituted into the failed attempt.

## Data and configuration

Measured labels are from the 61 Keio strain/substrate experiments in
[Mohammadbeygi et al., PLOS Computational Biology](https://doi.org/10.1371/journal.pcbi.1014469),
supplementary tables S8/S9. There are 36 no-growth, 17 partial-growth and 8
growth observations. Preparation and raw-OD checks are described in
[the data registry](v3-metabolic-data.md). Labels never enter the simulation;
each deletion is compared with wild type under the same carbon-source condition.
This set was selected through published auxotrophy predictions, so it cannot
establish population recall or specificity. Known model GPRs are retained for
this phenotype test; this is not a masked gene-function recovery experiment.

The evaluated configuration uses the unregulated iML1515 model supplied in the
authors' repository, carbon uptake 10, oxygen uptake 15, and ATP maintenance 6.86
(flux units mmol/gDW/h). All other import routes are closed unless explicitly
declared in the medium. These settings follow the acquired repository medium
and model defaults. They are **not the paper's complete configuration**: Methods
2.1–2.3 specify oxygen uptake 20, ATP maintenance 8.39, and substrate-dependent
regulatory reaction closures. The paper/repository discrepancy was identified
after this run; the separate protocol comparisons below now include those closures. The current result must not be called a
reproduction of the paper or a comparison with its full method.

## Results

All 61 cases executed, and all 61 matched wild types could grow. For the 44
clearly binary observations, 25 predictions agreed (56.8%). Partial-growth
observations remain separate: binary FBA does not predict those OD-based classes.

| Observed phenotype | Predicted growth | Predicted no growth | Unassessed |
| --- | ---: | ---: | ---: |
| Growth | 1 | 7 | 0 |
| Partial growth | 5 | 12 | 0 |
| No growth | 12 | 24 | 0 |

Of 43 predicted auxotrophs, 24 had measured no growth (55.8%); the denominator
includes the 12 partial-growth cases. All 43 no-growth predictions were
infeasible under the declared ATP maintenance requirement. This is reported
explicitly and is not interpreted as evidence for any protein function. In the
53 cases whose published metrics reproduced from raw OD, 24 of 38 binary
observations agreed; 24 of 40 predicted auxotrophs were confirmed. The full report
also gives the subset with optimal knockout solutions only.

The present model is insufficient to justify reliable function recovery from
feasibility alone. Other organisms and independent function-masking experiments
remain necessary; no conditions are tuned against these 61 labels.

## Registered protocol comparison

Slurm job 9736 ran all three frozen protocols through the production engine at
`e8bcc07`. Every protocol executed all 61 cases, all matched wild types grew,
and no prediction was unassessed. The original legacy predictions reproduced
exactly. Every saved simulation file was checked against its report SHA-256.

| Protocol | Binary agreement (44 cases) | Confirmed / predicted no-growth | Full report |
| --- | ---: | ---: | --- |
| Repository defaults, unregulated | 25/44 (56.8%) | 24/43 | [All cases](../../benchmark/v3/results/development-e8bcc07/auxotroph61-legacy-report.json) |
| Manuscript oxygen 20 and ATPM 8.39, unregulated | 25/44 (56.8%) | 24/43 | [All cases](../../benchmark/v3/results/development-e8bcc07/auxotroph61-manuscript_unregulated-report.json) |
| Manuscript conditions plus S2 regulatory closures | 26/44 (59.1%) | 25/44 | [All cases](../../benchmark/v3/results/development-e8bcc07/auxotroph61-manuscript_regulatory-report.json) |

The sole binary change was `bglX` on alpha-D-lactose: predicted growth became
predicted no growth, agreeing with the measured label. The 17 partial-growth cases
remain separate and unchanged (5 predicted growth, 12 predicted no growth).
The regulatory result still misses 11 measured no-growth cases and incorrectly
predicts no growth for 7 of the 8 measured growth cases. It does not support a
large improvement claim or reliable gene attribution.

Closures apply when their registered substrate has uptake enabled, to both
wild type and mutant; they remain separate from deleted-gene reactions and cannot
erase maintenance requirements. The raw source discrepancy remains: the S2 closure
sets were derived at oxygen 15 and default ATPM, then applied by the authors at
oxygen 20 / ATPM 8.39. This is a documented protocol comparison, not reproduction
of every original analysis or an independent test of the authors' model fitting.

## Independent solver replay

COBRApy 0.32.1 with GLPK independently replayed the same model, medium, GPR
deletions and registered closures (job 9738). All 61 cases agreed in each of
the three protocols: 183 matched wild-type/deletion comparisons, or 366 LPs.
The checker uses COBRApy's own GPR knockout operation and a different LP solver;
it does not call the production LP builder. Each oracle report binds the model,
case manifest, checker script and production report by SHA-256. Production report
hashes were checked against the committed artifacts before publication.

Full replays: [legacy](../../benchmark/v3/results/development-e8bcc07/oracle-legacy.json),
[manuscript unregulated](../../benchmark/v3/results/development-e8bcc07/oracle-manuscript_unregulated.json),
[manuscript regulatory](../../benchmark/v3/results/development-e8bcc07/oracle-manuscript_regulatory.json).
Agreement establishes implementation consistency for these configurations.
It does not improve the measured phenotype agreement above or validate gene
function assignments.

## Reproduction

```bash
python benchmark/v3/run_metabolic_benchmark.py \
  --cases /path/to/prepared/derived/cases.jsonl \
  --model /path/to/prepared/raw/iML1515.xml \
  --out /path/to/new-run-directory
```

The runner refuses an existing output directory or a model with a different
registered hash. Missing, failed and budget-skipped simulations remain unknown
and stay in the report. Tests verify phenotype differences between an essential
carbon route and redundant isozymes, unchanged simulations after label changes,
gene/medium identity, missing results and explicit partial-growth handling.
