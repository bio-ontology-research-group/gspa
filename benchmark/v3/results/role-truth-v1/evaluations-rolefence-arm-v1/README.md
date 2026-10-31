# Role-truth evaluations of the 21-genome arms under the role-fenced models (WP2e; jobs 15690304 and 15690691, 2026-10-29)

## What WP2e is

The learned integration models of job 10406 were trained on cohort claims computed against the cohort-fenced
Swiss-Prot reference. The three model genomes (M. tuberculosis H37Rv, E. coli K-12 MG1655, B. subtilis 168) use the
role-fenced reference, a different predictor contract, so those models abstained on every model-genome candidate
(see `../evaluations-learned-arm-v1/README.md`). WP2e retrains the `diamond_homology` arm on cohort claims computed
against the role-fenced reference (arm name `diamond_homology_rolefence`) so that the 18 cohort genomes and the three
model genomes share one contract and the models can score the model genomes.

Training ran on Shaheen III as job 15690304 (`run_rolefence_arm.py` stages prepare, export, train, evaluate; JAR 10402;
four contexts SEQONLY, UNIVERSAL, PATHWAY, ORGANISM, each with a logistic and a gbdt model; 65 unit tests passed in
`training/job-15690304/tests.log`; exit code 0). `training/` holds the five receipt JSONs (`protocol.json`,
`export.json`, `training.json`, `evaluation.json`, `registration.json`) and the job receipts (stage logs, start/finish
stamps, resources, input hashes). The cohort-split evaluation in `evaluation.json` is small (305 to 310 supported test
rows per context and kind) and conditional on the exported candidates; it makes no calibration or promotion claim.

The supervision universe differs from that of job 10406 by design and the difference is disclosed in
`training.json` (`supervision_comparison`); it is identical in all four contexts:

| context | shared examples | original_only (10406 universe) | new_only (role-fenced) | representation changes |
|---|---|---|---|---|
| SEQONLY | 3,424 | 838 | 366 | 423 |
| UNIVERSAL | 3,424 | 838 | 366 | 423 |
| PATHWAY | 3,424 | 838 | 366 | 423 |
| ORGANISM | 3,424 | 838 | 366 | 423 |

## How the arms were run

Twelve arms (SEQONLY/UNIVERSAL/PATHWAY/ORGANISM x none/logistic/gbdt) ran as Shaheen job 15690691 (JAR 10402,
`integrate --batch-manifest`, plans `context/plans-rolefence`, `--learned-model` pointing at the WP2e model of the
matching context and kind, none for the `none` arms) on the 18 cohort genomes plus the three model genomes. All twelve
exited 0 (`receipts/pack-15690691.log`, `receipts/<CTX>-<kind>/{exit-code,task,command,resources,...}.txt`).
The supervisor evaluated the integrated tables with `evaluate_role_truth.py --format integrated_tsv --score-field
score`, threshold 0.5, registry sha256 fa8e916a..., negative-only sha256 f2225d73..., exactly as in the earlier
archive. Layout of each `<CTX>-<kind>/` folder: `model-genomes/` = evaluator on the concatenated three model-genome
tables; `model-genomes-plus-sso/` = plus the Saccharolobus solfataricus cohort table GCA_000007005.1 (in-cohort
genome: its two ThrRS cases are not held out); `batch_stats.tsv` = per-genome runtime and heap. The large integrated
tables are not archived. `SHA256SUMS` covers every file in this folder (`sha256sum -c SHA256SUMS`).

## Table (a): WP2e evaluation on the cohort test split

Supported rows only (`profile_lacks_both_classes_in_train_and_calibration` rows abstain; the row set therefore differs
slightly between contexts, so values are not paired comparisons).

| context | kind | AUROC | log-loss | supported rows | total rows | test positive rate | n |
|---|---|---|---|---|---|---|---|
| SEQONLY | logistic | 0.740 | 0.419 | 310 | 429 | 0.145 | 310 |
| SEQONLY | gbdt | 0.738 | 0.433 | 310 | 429 | 0.145 | 310 |
| UNIVERSAL | logistic | 0.776 | 0.404 | 309 | 429 | 0.146 | 309 |
| UNIVERSAL | gbdt | 0.669 | 0.403 | 309 | 429 | 0.146 | 309 |
| PATHWAY | logistic | 0.772 | 0.428 | 308 | 429 | 0.146 | 308 |
| PATHWAY | gbdt | 0.708 | 0.395 | 308 | 429 | 0.146 | 308 |
| ORGANISM | logistic | 0.785 | 0.430 | 305 | 429 | 0.144 | 305 |
| ORGANISM | gbdt | 0.663 | 0.407 | 305 | 429 | 0.144 | 305 |

## Table (b): role-truth summary per arm

Strict top-1 is over the scored cases. Actives are counted over all 11 cases, negatives (`case negatives`) over the 13
negative members of the plus-SSO evaluation; `negative_only_per_case.jsonl` of every evaluation reports 0 of 7
negative-only loci called in all twelve arms.

| arm | scored, model genomes | strict top-1, model genomes | mean margin, model genomes | scored, plus SSO | strict top-1, plus SSO | mean margin, plus SSO | actives >= 0.5 (plus SSO) | case negatives >= 0.5 (plus SSO) |
|---|---|---|---|---|---|---|---|---|
| `SEQONLY-none` | 4 of 11 | 4 of 4 | 0.190 | 5 of 11 | 5 of 5 | 0.162 | 1 of 11 | 0 |
| `SEQONLY-logistic` | 4 of 11 | 4 of 4 | 0.356 | 5 of 11 | 5 of 5 | 0.298 | 2 of 11 | 0 |
| `SEQONLY-gbdt` | 4 of 11 | 4 of 4 | 0.255 | 5 of 11 | 5 of 5 | 0.211 | 2 of 11 | 0 |
| `UNIVERSAL-none` | 4 of 11 | 4 of 4 | 0.175 | 5 of 11 | 5 of 5 | 0.150 | 2 of 11 | 1 |
| `UNIVERSAL-logistic` | 4 of 11 | 4 of 4 | 0.304 | 5 of 11 | 5 of 5 | 0.259 | 1 of 11 | 0 |
| `UNIVERSAL-gbdt` | 4 of 11 | 2 of 4 | 0.107 | 5 of 11 | 3 of 5 | 0.094 | 0 of 11 | 0 |
| `PATHWAY-none` | 4 of 11 | 4 of 4 | 0.175 | 5 of 11 | 5 of 5 | 0.150 | 2 of 11 | 1 |
| `PATHWAY-logistic` | 4 of 11 | 4 of 4 | 0.292 | 5 of 11 | 5 of 5 | 0.248 | 1 of 11 | 0 |
| `PATHWAY-gbdt` | 4 of 11 | 2 of 4 | 0.134 | 5 of 11 | 3 of 5 | 0.118 | 0 of 11 | 0 |
| `ORGANISM-none` | 4 of 11 | 4 of 4 | 0.175 | 5 of 11 | 5 of 5 | 0.150 | 2 of 11 | 1 |
| `ORGANISM-logistic` | 4 of 11 | 4 of 4 | 0.294 | 5 of 11 | 5 of 5 | 0.250 | 1 of 11 | 0 |
| `ORGANISM-gbdt` | 4 of 11 | 2 of 4 | 0.097 | 5 of 11 | 3 of 5 | 0.088 | 0 of 11 | 0 |

Learned-prior coverage in the three model-genome integrated tables (counted by the supervisor from the integrated
tables, which are not archived: rows whose `priors_fired` contains a `learned.` entry). All scored candidates of the
scored registry cases are among the rescored rows (their scores differ from the `none` arm); a minority of rows carries no learned prior:

| context | rows in the three model-genome tables | rows carrying a learned prior, logistic | rows carrying a learned prior, gbdt |
|---|---|---|---|
| SEQONLY | 101,705 | 67,291 | 67,291 |
| UNIVERSAL | 101,705 | 62,444 | 62,444 |
| PATHWAY | 101,705 | 62,088 | 62,088 |
| ORGANISM | 101,705 | 61,649 | 61,649 |

## Table (c): per-case candidate scores, `model-genomes-plus-sso` evaluation

The registry marks GuaB2 (P9WKI7), FolP1 (P9WND1), GpmI (P39773), TalB (P0A870) and the ThrRS catalytic subunit
(AAK42622.1) as the active members; YhfR is O07617. Fsa (P78055) is not scored in any arm, including `none`, and
TalA has unknown truth. `n/s` = not scored.

| case | member | SEQONLY-none | SEQONLY-logistic | SEQONLY-gbdt | UNIVERSAL-none | UNIVERSAL-logistic | UNIVERSAL-gbdt | PATHWAY-none | PATHWAY-logistic | PATHWAY-gbdt | ORGANISM-none | ORGANISM-logistic | ORGANISM-gbdt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| PMID15240874-threonylation | AAK40712.1 | 0.1831 | 0.1126 | 0.1211 | 0.1831 | 0.0823 | 0.0777 | 0.1831 | 0.0814 | 0.0742 | 0.1831 | 0.0808 | 0.0712 |
|  | AAK42622.1 (active) | 0.2320 | 0.1809 | 0.1539 | 0.2320 | 0.1615 | 0.1223 | 0.2320 | 0.1567 | 0.1298 | 0.2320 | 0.1564 | 0.1227 |
| guab-imp-dehydrogenase-mtb | P9WKI3 | 0.2066 | 0.1433 | 0.1071 | 0.2066 | 0.1302 | 0.0994 | 0.2066 | 0.1271 | 0.0962 | 0.2066 | 0.1266 | 0.0992 |
|  | P9WKI5 | 0.1316 | 0.0592 | 0.0901 | 0.1316 | 0.0446 | 0.0773 | 0.1316 | 0.0452 | 0.0678 | 0.1316 | 0.0445 | 0.0737 |
|  | P9WKI7 (active) | 0.2626 | 0.2325 | 0.1365 | 0.2626 | 0.1877 | 0.0790 | 0.2626 | 0.1817 | 0.0911 | 0.2626 | 0.1812 | 0.0886 |
| folp-dihydropteroate-synthase-mtb | P9WNC9 | 0.1836 | 0.1132 | 0.1211 | 0.1836 | 0.1000 | 0.0994 | 0.1836 | 0.0984 | 0.0962 | 0.1836 | 0.0978 | 0.0992 |
|  | P9WND1 (active) | 0.2499 | 0.2102 | 0.1539 | 0.5109 | 0.1598 | 0.0706 | 0.5109 | 0.1554 | 0.0847 | 0.5109 | 0.1546 | 0.0773 |
| pgm-yhfr-gpmi-bsubtilis | O07617 | 0.1828 | 0.1122 | 0.1211 | 0.5047 | 0.1090 | 0.0916 | 0.5047 | 0.1070 | 0.0917 | 0.5047 | 0.1064 | 0.0919 |
|  | P39773 (active) | 0.5400 | 0.8392 | 0.5747 | 0.5400 | 0.7728 | 0.3487 | 0.5400 | 0.7487 | 0.3963 | 0.5400 | 0.7538 | 0.3170 |
| transaldolase-talb-fsa-ecoli | P0A870 (active) | 0.4180 | 0.5735 | 0.6199 | 0.4180 | 0.4850 | 0.2830 | 0.4180 | 0.4625 | 0.3174 | 0.4180 | 0.4664 | 0.2601 |
|  | P32669 | 0.1372 | 0.0641 | 0.1159 | 0.1372 | 0.0487 | 0.0639 | 0.1372 | 0.0489 | 0.0712 | 0.1372 | 0.0484 | 0.0653 |
|  | P78055 | n/s | n/s | n/s | n/s | n/s | n/s | n/s | n/s | n/s | n/s | n/s | n/s |
