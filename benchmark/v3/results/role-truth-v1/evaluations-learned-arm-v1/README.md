# Role-truth evaluations of the 21-genome learned-context arms (jobs 10658 and 10666, 2026-10-29)

Arms: `benchmark/v3/prepare_role_genome_context.py plan` outputs run under `/data/hohndor/gspa-v3/role-genome-arm-v1`
(JAR 10402, 10406 `diamond_homology` models, 18 cohort genomes + M. tuberculosis H37Rv GCF_000195955.2, E. coli K-12
GCF_000005845.2, B. subtilis 168 GCF_000009045.1). The `SEQONLY-*` arms are job 10658 (`arm.sbatch`); the nine
`UNIVERSAL-*`, `PATHWAY-*` and `ORGANISM-*` arms are job 10666 (`arm-v2.sbatch`, plans `plans-v2`), the rerun after
the layout wiring fix. Each folder: `model-genomes/` = evaluator on the concatenated three model-genome tables;
`model-genomes-plus-sso/` = plus the Saccharolobus solfataricus cohort table GCA_000007005.1 (in-cohort genome, so its
two ThrRS cases are not held out from the models); `batch_stats.tsv` = per-genome runtime/heap. The integrated tables
are not archived. Registry sha256 fa8e916a…, negative-only sha256 f2225d73…, evaluator `evaluate_role_truth.py
--format integrated_tsv --score-field score`, threshold 0.5. `SHA256SUMS` covers every file in this folder
(`sha256sum -c SHA256SUMS`).

"Candidates rescored" (model-genome candidates the learned model scored) is not carried by the evaluation folders and
is therefore not determined here; the job-10658 diagnosis (0 of 101,705 for both SEQONLY models, every candidate
`unseen_predictor_contract`) is the only measured value, and the identical model-genome scores in the table below
are consistent with the same abstention in the context arms. The called counts are over all 11 cases of the plus-SSO
evaluation; strict top-1 is over the scored cases.

| folder | learned model | candidates rescored | scored cases, model genomes | scored cases, plus SSO | strict top-1 (plus SSO) | mean margin, model genomes | mean margin, plus SSO | actives >= 0.5 | negatives >= 0.5 |
|---|---|---|---|---|---|---|---|---|---|
| `SEQONLY-none` | no learned model | not determined | 4 of 11 | 5 of 11 | 5 of 5 | 0.190 | 0.163 | 1 of 11 | 0 of 13 |
| `SEQONLY-logistic` | 10406 logistic | not determined | 4 of 11 | 5 of 11 | 5 of 5 | 0.190 | 0.163 | 1 of 11 | 0 of 13 |
| `SEQONLY-gbdt` | 10406 gbdt | not determined | 4 of 11 | 5 of 11 | 4 of 5 | 0.190 | 0.151 | 1 of 11 | 0 of 13 |
| `UNIVERSAL-none` | no learned model | not determined | 4 of 11 | 5 of 11 | 5 of 5 | 0.175 | 0.150 | 2 of 11 | 1 of 13 |
| `UNIVERSAL-logistic` | 10406 logistic | not determined | 4 of 11 | 5 of 11 | 5 of 5 | 0.175 | 0.153 | 2 of 11 | 1 of 13 |
| `UNIVERSAL-gbdt` | 10406 gbdt | not determined | 4 of 11 | 5 of 11 | 5 of 5 | 0.175 | 0.151 | 2 of 11 | 1 of 13 |
| `PATHWAY-none` | no learned model | not determined | 4 of 11 | 5 of 11 | 5 of 5 | 0.175 | 0.150 | 2 of 11 | 1 of 13 |
| `PATHWAY-logistic` | 10406 logistic | not determined | 4 of 11 | 5 of 11 | 5 of 5 | 0.175 | 0.153 | 2 of 11 | 1 of 13 |
| `PATHWAY-gbdt` | 10406 gbdt | not determined | 4 of 11 | 5 of 11 | 5 of 5 | 0.175 | 0.152 | 2 of 11 | 1 of 13 |
| `ORGANISM-none` | no learned model | not determined | 4 of 11 | 5 of 11 | 5 of 5 | 0.175 | 0.150 | 2 of 11 | 1 of 13 |
| `ORGANISM-logistic` | 10406 logistic | not determined | 4 of 11 | 5 of 11 | 5 of 5 | 0.175 | 0.153 | 2 of 11 | 1 of 13 |
| `ORGANISM-gbdt` | 10406 gbdt | not determined | 4 of 11 | 5 of 11 | 5 of 5 | 0.175 | 0.149 | 2 of 11 | 1 of 13 |

Per-case candidate scores of the `model-genomes-plus-sso` evaluation (the registry marks GuaB2, FolP1, GpmI, TalB and
the ThrRS catalytic subunit as the active members; Fsa (P78055) and TalA are not scored):

| case | member | SEQONLY-none | SEQONLY-logistic | SEQONLY-gbdt | UNIVERSAL-none | UNIVERSAL-logistic | UNIVERSAL-gbdt | PATHWAY-none | PATHWAY-logistic | PATHWAY-gbdt | ORGANISM-none | ORGANISM-logistic | ORGANISM-gbdt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GuaB (M. tuberculosis) | GuaB2 (active) | 0.2626 | 0.2626 | 0.2626 | 0.2626 | 0.2626 | 0.2626 | 0.2626 | 0.2626 | 0.2626 | 0.2626 | 0.2626 | 0.2626 |
|  | GuaB1 | 0.2066 | 0.2066 | 0.2066 | 0.2066 | 0.2066 | 0.2066 | 0.2066 | 0.2066 | 0.2066 | 0.2066 | 0.2066 | 0.2066 |
|  | GuaB3 | 0.1316 | 0.1316 | 0.1316 | 0.1316 | 0.1316 | 0.1316 | 0.1316 | 0.1316 | 0.1316 | 0.1316 | 0.1316 | 0.1316 |
| FolP (M. tuberculosis) | FolP1 (active) | 0.2499 | 0.2499 | 0.2499 | 0.5109 | 0.5109 | 0.5109 | 0.5109 | 0.5109 | 0.5109 | 0.5109 | 0.5109 | 0.5109 |
|  | FolP2 | 0.1836 | 0.1836 | 0.1836 | 0.1836 | 0.1836 | 0.1836 | 0.1836 | 0.1836 | 0.1836 | 0.1836 | 0.1836 | 0.1836 |
| PGM (B. subtilis) | GpmI (active) | 0.5400 | 0.5400 | 0.5400 | 0.5400 | 0.5400 | 0.5400 | 0.5400 | 0.5400 | 0.5400 | 0.5400 | 0.5400 | 0.5400 |
|  | YhfR (negative) | 0.1828 | 0.1828 | 0.1828 | 0.5047 | 0.5047 | 0.5047 | 0.5047 | 0.5047 | 0.5047 | 0.5047 | 0.5047 | 0.5047 |
| transaldolase (E. coli) | TalB (active) | 0.4180 | 0.4180 | 0.4180 | 0.4180 | 0.4180 | 0.4180 | 0.4180 | 0.4180 | 0.4180 | 0.4180 | 0.4180 | 0.4180 |
|  | TalC/FsaB | 0.1372 | 0.1372 | 0.1372 | 0.1372 | 0.1372 | 0.1372 | 0.1372 | 0.1372 | 0.1372 | 0.1372 | 0.1372 | 0.1372 |
| ThrRS (S. solfataricus; in-cohort) | catalytic subunit (active) | 0.2358 | 0.1536 | 0.1117 | 0.2358 | 0.1624 | 0.1295 | 0.2358 | 0.1616 | 0.1285 | 0.2358 | 0.1566 | 0.1259 |
|  | editing subunit | 0.1831 | 0.0966 | 0.1168 | 0.1831 | 0.0967 | 0.0732 | 0.1831 | 0.0952 | 0.0693 | 0.1831 | 0.0934 | 0.0784 |
