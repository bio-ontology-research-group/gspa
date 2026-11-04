# Off-growth second round: tier A (own CLEAN tables), tier B and tier C runs (2026-11-03)

Frozen plan: `tasks/v3-row12-offgrowth-improvement-2026-11-02.md`. No rule, constant, metric or engine code was changed. This report records what ran and the counts; it does not quote recovery, precision or harm values.

## 1. CLEAN candidate tables

The 47 CLEAN sidecar outputs (`clean/out/<tag>/<tag>.clean.tsv`) were converted with the repository function `benchmark/v3/offgrowth2/candidates.py: clean_sidecar_to_candidates` (driver `candidates-clean/convert.py`, run with `PYTHONPATH=benchmark/v3 /home/leechuck/tmp/gspa-v3-testenv/bin/python`). Output: `candidates-clean/<tag>.candidates.tsv` and `<tag>.candidates.receipt.json`, plus `candidates-clean/candidates.sha256` (sha256sum format, bare file names). Per-tag checks are in `candidates-clean/checks.tsv` (script `candidates-clean/check.py`).

| Check | Result |
|---|---|
| Candidate tables written | 47 |
| Function input rows equal sidecar receipt `emitted_calls` | 47 of 47 |
| Candidate rows + dropped rows equal input rows | 47 of 47 |
| Scores in (0, 1] | 47 of 47 |
| Protein ids all present among FASTA first tokens (`clean/inputs/all/<tag>.faa`) | 47 of 47 (0 ids missing) |
| Shaheen copy verified with `sha256sum -c` | yes (`/project/k10018/gspa-v3/offgrowth2-inputs-v1/candidates-clean/`) |

**Row counts are not equal to the sidecar EC-call counts.** The repository function drops two kinds of rows by its own rule (docstring: skip `score == 0`, keep only complete 4-level EC numbers):

| Rows | Count (all 47 tags) |
|---|---|
| Sidecar EC calls (input rows) | 260,519 |
| Dropped: score exactly 0 (CLEAN prints `0.0`) | 68,326 |
| Dropped: incomplete / preliminary EC (`n`-numbers such as `4.2.2.n1`) | 2,242 |
| Dropped: non-EC rows, duplicate pairs | 0, 0 |
| Written candidate rows | 189,951 |

The supervisor should confirm this is the intended behaviour; I did not change the function.

## 2. Shaheen jobs

Submitted from `login4`, from `/project/k10018/code/gspa-v3`, with `benchmark/v3/shaheen/offgrowth2-genome.sbatch --submit` (TIER/KIND/WALLTIME/GENOMES exactly as specified; one exclusive 192-core node per genome; nid00022, nid00497, nid01033 excluded). `kuq` before submission: /scratch 313,595 of 1,024,000 files; /project 799,558 of 1,000,000 files. All 47 jobs started at 10:14 and all COMPLETED; nothing was written to `failed/`; no resubmission.

Code check: the job's code path (`offgrowth2/runner.py`, `engine.py`, `medium.py`, `offgrowth/*`, `followups/*`, the sbatch) on Shaheen is byte-identical to the ws working tree and to what the tier A six jobs recorded in their `code.sha256`. The Shaheen copies of the scoring-side files (`report.py`, `rules.py`, `score.py`, `score_driver.py`, `panels.py`, tests) are older than ws; the jobs do not import them, and scoring ran on ws.

Observation: at start-up every job spent several minutes in the Python import step (`import cobra, optlang`, process waiting on Lustre I/O in `/project/k10018/sw/envs/gspa-v3`); the same import on the login node took 4 min 30 s while the 47 jobs were starting. All jobs went on to complete.

| Tier | Jobs | Completed | Node-hours (sacct Elapsed x 1 node) | Longest job (h) |
|---|---|---|---|---|
| A (own CLEAN) | 22 | 22 | 26.40 | 2.57 |
| B | 8 | 8 | 3.71 | 0.69 |
| C | 17 | 17 | 11.91 | 1.38 |
| Total | 47 | 47 | 42.02 (8,068 core-hours at 192 cores/node) | |

Tars: `/scratch/project/k10018/gspa-v3/out/offgrowth2-v1/genomes/<genome>.clean.outputs.tar` (+ `.sha256`), copied to `~/tmp/gspa-v3-offgrowth2/tier{A,B,C}-clean/` and verified against the Shaheen sha256 files (47 of 47). Raw sacct lines: `tiers-bc/sacct.tsv`.

## 3. Cases per tier

From each genome's `genome_results2.json.gz` (`runs[].case.case_source` for masked runs; `skipped_cases`). N2 fallback = N2 rows in `cost.tsv` with `measured_panel_fallback = 1`. Script: `tiers-bc/count_cases.py`.

| Tier | Genomes | Unmasked runs | Masked: Price-148 | Masked: Swiss-Prot | Masked: Fitness Browser | Masked total | Skipped cases | N2 rows, masked (fallback 1 / 0) | N2 rows, unmasked (fallback 1 / 0) |
|---|---|---|---|---|---|---|---|---|---|
| A (own CLEAN) | 22 | 22 | 148 | 32 | 643 | 823 | 0 | 0 / 823 | 0 / 22 |
| B | 8 | 8 | 0 | 75 | 0 | 75 | 0 | 75 / 0 | 8 / 0 |
| C | 17 | 17 | 0 | 0 | 318 | 318 | 0 | 1 / 317 | 1 / 16 |

Notes on counts (no interpretation):
- Tier A own-CLEAN has the same case counts by source as tier A six (148 / 32 / 643).
- Tier B: one masked case per labelled protein in `swissprot_experimental_ecs_cohort8.tsv` (79 label rows, 75 distinct proteins). GCA_000009985.1 has no rows in that label file, so it has 0 masked cases and only its unmasked run. Every tier B N2 run used the measured-panel fallback (no measured panel or org map is given for tier B).
- Tier C: 318 masked cases = 318 label rows in `fitness_browser_ec_reannotations_tierc.tsv`. The only N2 fallback runs in tier C are Xantho (1 masked, 1 unmasked).
- No job died at case loading in tier B or C.

Per-genome detail (job id, wall time, node, cases by source, skipped, N2 fallback):

| Tier | Genome | Job id | State | Wall time | Node | Masked: Price-148 | Masked: Swiss-Prot | Masked: Fitness Browser | Skipped cases | N2 fallback runs (masked / unmasked) |
|---|---|---|---|---|---|---|---|---|---|---|
| A | GCA_000006965.1 | 15782136 | COMPLETED | 02:34:23 | nid00171 | 18 | 13 | 51 | 0 | 0 / 0 |
| A | GCA_003443575.1 | 15782137 | COMPLETED | 01:28:48 | nid00194 | 9 | 0 | 38 | 0 | 0 / 0 |
| A | GCA_004402975.1 | 15782138 | COMPLETED | 00:17:05 | nid00202 | 1 | 0 | 0 | 0 | 0 / 0 |
| A | GCF_000015245.1 | 15782139 | COMPLETED | 01:02:18 | nid00225 | 7 | 0 | 28 | 0 | 0 / 0 |
| A | GCF_000016065.1 | 15782140 | COMPLETED | 00:23:16 | nid00268 | 1 | 1 | 4 | 0 | 0 / 0 |
| A | GCF_000018145.1 | 15782141 | COMPLETED | 00:49:15 | nid00270 | 2 | 1 | 21 | 0 | 0 / 0 |
| A | GCF_000020125.1 | 15782142 | COMPLETED | 01:38:50 | nid00271 | 11 | 0 | 43 | 0 | 0 / 0 |
| A | GCF_000146165.2 | 15782143 | COMPLETED | 01:13:08 | nid00280 | 1 | 12 | 27 | 0 | 0 / 0 |
| A | GCF_000166295.1 | 15782144 | COMPLETED | 00:37:16 | nid00286 | 2 | 1 | 13 | 0 | 0 / 0 |
| A | GCF_000203935.1 | 15782145 | COMPLETED | 00:56:58 | nid00288 | 4 | 4 | 21 | 0 | 0 / 0 |
| A | GCF_000236665.1 | 15782146 | COMPLETED | 00:20:43 | nid00292 | 1 | 0 | 3 | 0 | 0 / 0 |
| A | GCF_000237365.1 | 15782147 | COMPLETED | 00:24:55 | nid00382 | 1 | 0 | 6 | 0 | 0 / 0 |
| A | GCF_000327065.1 | 15782148 | COMPLETED | 01:20:17 | nid00439 | 10 | 0 | 34 | 0 | 0 / 0 |
| A | GCF_000698265.1 | 15782149 | COMPLETED | 02:00:46 | nid00443 | 17 | 0 | 51 | 0 | 0 / 0 |
| A | GCF_000832305.1 | 15782150 | COMPLETED | 01:02:18 | nid00444 | 5 | 0 | 26 | 0 | 0 / 0 |
| A | GCF_001307155.1 | 15782151 | COMPLETED | 02:01:12 | nid00447 | 13 | 0 | 56 | 0 | 0 / 0 |
| A | GCF_001307275.1 | 15782152 | COMPLETED | 02:02:56 | nid00448 | 16 | 0 | 53 | 0 | 0 / 0 |
| A | GCF_001440325.1 | 15782153 | COMPLETED | 01:46:27 | nid00572 | 6 | 0 | 53 | 0 | 0 / 0 |
| A | GCF_001623525.1 | 15782154 | COMPLETED | 01:45:33 | nid00575 | 10 | 0 | 49 | 0 | 0 / 0 |
| A | GCF_001625455.1 | 15782155 | COMPLETED | 00:55:48 | nid00576 | 6 | 0 | 23 | 0 | 0 / 0 |
| A | GCF_001633105.1 | 15782156 | COMPLETED | 00:54:10 | nid00879 | 3 | 0 | 24 | 0 | 0 / 0 |
| A | GCF_002797435.1 | 15782157 | COMPLETED | 00:47:53 | nid00880 | 4 | 0 | 19 | 0 | 0 / 0 |
| B | GCA_000008485.1 | 15782158 | COMPLETED | 00:20:37 | nid00984 | 0 | 4 | 0 | 0 | 4 / 1 |
| B | GCA_000008765.1 | 15782159 | COMPLETED | 00:31:24 | nid00985 | 0 | 12 | 0 | 0 | 12 / 1 |
| B | GCA_000009085.1 | 15782160 | COMPLETED | 00:33:15 | nid00986 | 0 | 14 | 0 | 0 | 14 / 1 |
| B | GCA_000009985.1 | 15782161 | COMPLETED | 00:16:06 | nid00987 | 0 | 0 | 0 | 0 | 0 / 1 |
| B | GCA_000012525.1 | 15782162 | COMPLETED | 00:25:52 | nid00989 | 0 | 8 | 0 | 0 | 8 / 1 |
| B | GCA_000027305.1 | 15782163 | COMPLETED | 00:41:37 | nid00990 | 0 | 20 | 0 | 0 | 20 / 1 |
| B | GCA_000046845.1 | 15782164 | COMPLETED | 00:30:59 | nid01024 | 0 | 12 | 0 | 0 | 12 / 1 |
| B | GCA_000147055.1 | 15782165 | COMPLETED | 00:22:31 | nid01120 | 0 | 5 | 0 | 0 | 5 / 1 |
| C | Bifido | 15782179 | COMPLETED | 00:21:25 | nid01134 | 0 | 0 | 5 | 0 | 0 / 0 |
| C | Btheta | 15782176 | COMPLETED | 01:10:04 | nid01131 | 0 | 0 | 38 | 0 | 0 / 0 |
| C | Burk376 | 15782175 | COMPLETED | 00:29:22 | nid01130 | 0 | 0 | 10 | 0 | 0 / 0 |
| C | Caulo | 15782171 | COMPLETED | 00:50:37 | nid01126 | 0 | 0 | 25 | 0 | 0 / 0 |
| C | Cola | 15782172 | COMPLETED | 01:17:50 | nid01127 | 0 | 0 | 43 | 0 | 0 / 0 |
| C | DvH | 15782177 | COMPLETED | 00:33:42 | nid01132 | 0 | 0 | 14 | 0 | 0 / 0 |
| C | Dyella79 | 15782173 | COMPLETED | 00:39:03 | nid01128 | 0 | 0 | 18 | 0 | 0 / 0 |
| C | HerbieS | 15782166 | COMPLETED | 01:10:39 | nid01121 | 0 | 0 | 38 | 0 | 0 / 0 |
| C | Keio | 15782167 | COMPLETED | 00:23:15 | nid01122 | 0 | 0 | 6 | 0 | 0 / 0 |
| C | Koxy | 15782169 | COMPLETED | 01:14:48 | nid01124 | 0 | 0 | 40 | 0 | 0 / 0 |
| C | Lysobacter_OAE881 | 15782182 | COMPLETED | 00:16:48 | nid01137 | 0 | 0 | 1 | 0 | 0 / 0 |
| C | Miya | 15782168 | COMPLETED | 00:30:01 | nid01123 | 0 | 0 | 11 | 0 | 0 / 0 |
| C | Pedo557 | 15782174 | COMPLETED | 01:23:02 | nid01129 | 0 | 0 | 45 | 0 | 0 / 0 |
| C | Putida | 15782170 | COMPLETED | 00:42:17 | nid01125 | 0 | 0 | 20 | 0 | 0 / 0 |
| C | SynE | 15782178 | COMPLETED | 00:16:49 | nid01133 | 0 | 0 | 1 | 0 | 0 / 0 |
| C | Xantho | 15782181 | COMPLETED | 00:16:53 | nid01136 | 0 | 0 | 1 | 0 | 1 / 1 |
| C | rhodanobacter_10B01 | 15782180 | COMPLETED | 00:18:05 | nid01135 | 0 | 0 | 2 | 0 | 0 / 0 |

## 4. Scoring on ws

All three ran `benchmark/v3/offgrowth2/score_driver.py` from `/mnt/data2/leechuck/gspa` (working tree of `development/v3.0.0`) with `OFFGROWTH2_SCORE_JOBS=8`, `TMPDIR=~/tmp`, strict sha256 checks on (no refusals). No `--q-fixed-threshold` was given (as in the tier A six run), so each run searched its own threshold. Each command is also saved as `score-tier{A,B,C}-clean/command.sh`; driver output in `driver.log`, `/usr/bin/time -v` in `time_v.txt`.

Tier A, own CLEAN:
```
cd /mnt/data2/leechuck/gspa
D=benchmark/v3/offgrowth2/data; T=~/tmp/gspa-v3-offgrowth2/score-tierA-clean
OFFGROWTH2_SCORE_JOBS=8 PYTHONPATH=benchmark/v3 /home/leechuck/tmp/gspa-v3-testenv/bin/python \
  -m offgrowth2.score_driver --tars-dir ~/tmp/gspa-v3-offgrowth2/tierA-clean --inputs-dir ~/tmp/gspa-v3-offgrowth2 \
  --staging-dir $T/staging --out-dir $T/out --kind clean --tier A \
  --labels $D/price148_labels.tsv benchmark/v3/offgrowth/data/swissprot_experimental_ecs_price22.tsv $D/fitness_browser_ec_reannotations_price22.tsv \
  --threshold-labels $D/price148_labels.tsv benchmark/v3/offgrowth/data/swissprot_experimental_ecs_price22.tsv \
  --test-labels $D/fitness_browser_ec_reannotations_price22.tsv \
  --group-tsv $D/tier_a_groups.tsv \
  --matches-tsv benchmark/v3/results/stage2b-followups-v1/experiment_panel_matches.tsv \
  --org-map-tsv benchmark/v3/results/stage2b-followups-v1/price22_fitness_browser_organisms.tsv \
  --bootstrap-resamples 2000
```

Tier B:
```
cd /mnt/data2/leechuck/gspa
D=benchmark/v3/offgrowth2/data; T=~/tmp/gspa-v3-offgrowth2/score-tierB-clean
OFFGROWTH2_SCORE_JOBS=8 PYTHONPATH=benchmark/v3 /home/leechuck/tmp/gspa-v3-testenv/bin/python \
  -m offgrowth2.score_driver --tars-dir ~/tmp/gspa-v3-offgrowth2/tierB-clean --inputs-dir ~/tmp/gspa-v3-offgrowth2 \
  --staging-dir $T/staging --out-dir $T/out --kind clean --tier B \
  --labels $D/swissprot_experimental_ecs_cohort8.tsv \
  --threshold-labels $D/swissprot_experimental_ecs_cohort8.tsv \
  --test-labels $D/swissprot_experimental_ecs_cohort8.tsv \
  --bootstrap-resamples 2000
```

Tier C:
```
cd /mnt/data2/leechuck/gspa
D=benchmark/v3/offgrowth2/data; T=~/tmp/gspa-v3-offgrowth2/score-tierC-clean
OFFGROWTH2_SCORE_JOBS=8 PYTHONPATH=benchmark/v3 /home/leechuck/tmp/gspa-v3-testenv/bin/python \
  -m offgrowth2.score_driver --tars-dir ~/tmp/gspa-v3-offgrowth2/tierC-clean --inputs-dir ~/tmp/gspa-v3-offgrowth2 \
  --staging-dir $T/staging --out-dir $T/out --kind clean --tier C \
  --labels $D/fitness_browser_ec_reannotations_tierc.tsv \
  --threshold-labels $D/fitness_browser_ec_reannotations_tierc.tsv \
  --test-labels $D/fitness_browser_ec_reannotations_tierc.tsv \
  --matches-tsv $D/tierc_measured_panel.tsv \
  --bootstrap-resamples 2000
```

| Tier | Exit | Wall time | Peak RSS (GB) | arms.tsv | cases.tsv | cost.tsv | decisions.tsv | pathway_calls.tsv | orphan_steps.tsv |
|---|---|---|---|---|---|---|---|---|---|
| A (own CLEAN) | 0 | 4 min 55 s | 1.90 | 948 | 17,283 | 17,745 | 123,636 | 178,492 | 7,430 |
| B | 0 | 40 s | 0.62 | 316 | 1,575 | 1,743 | 52,101 | 12,385 | 2,478 |
| C | 0 | 2 min 33 s | 1.30 | 316 | 6,678 | 7,035 | 172,734 | 57,069 | 5,531 |

(Data rows, header excluded.) Each `out/` holds `arms.tsv, cases.tsv, cost.tsv, decisions.tsv, pathway_calls.tsv, orphan_steps.tsv, summary.json, report.md, staging_receipt.json` and `per-genome/`.

## 5. Things for the supervisor to check

1. Candidate row counts differ from sidecar EC-call counts because `clean_sidecar_to_candidates` drops score-0 rows and preliminary `n` EC numbers (section 1).
2. Two proteins were run through CLEAN from sanitized sequences (AAC23084.2 in GCA_000027305.1, Xcc-8004.1.1 in Xantho; see `clean/REPORT.md`); their candidate rows are included unchanged.
3. Tier B: GCA_000009985.1 contributes no masked case (no label rows); all tier B N2 runs use the measured-panel fallback.
