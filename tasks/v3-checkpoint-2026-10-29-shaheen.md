# Checkpoint 2026-10-29 (Shaheen phase) — resume here after context compaction

Written 2026-10-29 15:30 Asia/Riyadh by the supervisor (Fable, session on ws, repo `/mnt/data2/leechuck/gspa` =
`~/Public/software/gspa`). Live log: `tasks/v3-supervision-plan-2026-10-28.md` (append-only, timestamps from `date`).
Rules unchanged: frozen unimatrix jobs 10297/10298/10410/10432 untouched; no commits; negative results reported;
sub-agents via `agy-auto` (quota exhausted until ~18:30) or the Claude Agent tool; pi retired. Shaheen access: one
ControlMaster on ws kept alive by `~/tmp/shaheen-keepalive.sh` (log `~/tmp/shaheen-keepalive.log`); if it drops, ask
Robert for one OTP and run `~/.claude-max/skills/shaheen/scripts/shaheen-connect <code>`.

## Shaheen layout (k10018)
- Mirror of unimatrix bundles: `/project/k10018/gspa-v3/<bundle>` (read-only on compute, bound at
  `/data/hohndor/gspa-v3` inside containers). Writable outputs: `/scratch/project/k10018/gspa-v3/out/<bundle>`;
  per-job work: `/scratch/project/k10018/gspa-v3/work/<jobid>` (removed on success). **Inode quota 1,024,000 files
  per user on /scratch: extract packages once per node job; check `kuq` before fan-out.**
- Code slice: `/project/k10018/code/gspa-v3` (re-sync with the rsync in the log entry of 14:3x after local edits).
- Wrappers: `benchmark/v3/shaheen/{rolefence-arm,arm-array,arm-pack,source-context-resume,wp5_recovery,wp4-fenced-reference,wp4-gapsmith-search}.sbatch`
  (copies in `/project/k10018/gspa-v3/shaheen/`). arm-pack env: BUNDLE, ARMS, ARM_COLS, NAME_COLS (use ranges like
  `3-5`, never commas, in `--export`), MANIFEST_COL, COMMAND_COL, PARALLEL, EXTRA_EXCLUDES, TASKS.
- Image `/project/k10018/sw/images/debian-trixie` (+libexpat); envs `sw/envs/{gspa-v3,gspa-tools,diamond-2.0.15,diamond-2.1.11}`;
  JDK `sw/jdk/current`; gapsmith binary `/project/k10018/gspa-v3/tools/gapsmith-linux` (sha 8bb3cddd…).

## Jobs in flight (check `squeue -u hohndor`)
- 15691149 `gspa-resume` (24 remaining 10406 replay arms, 12 at a time, `source-context-resume.sbatch` MODE=arm-pack,
  out `/scratch/project/k10018/gspa-v3/out/source-context-resume-v2/{jobs,results}`). When all arms exit 0: submit
  `MODE=final` (same wrapper, `--export=ALL,MODE=final`) → `results/replay/verification.json`, analysis, final receipt;
  then mirror `results/` (without training) to ws and run the evaluation summariser / doc update (WP1 row 1/4/5).
- 15690710 `gspa-wp5-recovery` (stages reference → cobra on the existing bundle; out
  `/scratch/project/k10018/gspa-v3/out/wp5-recovery-v1/{bundle,runs,oracle,jobs/15690710}`). Then on unimatrix:
  qwen + transcripts stages (`benchmark/v3/wp5_recovery.sbatch`, needs the frozen bundle copied there); then locally
  `evaluate_wp5_recovery.py join → calibrate → evaluate --policy-sha256 → report`.
- WP4 searches: driver `~/tmp/gspa-wp4-full-v1/driver-search.sh` (log `driver-search.log`) writes manifests
  (`~/tmp/gspa-wp4-search-v1`, staged to `/project/k10018/gspa-v3/wp4-search-v1`) and submits 8 jobs
  `gspa-wp4-search-<set>-<sk>` (out `/scratch/project/k10018/gspa-v3/out/wp4-search-v1/<set>/<sk>/observations`).

## WP4 chain after the searches (design §10 steps 3–5; CLIs from the Task R/Q reports)
1. Mirror observations to ws: `rsync -a shaheen:/scratch/project/k10018/gspa-v3/out/wp4-search-v1/ ~/tmp/gspa-wp4-search-v1/`
   (the roots are then `~/tmp/gspa-wp4-search-v1/<set>/<sk>`).
2. `prepare_wp4_batches.py gaps --ec2go <ec2go> --out ~/tmp/gspa-wp4-full-v1/batches/gaps --search-root cohort=… --search-root peers=… --search-root removal=… --search-root fragment=… --peers-dir ~/tmp/gspa-wp4-full-v1/peers --perturbations <each removal/<g>> --perturbations <each fragment/<g>> [--cohort-replay ~/.cache/gspa-v3-data/cohort-pathway-replay-v2]`
   (one `--search-root NAME=ROOT` per set × superkingdom as Task S names them, e.g. `peers-archaea=`; ec2go path as in
   the pilot gaps receipt `~/tmp/gspa-wp4-pilot-v1/batches/gaps/receipt.json`).
3. `prepare_wp4_batches.py context --cohort-inputs ~/.cache/gspa-v3-data/cohort-context-modes-v2-inputs --cohort-root <fenced cohort root> --search-root peers=… --sequence ~/tmp/gspa-wp4-full-v1/sequence --out …/batches/context` once for cohort+peers, and once per removal/fragment output with `--search-root removal=… --metadata PERT/removal/<g>/cohort_metadata.json`.
4. `prepare_wp4_batches.py plan --blocks all --peers-tsv benchmark/v3/results/peer-selection-v1/peers.tsv --plans-10406 ~/tmp/gspa-role-genomes-v1/10406-plans --peers-dir … --sequence … --layouts ~/.cache/gspa-v3-data/training-cohort/layout --models-root /data/hohndor/gspa-v3/source-context-full-v1/jobs/10406/checkpoint/training/diamond_homology --context … --gaps … --perturbations PERT/removal/<g>=CONTEXT_REMOVAL … --perturbations PERT/shuffle --perturbations PERT/duplicate/<g> --cohort-claims-dir ~/tmp/gspa-wp4-pilot-v1/cohort-claims-10406 --path-map /home/leechuck/tmp/gspa-wp4-pilot-v1/cohort-claims-10406=/data/hohndor/gspa-v3/source-context-full-v1/jobs/10406/checkpoint/training/claims/diamond_homology --path-map /home/leechuck/tmp/gspa-wp4-full-v1=/data/hohndor/gspa-v3/wp4-full-v1 --path-map /home/leechuck/tmp/gspa-wp4-search-v1=/data/hohndor/gspa-v3/wp4-search-v1 --out ~/tmp/gspa-wp4-full-v1/batches/plans`
   (630 rows / 420 runs expected; D2: learned arms may be filtered with `--kinds none` first if WP2e models are not
   used; the 10406 models are the frozen ones per §6).
5. `audit_wp4_leakage.py --receipts <all receipt dirs> --fenced-report <bacteria report> --fenced-verify … --targets all_proteins.faa --out leakage.json` — gate.
6. Stage `batches/` to `/project/k10018/gspa-v3/wp4-full-v1/batches` (+ observations roots under `wp4-search-v1`),
   then arm waves with `arm-pack.sbatch`: `BUNDLE=wp4-full-v1 ARMS=batches/plans/arms.tsv ARM_COLS=3 NAME_COLS=3-5
   MANIFEST_COL=11 COMMAND_COL=12 PARALLEL=12 TASKS="<row indices>"` (waves per §10 step 4: R0 + parity first).
   Note: arms.tsv row = (block,target,regime,ctx,kind,variant,…); rows sharing one run_dir are one integration.
7. Mirror run outputs; `analyze_wp4_recurrence.py seal --arms <per-block arms.tsv> --runs … --out …` (per block;
   `--claims-dir` for cohort and sulfolobaceae), `score --labels ~/.cache/gspa-v3-data/training-cohort/labels/labels.tsv --seal … --out …`, `report`; archive `benchmark/v3/results/wp4-recurrence-v1/`; update row 6.

## Other open items
- Role-truth: archived `evaluations-rolefence-arm-v1` + doc section done; GuaB1 GMPR / cutinase cases pending.
- WP3 Price: unimatrix 10432 after 10410 (~2026-10-30); then `evaluate_price_reference` report.
- Release acceptance (row-by-row evidence) at the end; AllTheBacteria/EQ production per `tasks/v3-shaheen-orix-plan-2026-10-29.md`.

## Update 2026-10-29 16:40 (after the search launch)
- WP4 searches: jobs 15691729–15691736 (regenerated manifests with `/lustre2/project/...` code path; first hand
  submission 15691712–19 failed preflight, kept under `out/wp4-search-v1/failed-15691712-19/`). At 16:30 seven roots
  complete (cohort 10+8, fragment 60+48, peers 56+49, removal/bacteria 150; zero nonzero exits), removal/archaea 292/432.
  Observations ≈ 55 M / 9 files per genome, ≈ 45 G in total: they stay on Shaheen (ws disks at 94 %/96 %).
- WP4 chain now runs on Shaheen: `benchmark/v3/shaheen/wp4-batches.sbatch` (STAGES tests gaps context plan audit;
  container with the mirror at /data, scratch `out/wp4-search-v1` bound over `/data/hohndor/gspa-v3/wp4-search-v1`,
  writable `out/wp4-full-v1-out` at `/data/hohndor/gspa-v3/wp4-full-v1-out`, code slice at the pilot repo path,
  `--path-map wp4-full-v1-out=wp4-full-v1`). Inputs staged to the mirror and verified (stage-v4.src.sha256, 1,497 files):
  cohort-pathway-search-v1, cohort-pathway-search-archaea-v2, cohort-context-modes-v2-inputs, training-cohort,
  cohort-pathway-replay-v2; `ec2go.txt` in the code slice. After the job: copy `out/wp4-full-v1-out/batches` to
  `/project/k10018/gspa-v3/wp4-full-v1/batches` and `out/wp4-search-v1/<set>/<sk>/observations` into the mirror roots,
  then `arm-pack.sbatch BUNDLE=wp4-full-v1 ARMS=batches/plans/arms.tsv ARM_COLS=3 NAME_COLS=3-5 MANIFEST_COL=11 COMMAND_COL=12`.
- Blocker for `context`: Task S roots are per (set, superkingdom) but the fenced context requires metadata population ==
  union of observed genomes. Task T (Claude Agent, Opus; spec `~/tmp/gspa-v3-agy/task-T/spec.md`) adds
  `context --restrict-to-metadata` (+ `audit_subset`, receipt `restriction`). The wrapper already passes the flag.
  Until Task T lands run `STAGES="tests gaps"` only. agy Opus quota exhausted until ≈ 18:45.
- Resume: 15691762/15691795 failed (wrapper bugs: host `--logs` path; masked exit; unprobed `--resume-changed-source`
  declaration). Fixed wrapper deployed; arm-pack 15691804 running. Then `MODE=final`.
- WP5 15690710 still running (stages reference → cobra).
