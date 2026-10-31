# Plan: Shaheen III (CPU) and Orix (GPU) for GSPA completion, AllTheBacteria + Empty Quarter annotation, and metabolism discovery

Drafted 2026-10-29 10:58 Asia/Riyadh by the supervisor (Fable) at Robert's request. Status: **proposal, nothing
submitted**. The stop order of 10:51 stands; the frozen Price chain on unimatrix01 is untouched.

Budget given: 10,000,000 core-hours on Shaheen III (CPU); Orix for GPU-heavy work only (CPU-only or preparation
work under a GPU allocation is a compliance flag there).

## 0. Facts the plan rests on

Shaheen III CPU partition (docs.hpc.kaust.edu.sa, checked today): 4,608 nodes, 2 × AMD EPYC 9654 (192 cores,
384 threads), 384 GB DDR5 per node; **exclusive node allocation, 192 core-hours charged per node-hour whatever
is used**; `workq` 24 h walltime, `72hours` queue on request with justification; Lustre `/scratch` (32 PB) and
`/project`; `/home` is not mounted on compute nodes; compute nodes are on a private network (no internet;
login nodes have it); Singularity/Apptainer is supported. Login is OTP only: Robert opens one interactive
session, the agent reuses the ControlMaster socket (`~/.ssh/config`, `ControlPersist 12h`). 10 M core-hours =
52,083 node-hours = 2,170 node-days, e.g. 100 nodes for three weeks.

Orix: 64 H100 (80 GB) + 64 H200 (141 GB), 128 vCPU and 1.5 TiB per node, 3-day walltime, `pi-hohndor` reserved
half/full node, `freecycle` preemptible (GPU type must be named), no modules (Pixi/conda in `$HOME`), heavy
data under `/mnt/data/u/hohndor`.

Cost measurements we already have (ibex, Empty Quarter run, 17,055 genomes, `benchmark/emptyquarter_mags/`):
channels once per mmseqs linclust representative (90 % id / 90 % cov, 2.3× reduction): 18.6 k CPU-h for
~23 M representatives ≈ **2.9 CPU-s per representative** (InterProScan dominates); gapsmith model median
86 s wall on 4 CPUs (≈0.1 CPU-h per genome); Stage 0 + 1 + 2A + 2B median 174 s on 4 CPUs (≈0.19 CPU-h per
genome). AllTheBacteria v0.2: 2,440,377 assemblies; the earlier GSPA work covered an alignment survey of
454,691 genomes (1.78 B CDS) and ran the real integrator only on 20 genomes per batch, with one collapsed
evidence channel and interpolated embeddings (paper validation report, item 5). No AllTheBacteria genome has
yet been annotated with the production workflow (per-channel claims, gapsmith model, Stages 1–2B).

Genoa cores are faster than the ibex cores the figures come from; the plan uses the ibex figures as an upper
bound and a pilot to replace them.

## 1. Budget split (core-hours, Shaheen) and GPU-hours (Orix)

| Work package | Shaheen core-h (estimate) | Orix GPU-h | Notes |
|---|---:|---:|---|
| A. GSPA completion and all benchmark experiments | 1.0–1.5 M | ~300 | Moves the unimatrix queue; rerun of frozen comparisons under the accepted JAR |
| B1. AllTheBacteria production annotation | 2.5–4.5 M | ~800 | Dominated by InterProScan on representatives; pilot decides |
| B2. Empty Quarter rerun under v3.0.0 (+ neural channels) | 0.1 M | ~100 | Cheap; inputs already staged on ibex |
| C. Metabolism discovery (FBA / pFBA / FVA / knockouts / communities) | 0.3–0.8 M | 0 | LP work, HiGHS; scales with condition set |
| Reserve (pilot misestimates, reruns after bug fixes, Lustre/packing inefficiency) | ~2.5 M | | Not pre-allocated |

Everything is committed in the order A → B1 pilot → C pilot → B1 production → C production → B2, with gates
(section 6). If the B1 pilot shows the uncompressed channel bound (≈6 M core-h for 7.3 B CDS without
clustering) rather than the clustered one, the fallback in section 3.3 applies before any production wave.

## 2. Work package A: finish GSPA and the benchmark experiments on Shaheen

What moves: everything that was running or queued on unimatrix01 when the stop order came, plus the reruns
that acceptance needs. The frozen artefacts (10406 fits, Price collection) stay where they are and are read,
never regenerated in place.

1. **Resume replay of 10406** (`run_source_context_remote.py resume`, bundle `source-context-resume-v1`):
   on Shaheen the 720 production genome runs of the replay spread over one node (192 cores) instead of 12
   CPUs; expected < 24 h, i.e. inside `workq`. Inputs: the 79 G training restore, the sealed package. The
   NFS EIO class of failure disappears (Lustre), the seal/receipt logic is unchanged.
2. **WP2e role-fence retraining** (`run_rolefence_arm.py`, `rolefence-arm-v1`): 4 contexts × 18 genomes
   export + fits; one node, < 24 h.
3. **21-genome joint context arms** (ORGANISM logistic/gbdt, PATHWAY gbdt still open) and archiving of the
   evaluations into `benchmark/v3/results/role-truth-v1/`.
4. **WP4 cohort recurrence** (design `tasks/v3-wp4-recurrence-design-2026-10-29.md`): fenced gapsmith
   reference build, peer acquisition for every block, R0/R0-SC/R1/R2 × 3 contexts × 3 kinds, perturbations
   (fragmentation, near-duplicate, clade). Each arm is one `gspa integrate --batch-manifest` call; one node
   runs ~40 arms concurrently. Budget ≈ 0.3 M core-h because of the perturbation grid.
5. **WP5 masked recovery** (design `tasks/v3-wp5-masked-recovery-design-2026-10-29.md` §8): same shape.
6. **Price full matrix under the accepted JAR**: the frozen collection (132 cells, 1.53 B claims) was made
   with the original production JAR; acceptance row 1 needs the same matrix under the release candidate.
   Large cells need 192 GiB heap (fits a 384 GB node, one cell per node, ~14 h); medium cells pack 8 per node.
   ≈ 0.25 M core-h. Six external predictors' GPU inference (CLEAN, GraphEC, EnzBert, DPZ, …) for the 22
   genomes: Orix, a few GPU-hours.
7. **Nutrient phenotypes, community, absence, masked-assignment benchmarks** at their full catalogues
   (today they run bounded subsets because of unimatrix capacity): 6,864 conditions × solver arms, the
   1,320 coculture cases, Sco-GEM endpoints. LP only; ≈ 50 k core-h.
8. **Full JVM gate + Python suites** per candidate JAR on one node (replaces the 10402/10526 build-and-test
   jobs on node004).

Operating rules carried over unchanged: never alter the frozen jobs/artefacts, unconfirmed labels are never
negatives, no release claim without empirical evidence per acceptance row, report negative results. The
cluster-specific bits (sbatch headers, scratch layout, `/data` paths) become a thin `shaheen` profile next to
the existing `unimatrix` one; the harness scripts themselves do not change.

Orix part of A: ESM-2 3B embeddings and the neural sidecar heads (DeepGOPlus-ESM2, ProteInfer, CLEAN,
ESM-2 centroids) for the benchmark genomes and the reaction centroids; the learned-context features that
need embeddings. Container: Pixi env with torch cu12 built on the login node without `-march=native`.

## 3. Work package B1: annotate all of AllTheBacteria with the production workflow

Scope decision proposed: **all 2,440,377 assemblies of v0.2** (HQ subset 1.9 M as the fallback). Output is
the same per-genome product as the Empty Quarter run (Stage 0 table, Stage 1/2A/2B table, actions, RLGC
suggestions, gapsmith model, gap report, score), plus the v3 context stages once they are accepted. ATB
remains excluded from benchmark acceptance (goal text); it is the production application.

### 3.1 Staging (login/datamover nodes and ibex, no compute charge)

- Assemblies: ATB v0.2 batch tarballs (~1.6 TB xz) to `/scratch/<project>/atb/assemblies/` via the Shaheen
  datamovers (or from the ibex copy under `/ibex/scratch/projects/c2014/hohndor/allthebacteria` through
  DataWaha). Bakta annotations where ATB released them; pyrodigal 3.7.1 otherwise (EQ rule).
- Reference data, one read-only tree: Swiss-Prot experimental-GO DIAMOND db, Pfam-A, InterProScan 5.61 data,
  gapsmith references (bacterial and archaeal), ModelSEED universal network + media, taxon constraints, GO
  and EC2GO, GTDB-Tk-free taxonomy from ATB's own species calls (GTDB r220 names) mapped to GSPA's taxon
  table with the EQ mapper (`00_select.py` logic).
- One Apptainer image: Java 21 + the accepted GSPA JAR, pyrodigal, mmseqs2, DIAMOND 2.1.16, HMMER 3.3.2,
  InterProScan 5.61, gapsmith 0.1.0 (with the medium-CSV fix), the Python harness env (scipy/HiGHS, cobra,
  libsbml). Built on ws (x86-64-v3, no native flags), tested on one Shaheen node.

### 3.2 Pipeline (same design as EQ, re-shaped for exclusive 192-core nodes and Lustre)

| Step | What | Packing | Cost model |
|---|---|---|---|
| S0 catalogue + assembly dedup | species, quality, identical-proteome hashing; batches of 1,000 genomes by species | login/pp nodes | 0 |
| S1 proteins | Bakta import or pyrodigal | 192 genomes concurrently per node | 2.44 M × ~30 CPU-s ≈ 20 k core-h |
| S2 cluster | linclust 90/90 **per species shard** (not one pool of 7 B proteins); shards > 100 M proteins split by genus/ANI cluster; cross-shard second pass on shard representatives | one shard per node, 384 GB | ~30 k core-h |
| S3 channels | DIAMOND + Pfam + InterProScan on representatives, chunks of 50,000 | 24 chunks × 8 CPUs per node, 24 h jobs, auto-requeue by manifest | reps × 2.9 CPU-s: 500 M reps → 0.4 M; 1 B reps → 0.8 M |
| S4 scatter | claims to members, per batch | per batch | 10 k |
| S5 gapsmith | model per genome (`-t Archaea` where needed) | 48 × 4 CPUs per node | 2.44 M × 0.1 ≈ 0.25 M |
| S6 integrate | Stage 0 (`--enable-priors none`) and Stages 1+2A+2B via `integrate --batch-manifest`, 96 genomes per JVM, shared ontology load | 12 JVMs × 16 threads per node | 2.44 M × 0.19 ≈ 0.47 M, likely less with the shared load |
| S7 aggregate | per genome row, per species/phylum summaries, dark-fraction tables that replace the heuristic ones (S7/S17/S18) | pp nodes | small |

Expected total with 500 M–1 B representatives: 1.2–1.7 M core-h; with the ×1.3 packing/Lustre factor and
reruns, the 2.5–4.5 M range above.

Lustre rules: no per-genome directories at 2.4 M scale (that is 50 M inodes). One directory per 1,000-genome
batch; per-genome outputs inside batch-level tar or sharded `tsv.gz`; progress in a per-batch SQLite manifest
(DONE rows), not DONE files. Each node-job pulls work from the manifest, so a 24 h wall clock cut loses at
most one genome per slot; `--requeue` plus manifest resume.

### 3.3 Fallbacks if the pilot says channels are too expensive

In order: (i) InterProScan member databases reduced to Pfam, TIGRFAM, CDD, Hamap for non-representative
genomes while species representatives keep the full set; (ii) cluster at 80 % identity for the channel step
only (claims are entries and GO terms, not coordinates; the EQ Jaccard test 0.9985–1.0 is repeated at 80 %
before adopting); (iii) HQ subset only. Each fallback is disclosed in the output README.

### 3.4 Orix part of B1

GPU-only steps, run on the representative FASTA shipped from Shaheen (≈150 GB for 500 M sequences):
ESM-2 3B embeddings → GSPA-Seq / DeepGOPlus-ESM2 heads and ESM-2 centroid channel (this closes the "no ESM-2
forward pass was run" gap of the paper); CLEAN EC predictions; ProstT5 3Di strings for a FoldSeek structural
channel on the dark fraction (FoldSeek search itself is CPU and goes back to Shaheen). Throughput on H100
with bf16 batching is several hundred sequences per second for ESM-2 3B, so 500 M representatives ≈ 300–500
GPU-h per model; ~800 GPU-h all in, i.e. a few days on the `pi-hohndor` node. Only claims (protein, term,
score) come back, never raw embeddings (2.5 TB).

Stage 3 (LLM agent) is **not** part of the B1 production wave: acceptance row 11 has not shown an agent
advantage. If it does, Stage 3 runs on the ~12 k ATB species representatives and the EQ genomes only, with a
vLLM server on an Orix H200 node (Qwen3.8-27B for parity with the unimatrix endpoint, or a larger open model
if Robert prefers), ≈ 15 B tokens ≈ 1,000–1,500 GPU-h. That is a separate gate.

## 4. Work package B2: Empty Quarter under the released workflow

The 17,055 EQ genomes already have proteins, clusters, channel claims and gapsmith models on ibex. Once
v3.0.0 is accepted: copy the staged inputs to Shaheen (or run on ibex as before; 25 k core-h either way),
rerun S6 with the release JAR and the context stages (cohort recurrence within each sample's bins is exactly
the WP4 setting), add the Orix neural channels for the ~23 M representatives (~100 GPU-h), and refresh
`benchmark/emptyquarter_mags/results/`. Also run the never-executed `run.sh compare` (gapsmith vs gapseq).

## 5. Work package C: metabolism discovery under multiple conditions

Inputs: one gapsmith model per genome from B1/B2 (SBML-FBC), gap-filled per gapsmith's own medium
prediction, plus the GSPA Stage 2B flux-certified EC assignments. Solver: the existing stack
(`scripts/metabolic_experiments.py`, `service/metabolic_conditions.py`, the COBRApy checkers in
`benchmark/v3/check_*_cobra.py`), HiGHS LPs, every record conditional on model/medium/bounds, no automatic
promotion (the rules of `docs/scenarios/metabolic-experiments.md` stay).

Condition grid, tiered by cost:

| Tier | Conditions | Per genome | Population | Core-h |
|---|---|---|---|---|
| C1 growth screen | 17 gapsmith media (ALLmed, LB, LB+glc, TSB, M9 variants, MM glucose ± thiamine, MM anaerobic acetate / acetate+H2 / CO2+H2, gut, seawater, autotrophic, YCFA, FT) × O2 uptake 0 / 2 / 20 mmol gDW⁻¹ h⁻¹ → 51 conditions; FBA + pFBA | ~100 LPs, 20–60 s | all 2.46 M | 20–40 k |
| C2 source panels | Biolog-style single carbon (~190), nitrogen (~95), sulfur (~35), phosphorus (~60) sources on the minimal base; pFBA each | ~380 LPs, 1–3 min | all 2.46 M | 50–120 k |
| C3 pathway usage | pFBA flux per pathway/module under C1 conditions; predicted secretions (max/min secretion at a growth floor, classes as in the scenario doc); fermentation products under O2 = 0 | included in C1 solves | all | — |
| C4 essentiality + FVA | single-gene knockouts (all GPR genes) and FVA at 90 % growth under 6 representative media | 2,000–5,000 LPs, 5–15 min | species representatives (~12 k ATB) + all 17 k EQ | 10–20 k |
| C5 auxotrophy and supplementation | the 20 amino acids, vitamins and nucleobases as supplements on minimal media; matches the nutrient-dependence machinery already in the repo | ~60 LPs | all | 20 k |
| C6 community cross-feeding | pairwise and sample-level coupled models (`check_community_cobra.py` formulation, measured-uptake convention) for the 542 EQ samples (bins of one sample) and for ATB co-occurring species sets chosen by Robert | 10²–10⁴ LPs per set | EQ samples; ATB sets on request | 50–200 k |

Deliverables: per genome a `conditions.tsv` (condition, growth, pFBA objective, secretions, limiting
exchanges) and `sources.tsv`; per species and phylum the predicted metabolic phenotype matrix; a "discovery"
table of genomes whose predicted phenotype differs from their species majority, traced to the GSPA
assignments that create the difference (these are the candidates for experimental follow-up and for the
agent-guided experiments of requirement 12). Validation of C on the known-phenotype sets we already hold
(Baichman-Kass nutrient phenotypes, the Sco-GEM endpoints, BacDive/KEGG growth media for ATB type strains)
runs inside A, before C production.

## 6. Gates and timeline (calendar weeks from access)

| Week | Milestone | Gate |
|---|---|---|
| 0 | Robert: Shaheen project/account id, storage quota (ask 150 TB scratch, 50 TB project), OTP session; Orix `pi-hohndor` confirmed | — |
| 1 | Container built and tested on one node; **pilot 10,000 ATB genomes** (5 nodes, 24 h) and the EQ 1,000-genome replay; Orix embedding throughput measured | G0: measured core-h per genome × 2.44 M ≤ 4.5 M, else section 3.3 |
| 1–4 | WP A items 1–8 | G1: acceptance rows with evidence; v3.0.0 candidate JAR frozen |
| 2–3 | C1/C2 pilot on the 10 k pilot genomes; validation on known phenotypes | G2: solver time per genome and agreement on known phenotypes reported |
| 3–8 | B1 production waves (S1–S7), ~100 nodes sustained | per-wave receipts, inode and quota checks |
| 4–9 | B1 Orix channels on representatives; C1–C5 production | — |
| 9–10 | B2 EQ rerun; C6 communities; aggregation, dark-fraction and phenotype tables; paper sections | G3: release v3.0.0 only if every acceptance row holds |

Supervision model stays as now: Fable supervises and validates; pi-opencode (kimi-k3.1-pro) or agy
(Opus/Sonnet) implement scripts; every cluster job has a receipt; negative results are reported.

## 7. Decisions needed from Robert

1. Shaheen project/account to charge, and whether to request the `72hours` queue (useful for the large
   Price cells and WP4 grids; everything else fits 24 h).
2. ATB scope: all 2.44 M assemblies (proposed) or the HQ 1.9 M subset.
3. InterProScan: full member set for every representative (proposed, pilot permitting) or the reduced set
   for non-species-representatives from the start.
4. Stage 3 agent at scale: gated on acceptance row 11 (proposed), and which model to serve on Orix.
5. Condition grid for C (proposed tiers C1–C6); any specific media or environments to add (e.g. Empty
   Quarter soil chemistry as a custom medium).
6. Where the ~50–100 TB of ATB outputs live long-term (DataWaha?).
7. Whether to restart the stopped unimatrix work now or hold everything until Shaheen access (the frozen
   Price chain is still running on unimatrix01 either way unless you say otherwise).

Sources for the Shaheen figures: https://docs.hpc.kaust.edu.sa/policy/shaheen3.html ,
https://docs.hpc.kaust.edu.sa/systems/shaheen3/compute_node.html ,
https://docs.hpc.kaust.edu.sa/soft_env/job_schd/slurm/shaheen_cpu_jobs.html ,
https://docs.hpc.kaust.edu.sa/systems/shaheen3/login_node.html .
