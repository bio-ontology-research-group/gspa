# P11: AllTheBacteria pan-phylum analysis from GSPA integrator output

Host: KAUST Ibex. Stage A–C inputs stay under `S=/ibex/scratch/projects/c2014/hohndor/allthebacteria`
(read only). Stage D outputs go to `R=/encrypted/e3008/allthebacteria/repro_2026-10-19/`; they lived
under `$S/repro_2026-10-19/` until 2026-10-20, when the c2014 project reached its 145 TB quota mid-run
and the tree was moved (rsync array 54046232, source removed after verification). Driver: `run.sh`. Output files are described in `OUTPUTS.md`. Checksums of
the inputs are in `manifest_rows.tsv`, in the column layout of `../../manifest/data_manifest.tsv`.

## Why this pipeline exists

The first ATB tables (`results/allthebacteria/insight*.tsv` and `atb_insights_summary.json`, made
by `aggregate_atb_insights.py`) summarise the per-genome `delta_events.tsv` files. Those files come
from a Python heuristic in `stageC_whole_genome_integrate.py` (about lines 540–655). The heuristic
applies fixed increments to the Stage B scores: operon partner +0.42·w, KEGG gap +0.24, complex or
module +0.20, essential +0.22, singleton ×0.35, and zero for taxon- or anatomy-forbidden terms.
**None of those numbers are GSPA output.** Stage C ran the real JVM integrator on only the first 20
genomes (sorted sample id) of each batch.

P11 recomputes every reported quantity from the JVM outputs alone.

## Stages: what each one really computes

| Stage | Script | What it computes | Caveats |
|---|---|---|---|
| A | `stageA_extract_and_linclust.py` | pyrodigal ORFs per genome, co-directional neighbour pairs (`operons.tsv`), exact-hash dedup, then MMseqs2 `easy-linclust --min-seq-id 0.9 -c 0.8` over the batch | Batch-local clustering. Every protein inherits the scores of its representative |
| B | `stageB_gspa_seq_inference.py` | Per representative: DIAMOND blastp `--sensitive -k 25` (default E ≤ 1e-3) against `ref/train_db.dmnd`, the CAFA (Kaggle) training Swiss-Prot subset. This gives `rep_max_pident.tsv` and DIAMOND-kNN GO transfer. Also InterPro-LR + Pfam2GO, a 16-kernel DeepGOCNN (CPU), three ESM-2 3B heads (kNN, BCE, C-HMCNN), and the GSPA-Seq stacker, closed upward over GO. Output: `rep_claims.jsonl` | **No ESM-2 3B embeddings were computed for any ATB representative.** Sequences that occur in `wg_all7`/`train` use cached embeddings. All other representatives get a bitscore-weighted interpolation of the embeddings of their DIAMOND hits (`train_esm2_3b.npz`), so zero-alignment orphans have no ESM channel at all (the live-extraction branch needs a GPU node and was not used). **No hmmsearch ran.** The InterPro/Pfam channel copies domain lists from Swiss-Prot hits at ≥ 35% identity, so it is also empty for orphans. In practice dark-matter proteins rely on the CNN channel alone. |
| C | `stageC_whole_genome_integrate.py` | (1) `claims.jsonl` per genome: top 35 GO terms + all ECs per protein, one source `esm2-go` with `n_supporting = 1`. (2) `gaps.jsonl`: KEGG pathways with between 40% and 100% (exclusive) of their GO terms covered at ≥ 0.45; the missing reactions are listed, at most 350. These are **GSPA-side KEGG-coverage gap claims, not gapseq**. (3) JVM `integrate --batch-manifest` on 20 genomes per batch. (4) heuristic `delta_events.tsv` for all genomes | (3) covers about 2% of the genomes. (4) is not GSPA output. Some early wave-2 batches ran an older Stage C: see "Legacy batches" |
| D (new) | `stageD_java_integrate_batch.py`, `stageD_java_integrate_all.sbatch` | The identical JVM call on **every** genome of a completed batch. `claims.jsonl` and `operons_chains.tsv` are regenerated with code copied verbatim from Stage C, in chunks of 96 genomes per JVM with a shared ontology load. Genomes Stage C already integrated with the same config are gzip-copied (`--adopt-existing`) | — |
| Aggregation (new) | `aggregate_atb_java.py` (`summarize` per batch, `combine`) | Reads only `noisy_or.tsv` (Stage 0: `max_iter 0`, no priors, no Stage 2) and `gspa_full.tsv` (`max_iter 3`, priors essentiality, coherence, consistency, anatomical, gap_filling, genomic_context, then `--enforce-stage2 --enforce-consistency-mode relax --lite`), together with `rep_max_pident.tsv`, the manifests and `gaps.jsonl` | Definitions are in the script docstring and in `OUTPUTS.md` |

The JVM call, identical in Stage C and Stage D:
`java -Xmx46g -jar gspa-unified.jar integrate --batch-manifest M --batch-threads T --go-owl ref/go.owl
--pathways ref/kegg_pathways.tsv --modules ref/kegg_modules.tsv --ec2go ref/ec2go.txt
--enable-priors essentiality,coherence,consistency,anatomical,gap_filling,genomic_context
--enforce-stage2 --enforce-consistency-mode relax --lite`. The jar is v2.1.0, sha256
`e415d9c6e13fbf96382ba7d798074b2ca4ec2991499de6714a5e4da2b19f72a7`, the same file as the local
`gspa-unified.jar`.

### Facts about the integrator that matter for the numbers

- `priors_fired` holds the per-prior boost **evaluated at the final refinement state**, not the
  cumulative shift. A prior that fired only in an earlier, damped iteration leaves a residual.
  For example, eukaryotic CC terms carry a −0.085 log-odds residual with an empty `priors_fired`.
- The Stage 2A relax step drops violating rows entirely. The aggregator classifies those rows with
  the jar's own constraint tables (`dump_forbidden_terms.groovy`).
- The taxon constraints act only when the manifest taxid is in the jar's bundled NCBI hierarchy,
  which holds 1,302 taxa. Only 43 of the 6,460 distinct manifest taxids are in it (`NCBITaxon:2`
  and a few model species). For every other genome the organism-level forbidden set is empty, and
  only genome-wide co-annotation SAT conflicts can fire.
- Every claim has `n_supporting = 1`, because Stage B stacks all channels into one source. Rescue
  paths that need ≥ 2 supporting claims, or a candidate with likelihood in [0.18, 0.50] matching a
  gap term, therefore rarely engage.

## Verification

- **Byte-compatibility (job 53986478, `debug`, batch_0001).** Stage D re-ran the 20 genomes that
  Stage C integrated. All 40 of 40 outputs (`noisy_or`, `gspa_full`) and all 40 of 40 regenerated
  inputs (`claims.jsonl`, `operons_chains.tsv`) are byte-identical to Stage C (`cmp`). The Stage D
  arrays depend on this job (`afterok`).

## Legacy batches (excluded from the main tables)

- `w2_0001`, `w2_0009`, `w2_0017` (6,854 genomes): an older Stage C ran the JVM on every genome but
  with empty claims (`n_claims = 0`), so the outputs are empty. Excluded.
- `w2_0013`, `w2_0018`, `w2_0019`, `w2_0020`, `w2_0021`, `w2_0024` (3,312 genomes): `max_iter 5`, and
  pairwise `operons.tsv` instead of chains. Summarised separately as `*_java_legacy`. Stage D
  re-runs all of these batches with the final config.

## Runtime of the JVM integrator (Task 1)

Source: `survey/java_runtime.tsv` (job 53986011). It holds every Stage C `java_batch_stats.tsv` row
of the final config, 4,740 genomes, joined with ORF counts. Stage C ran 16 single-threaded genome
jobs on 16 cores, so wall seconds per genome are close to CPU seconds.

| | median | p90 | mean |
|---|---|---|---|
| `gspa_full` row (s) | 223 | 417 | 248 |
| `noisy_or` row (s) | 8.5 | 18.5 | 10.1 |
| both (s / genome) | 233 | 430 | 258 |

Runtime grows with genome size: Spearman ρ = 0.65 against n_proteins (n = 4,740). A linear fit
gives ≈ 80 s + 0.045 s × n_proteins. By ORF stratum, the median total is 104 s (<1,500 ORFs),
168 s (1,500–2,499), 210 s (2,500–3,499), 285 s (3,500–4,499), 326 s (4,500–5,499) and 363 s
(≥5,500). The Stage D debug test measured 20 genomes in 529 s on 8 threads, ≈ 212 CPU-s per genome.

Cost of running the JVM on **all** genomes of the completed batches, at 258 CPU-s per genome:

- wave 1: 43 batches, 20,683 genomes, ≈ 1,480 CPU-h
- wave 2 completed at 2026-10-19 09:30: 204 batches, 449,653 genomes, ≈ 32,200 CPU-h
- total ≈ 33,700 CPU-h. One JVM shares the ontology load across 96 genomes, adding about 15 s per
  chunk.

## Slurm jobs (2026-10-19)

| Job | What | Notes |
|---|---|---|
| 53986011 | survey (`stageD_survey_batches.py`) | `survey/batch_inventory.tsv`, `java_runtime.tsv` |
| 53986494 | sha256 of the inputs | `checksums.tsv` → `manifest_rows.tsv` |
| 53986478 | byte-compat test (debug) | passed |
| 53986479 | Stage D, wave 1, `--array=1-43%20`, 8 CPU / 60 GB / 16 h | list `batch_lists/wave1.txt` |
| 53986492 | Stage D, wave 2, lines 1–60 of `batch_lists/wave2_done_2026-10-19a.txt`, `%12`, 16 CPU / 64 GB / 24 h | — |
| 54001443 | Stage D, wave 2, lines 61–197, `%12`, starts after 53986492 ends (`afterany`) | 7 or more wave-2 batches that finished after the list was written are not included |
| 53986492 + 54001443 outcome | 12 of 60 tasks of the first array completed; the c2014 quota filled at 04:19 on 2026-10-20, the running tasks died (Python exit 120, unflushable stdout) and every later task failed at start (Slurm signal 53, log file not writable) | 55 of 240 batches had partials; 12 wave-2 batches had partial Java output, which the batch script resumes |
| 54046232 | 78-task rsync of `repro_2026-10-19` from c2014 to `/encrypted/e3008/allthebacteria/` (177 GB, 124,770 files; dry-run comparison 54047473 showed no file differences; source removed by 54047617) | stale `w2_*/tmp` chunk files (39 GB) dropped from the copy; the batch script regenerates them |
| 54047084 | 694-task array gzipping every per-genome `claims.jsonl` of the Stage C complete wave-1 and wave-2 batches in place (pigz -6, integrity test); nothing downstream reads these files | frees about 17 TB on c2014 |
| 54047618 | `run.sh stageD-resume 20`: Stage D on the 185 wave-2 batches without a partial (`batch_lists/wave2_remaining_2026-10-20.txt`), 16 CPU / 64 GB / 1 day per task, 20 concurrent, outputs under the e3008 tree | 12 of them resume from partial Java output |
| 54100832 → 54232400 → 54232401 | `run.sh stageD-autoresume 54047618 60`: chain job after the resume array. The 4,000-genome batches run 20 to 29 h on 16 CPUs (96 genomes per 40 min) and hit the 24 h wall; the job resubmits every batch without a partial (the finished genomes are skipped) and chains itself until none is left (`batch_lists/autoresume_chain.tsv`) | 54047618 ended 2026-10-22 07:38 with 179 completed and 6 timed-out tasks (84–96 % of their genomes written); 54232400 finished the 6 batches by 11:36; 54232401 found no batch without a partial. 240 partials, no errors in the logs |
| 54492314 → … | `run.sh stageD-extend 50` (2026-10-24 20:30): Stage D on the 454 wave-2 batches whose Stage C finished after the 2026-10-19 list was drawn (`batch_lists/wave2_done_2026-10-24.txt`). The per-user QOS caps submitted jobs at 2,000 including array tasks (a first attempt, 54492141, failed with `QOSMaxSubmitJobPerUserLimit` with 1,644 jobs queued), so the autoresume job now submits only as many batches as fit and the chain feeds the rest after each array | running; about 90,000 CPU-h expected; combine with `run.sh combine $R/partials_stageD _stageD` when the chain reports completion |
| 54246242 | `run.sh combine $R/partials_stageD _stageD` → `results_stageD/`, fetched into `results/allthebacteria/` | 454,691 genomes from 240 partials, 2,979 duplicate copies of wave-1 representatives dropped |
| 53986460 / 53986461 | run (a) summarise, array version | slow to schedule, pending tasks cancelled |
| 53994239 | run (a) summarise, serial version (2 CPU) | completed. Partials in `partials_stageC/` (230 batches) and `partials_stageC_legacy/` (6) |

Batch lists live in `$R/batch_lists/`, one `wave<TAB>batch` per line: `wave1.txt`
(43), `wave2_done_2026-10-19a.txt` (197 batches complete when it was written), `existing_current.txt`
(230), `existing_legacy.txt` (6).

## Results

### Stage D: every genome of all completed batches (2026-10-22)

454,691 genomes in 240 batches (43 wave-1 + 197 wave-2), 37 phyla, 11,633 species, 1,779,514,422 ORFs;
each wave-1 species representative that reappears in a wave-2 batch is counted once (2,979 copies dropped).
Files: `*_java_stageD.tsv`, `insight3b_phylum_summary_java_stageD.tsv`, `atb_java_summary_stageD.json`.
Supplementary Tables S16 (survey) and S17 (integrator) are generated from these files.

| Regime | % of ORFs | rescued proteins | % of regime | % of regime unannotated by noisy-OR |
|---|---|---|---|---|
| dark < 30% id | 34.33 | 70,628 | 0.01 | 0.01 |
| (zero alignment, a subset of dark) | 29.51 | 14 | 0.00 | — |
| twilight 30–50% | 20.76 | 311,765 | 0.08 | 0.38 |
| moderate 50–70% | 7.76 | 12,435 | 0.01 | 0.08 |
| close ≥ 70% | 37.14 | 224,301 | 0.03 | 0.20 |

- 619,129 proteins (0.03% of ORFs) and 5,230,969 protein–GO pairs crossed 0.50; no protein lost all its
  ≥ 0.50 terms.
- Rescued proteins by prior (a protein can count under several): essentiality 478,313, coherence 129,053,
  stage2_relax 62,140, genomic_context 13,402, anatomical 2, gap_filling 0, consistency 0.
- Suppressed pairs, 47,506,880 in total: genome-wide SAT co-annotation conflict 26,363,207; taxon-forbidden
  18,985,742; other 2,021,859; removed by Stage 2A with no forbidden-table hit 136,072. 245,342 genomes
  (54.0%) have a taxid inside the jar's hierarchy.
- 1,897,198 KEGG-coverage gap claims; the gap_filling prior fired on 0 rows.
- Dark fraction by phylum: 28.9% in Pseudomonadota (261,841 genomes, 58% of the collection) to 65.2% in
  Bacteroidota; the species-representative snapshot (7,009 genomes, below) gives 49.5%.
- Spearman ρ against n_orfs (n = 454,691): rescued proteins per ORF −0.20; rescued pairs per ORF −0.81;
  dark fraction −0.51 overall (−0.34 below the median of 4,326 ORFs, +0.62 above).
- Median JVM time per genome inside the 96-genome batch runs: 216 s (Stage 0 + full run).
- The rates agree with run (a) and the snapshot (below), which are subsets of this run.

### Run (a): JVM outputs that Stage C already wrote (2026-10-19)

4,592 genomes in 230 batches (43 wave-1 + 187 wave-2), 16 phyla, 767 species, 17,986,002 ORFs.
Thirty wave-1 species representatives that also appear in a wave-2 batch were counted once. These
are the first 20 genomes by sorted sample id in each batch, so they are not a random sample. No
number below is extrapolated. Files: `*_java.tsv`, `atb_java_summary.json`.

| Regime | % of ORFs | rescued proteins | % of regime | % of regime unannotated by noisy-OR |
|---|---|---|---|---|
| dark < 30% id | 39.59 | 833 | 0.01 | 0.01 |
| (zero alignment, a subset of dark) | 33.83 | 0 | 0.00 | — |
| twilight 30–50% | 25.04 | 3,250 | 0.07 | 0.37 |
| moderate 50–70% | 9.54 | 177 | 0.01 | 0.10 |
| close ≥ 70% | 25.83 | 1,434 | 0.03 | 0.20 |

- In total, 5,694 proteins (0.03% of ORFs) and 54,034 protein–GO pairs crossed 0.50.
- Rescued proteins by prior (a protein can count under several): essentiality 4,231, coherence 1,256,
  genomic_context 248, stage2_relax 356, gap_filling 0, consistency 0, anatomical 0, unattributed 0.
- Rescued pairs by prior: essentiality 33,325, coherence 18,868, genomic_context 6,369,
  stage2_relax 581.
- Suppressed pairs, 568,068 in total: genome-wide SAT co-annotation conflict 350,588; taxon-forbidden
  184,327; removed by Stage 2A with no forbidden-table hit 1,945; other 31,208; anatomy-forbidden 0;
  complex singleton pruning 0. Only 1,966 of the 4,592 genomes have a taxid inside the jar's
  hierarchy.
- No protein lost all its ≥ 0.50 terms.
- The genomes carry 19,608 KEGG-coverage gap claims in total, but the gap_filling prior fired on 0
  rows.
- Complex coherence against singleton pruning (`insight2_…_java.tsv`): coherence crossed 0.50 on
  complex terms in only 5 phylum × term cells, mostly GO:0045254 (pyruvate dehydrogenase complex).
  Stage 2A never pruned an obligate-heteromer singleton.
- Genome size (`insight6_…_java.tsv`), Spearman ρ against n_orfs, n = 4,592:
  - rescued proteins per ORF: −0.275 (−0.17 below the median of 4,252 ORFs, −0.32 above). The rate
    falls monotonically: 0.068% (<1,500 ORFs) down to 0.023% (≥5,500).
  - rescued pairs per ORF: −0.845.
  - dark fraction: −0.33 overall, but −0.16 below the median and +0.50 above. Stratum means are
    63%, 47%, 50%, 35%, 33% and 45%, a U shape in the homology landscape, not in rescue.
  - gaps.jsonl per genome: +0.52.

The legacy configuration (3,312 genomes, `*_java_legacy`) gives the same picture: rescued
proteins are 0.09% of dark and 0.35% of twilight; 97,069 rescued pairs; 539,945 suppressed pairs.

### Stage D snapshot (whole batches, 2026-10-19 14:04)

`*_java_stageD_snapshot.tsv` and `atb_java_summary_stageD_snapshot.json` cover the 15 wave-1
batches finished by then: batch_0000_pilot, 0001–0005 and 0007–0015. That is 7,009 genomes (all
genomes of those batches), 29 phyla, 3,965 species and 26.1 M ORFs.

| Quantity | Value |
|---|---|
| dark | 49.5% of ORFs |
| zero alignment | 41.7% of ORFs |
| twilight | 34.4% of ORFs |
| rescue rate, dark | 0.02% |
| rescue rate, twilight | 0.14% (0.40% of twilight proteins unannotated by noisy-OR) |
| rescue rate, zero alignment | 0 |
| rescued proteins / pairs | 15,603 / 86,073 |
| suppressed pairs | 875,786 |
| ρ, rescued proteins per ORF vs n_orfs | −0.19 |
| ρ, dark fraction vs n_orfs | +0.10 |

Across whole batches the dark fraction is flat across strata (45–51%). The U shape in run (a) is
therefore likely an artefact of the first-20-by-id subsample.

### Against the heuristic tables

`heuristic_vs_java_same_genomes_java.tsv` compares the two on the 1,138 genomes present in both
the heuristic `insight3` table and run (a).

| Quantity | heuristic | GSPA JVM |
|---|---|---|
| rescue events / threshold-crossing pairs | 17,697,360 | 14,551 (≈ 1,200×) |
| operon (B6) / genomic_context | 14,305,677 | 1,582 |
| coherence (B2) / coherence | 3,383,895 | 5,055 |
| gap fill (B5) / gap_filling | 7,401 | 0 |
| essential (B1) / essentiality | 387 | 8,809 |
| suppressions | 37,860,429 | 167,142 |
| anatomical (B4) / anatomy | 37,197,462 | 0 |
| singleton (B2) / singleton prune | 20,881 | 0 |

The headline heuristic numbers, 42.9% of dark proteins and 77.0% of twilight proteins "rescued" over
72,533 genomes, have no counterpart in GSPA output. The JVM integrator rescues 0.01–0.02% of dark
proteins and 0.07–0.14% of twilight proteins. The homology landscape itself agrees: dark is 40–50%
of ORFs, against 41% in the heuristic table.

## Re-running

```bash
bash run.sh deploy
bash run.sh test                         # must print RESULT ... differ=0
bash run.sh stageD-wave1
bash run.sh stageD-wave2 61 197          # the rest of the completed wave-2 list
bash run.sh summarize-existing           # run (a)
bash run.sh combine $R/partials_stageC ""            ; bash run.sh fetch ""
bash run.sh combine $R/partials_stageC_legacy _legacy; bash run.sh fetch _legacy
# run (b), once Stage D has finished (each task writes its partial to partials_stageD/):
bash run.sh combine $R/partials_stageD _stageD       ; bash run.sh fetch _stageD
```

To add wave-2 batches that finished after 2026-10-19, append them to a new list and call
`stageD-wave2 FROM TO LIST`.
