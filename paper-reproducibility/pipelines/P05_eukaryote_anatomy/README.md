# P05: eukaryote anatomical / developmental constraint stress test (Table 4, Fig. 4c-d)

Entry point: `benchmark/unified/pillar3/run_pillar3_eukaryotic_anatomy.py` (stdlib Python 3)
Driver: `run.sh` in this folder
Host: unimatrix01, partition `debug`. node007 has no `/data` mount, so the driver excludes it.
Remote layout: `/data/hohndor/gspa/proteomes/pillar3_unified/repro_2026-10-19/P05_eukaryote_anatomy/`
(the scripts), with every output under `work/`.
Integrator: `/data/hohndor/gspa/bin/gspa-unified.jar` (v2.1.0, sha256
`e415d9c6e13fbf96382ba7d798074b2ca4ec2991499de6714a5e4da2b19f72a7`). Input checksums are in
`manifest_rows.tsv`.

## Reproduce

```bash
bash run.sh stage                                         # local: copy the scripts to unimatrix01
ssh unimatrix01 'bash .../P05_eukaryote_anatomy/run.sh submit'   # default: -c 16 --mem=90G, 6 JVMs in parallel
# when the debug nodes are busy, the paper run (job 9250, node003, 80 min) used:
#   SBATCH_EXTRA="-c 8 --mem=15G" WORKERS=1 XMX=11g bash run.sh submit
bash run.sh fetch                                         # local: copy the 3 TSVs into results/
```

The run takes 90 JVM runs (6 species x 5 configurations, on the real perturbed claims, the real
unperturbed claims and the synthetic claims), then the evaluation. Every output has a stamp file
that records the claims sha256 and the exact command. An output is reused only when its stamp
matches. The earlier script reused any existing output, and that is how its table ended up
mixing stale integrator outputs (see below).

## Design (`--input-source real`, the default)

- **Claims.** The real blind predictor claims are
  `multikingdom_batch/blind_claims/{species}_claims.jsonl`. They come from DIAMOND blastp against
  `reference_loo9` (e-value 1e-5, top 5 hits; hits from the same species or genus are excluded;
  raw_score = pident/100) and from Pfam->GO via pfam2go (raw_score 0.78).
- **Cohort.** The cohort is every protein with at least one real claim. All six species have fewer
  than 1,500 such proteins, so all of them are used:

  | Species | Proteins | Real claims |
  | --- | --- | --- |
  | S. cerevisiae | 964 | 15,375 (14,132 DIAMOND + 1,243 Pfam) |
  | C. elegans | 1,079 | 32,448 (30,757 + 1,691) |
  | D. melanogaster | 1,081 | 31,498 (29,688 + 1,810) |
  | P. patens | 147 | 1,706 (1,443 + 263) |
  | A. thaliana | 1,080 | 20,022 (18,263 + 1,759) |
  | M. musculus | 1,195 | 40,729 (39,016 + 1,713) |

- **Injection** (seed 42). One cross-clade developmental BP claim (DIAMOND, raw_score
  U(0.82, 0.95)) is added to 300 randomly chosen cohort proteins. P. patens gets 44, which is 30%
  of its cohort; the rule applies when a cohort has fewer than 600 proteins. The term pool is the
  same as in the synthetic design: developmental BP terms that are linked to UBERON/PO/CL in the
  GO-anatomy bridge and occur in the other five species' Swiss-Prot truth, plus the fixed mammal,
  plant and arthropod lists. Unlike the synthetic design, the injection only adds claims. The
  protein's own real developmental claims are kept. A draw that duplicates an existing claim of that
  protein is redrawn; this happened once.

  How many injected terms are forbidden for the species under the jar tables:

  | Species | Forbidden of injected | Already true |
  | --- | --- | --- |
  | S. cerevisiae | 291 / 300 | 1 |
  | C. elegans | 210 / 300 | 1 |
  | D. melanogaster | 214 / 300 | 5 |
  | P. patens | 40 / 44 | 0 |
  | A. thaliana | 282 / 300 | 0 |
  | M. musculus | 125 / 300 | 1 |

- **Attenuation** (seed 42). Each real claim pair is checked for being true (in the Swiss-Prot truth
  closure), allowed for the species, and developmental (bridge-linked BP). Eligible pairs have a
  real max raw_score above 0.62, so that attenuating them actually lowers the score. 20% of the
  eligible pairs are lowered to min(score, U(0.55, 0.62)). Real non-mammal proteomes carry very few
  such claims, so the developmental-rescue test is only informative for mouse, and partly for fly:

  | Species | True developmental pairs | Eligible | Attenuated |
  | --- | --- | --- | --- |
  | S. cerevisiae | 3 | 1 | 1 |
  | C. elegans | 84 | 22 | 3 |
  | D. melanogaster | 217 | 83 | 14 |
  | P. patens | 5 | 0 | 0 |
  | A. thaliana | 64 | 18 | 3 |
  | M. musculus | 1,211 | 1,141 | 224 |

- **Configurations.** E1..E5 use the same `gspa integrate` flags as before. They are listed in the
  script docstring. E2 includes `--enforce-stage2 --enforce-consistency-mode=remove`, as in the
  earlier script.
- **Cohorts run.**
  - The perturbed cohort ("injected").
  - The unperturbed real claims ("endogenous"), run through E1..E5 as well. Endogenous E4 is a real
    run, not the earlier formula.
- **Truth.** Unmodified Swiss-Prot GO annotations (`multikingdom_batch/blind_truth/{species}_truth_all.tsv`),
  closed over is_a + part_of.
- **Violations.** An output row counts as a violation when its posterior is at least 0.50 and its
  term is forbidden for the species. Forbidden terms are the direct GO taxon constraints plus the
  GO->anatomy->taxon bridge, closed downward over GO. They are computed from the constraint tables
  extracted from the same jar the JVM runs.
- **Metrics.**
  - Micro F-max (`cafa_fmax`, the earlier definition) and CAFA protein-centric F-max, each on all
    cohort proteins and on the injected subset.
  - Mean IC of raw terms with posterior >= 0.50.
  - For the endogenous cohort, protein-centric F-max over all 1,200 Swiss-Prot panel proteins (155
    for P. patens). This is the definition used in `multikingdom_per_genome_eval.tsv`. The
    evaluator reproduces that file's noisy-OR values to 4 decimals for all six species.

## `--input-source synthetic` (the earlier design)

This mode regenerates the earlier claims, built from each protein's Swiss-Prot truth with
U(0.80, 0.95). The six regenerated claims files are byte-identical to the files in
`pillar3_unified/` (checked in `run_info.json`).

- **Constraint tables.** It uses the earlier copy in `pillar3_unified/resources`. That copy lacks
  one bridge row that the jar has: GO:0048513 "animal organ development" -> UBERON:0000062
  "organ" (only_in Metazoa).
- **Endogenous columns.** These use the earlier procedure, flagged in `endog_e4_method`:
  - F-max parsed from `multikingdom_per_genome_eval.tsv`.
  - E4 F-max = max(E1, E5 - 0.0008).
- **What reproduces exactly** with the current jar:
  - all B1 and E1 numbers;
  - every S. cerevisiae number;
  - all endogenous columns.
- **What does not reproduce.** For the other five species, the earlier `pillar3_unified/{species}_e{2..5}.tsv`
  outputs differ from what the current jar produces on byte-identical claims. Examples:
  C. elegans E3 differs in 1,532 rows; mouse E2 in 55,369 rows. Neither `gspa-unified.jar` nor
  `gspa-pillar3.jar` reproduces them, and rerunning from the earlier working directory changes
  nothing. They are stale outputs from an earlier build that the earlier script reused because they
  already existed.

  Effect on the synthetic rows:

  | Quantity | Earlier table | Current jar |
  | --- | --- | --- |
  | E2 violations | 11 / 185 / 420 / 0 / 0 / 0 | 11 / 392 / 866 / 5 / 25 / 13 |
  | E3 violations | 16 / 208 / 564 / 1 / 0 / 0 | 16 / 215 / 572 / 1 / 12 / 13 |
  | Injected-subset E4 / E5 F-max, mean | 0.9629 / 0.9629 | 0.9657 / 0.9440 |

  E4 and E5 violations are 0 under both.
