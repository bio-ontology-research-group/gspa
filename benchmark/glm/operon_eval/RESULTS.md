# Direct operon-prediction quality eval

The phase-1 F-max ablation scored *GO predictions* against GOA, not the
operons themselves. This evaluation closes that gap by scoring each
operon caller's predictions against two independent operon ground
truths, on multiple genomes.

## Ground truth sources

### 1. gLM-shipped E. coli annot (RegulonDB-lineage)

`gLM/repo/data/ecoli_operon_data/operon.annot` — 4,305 E. coli K-12
gene names with operon assignments (2,554 in named operons). Ships
with the gLM repo; the gLM paper trained its operon LogReg on this
exact file via 5-fold CV.

Used by `eval_operons.py`.

### 2. ODB4 (Operon Database v4)

<https://operondb.jp/> — multi-genome curated operon database.
9,480 known operons across many organisms. Per-genome coverage of our
benchmark set:

| Genome (taxid) | ODB4 known operons |
|---|---:|
| E. coli K-12 MG1655 (511145)        | 2,845 |
| B. subtilis 168 (224308)            |   658 |
| H. pylori 26695 (85962)             |   413 |
| P. aeruginosa PAO1 (208964)         |    33 |
| (mtb / mgenitalium / mjannaschii / ecolo157) | 0–4 — **excluded** |

Downloaded raw at `odb4_known.txt`. Used by `eval_odb4.py`.

DOOR3's URL (`intelligentoffloading.com/DOOR3/`) was unreachable as of
2026-05-05. RegulonDB is E. coli-only and was not added separately
since gLM's shipped annot covers the same data with locus-tag bridging
already done.

## Evaluation protocol

We score the **adjacent-pair "same operon" classification**: for genes
*a* and *b* sitting consecutively on the same contig and same strand,
each caller predicts "co-operonic" or "not".

### Methodology fix (2026-05-06)

The first version of this eval restricted the candidate set to "pairs
where BOTH genes are in some ODB4 known operon." That was a
sample-selection trap: ODB4 surveys known operons, so within that
subset 92–99% of adjacent same-strand pairs are co-operonic, and any
permissive classifier scored F1 ≈ 0.96–0.99 just by mirroring the
prior. The "0.97 F1 for the heuristic" result that produced was an
artifact, not signal — if a one-line heuristic could really score
F1 ≈ 0.97 there'd be no specialized operon-prediction literature.

The corrected eval uses **all** adjacent same-strand same-contig pairs
as candidates. This accepts that ODB4 has incomplete coverage and that
some "negatives" are actually unsurveyed positives, which hurts
permissive callers slightly more than conservative ones. The
trade-off is unavoidable without a fully-curated ground truth.

The genome set is also restricted: paeruginosa (1.7% positive rate
under all-cand mode → eval is pure noise) and bsubtilis (21% positive
rate, suggesting ~80% of "negatives" are unsurveyed positives) are
excluded. Only ecoli (57.5%) and hpylori (64.6%) have ODB4 coverage
dense enough to give a meaningful eval.

`eval_odb4.py` reports both candidate-set modes (`all` and
`restricted`) so the discrepancy is visible.

### Locus-tag bridging

- `ecoli`: ODB4 uses `b0001`, my GFF `locus_tag` matches directly.
- `bsubtilis`: ODB4 uses `BSU34960` (no underscore), my GFF
  `old_locus_tag=BSU34960` matches; `locus_tag=BSU_34960` is the
  modern form.
- `hpylori`: ODB4 uses `HP1072`, my GFF `old_locus_tag=HP0001%2CHP_0001`
  bridges (URL-encoded comma).
- `paeruginosa`: both use `PA####` directly.

## Results — ODB4, all-candidate set (corrected)

Per-genome (only ecoli + hpylori have usable ODB4 coverage):

| Genome | pos rate | caller | TP | FP | FN | TN | Prec | Rec | F1 |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|
| ecoli | 57.5% | **heuristic** | 1714 | 1024 |   20 |  257 | 0.626 | 0.988 | **0.767** |
| ecoli | | gLM | 1051 |  871 |  683 |  410 | 0.547 | 0.606 | 0.575 |
| ecoli | | gLM2 |  614 |  274 | 1120 | 1007 | 0.691 | 0.354 | 0.468 |
| hpylori | 64.6% | **heuristic** |  719 |  309 |    4 |   87 | 0.699 | 0.994 | **0.821** |
| hpylori | | gLM |  490 |  214 |  233 |  182 | 0.696 | 0.678 | 0.687 |
| hpylori | | gLM2 |  562 |  257 |  161 |  139 | 0.686 | 0.777 | 0.729 |

### Macro mean (ecoli + hpylori)

| caller | <P> | <R> | <F1> |
|---|---:|---:|---:|
| **heuristic** | 0.663 | **0.991** | **0.794** |
| gLM           | 0.621 | 0.642     | 0.631 |
| gLM2          | 0.689 | 0.566     | 0.598 |

These numbers are now in the range of published operon-prediction
tools — DOOR ~0.85, ProOpDB ~0.80, Operon-mapper 0.80–0.85, the gLM
paper's reported PR-AUC ~0.67. The heuristic's F1 ≈ 0.79 is
respectable but doesn't match specialized methods like DOOR; it
benefits from bacterial genome compaction (co-operonic genes really
are overwhelmingly same-strand and ≤ 300 bp apart) but the precision
is a real ceiling on what the trivial rule can reach.

### Excluded genomes

- **paeruginosa**: 33 ODB4 known operons → 1.7% gold positive rate
  among all adjacent same-strand pairs → all callers score F1 ≈ 0.04.
  ODB4 coverage too sparse; the eval is dominated by unsurveyed pairs
  treated as negatives.
- **bsubtilis**: 658 ODB4 known operons but 21% positive rate, meaning
  ~79% of "negatives" are likely operons ODB4 hasn't surveyed. All
  callers' F1 sits at 0.36 — driven by the FP rate of unsurveyed
  positives, not by detection skill. Numbers reported in the eval log
  for transparency but don't reflect operon-prediction quality.

### The biased-candidate-set numbers (for transparency)

The original (botched) eval mode is left in the script as
`mode='restricted'` and produces:

| caller | <P> | <R> | <F1> |
|---|---:|---:|---:|
| heuristic | 0.954 | 0.993 | **0.973** |
| gLM       | 0.975 | 0.697 | 0.812 |
| gLM2      | 0.968 | 0.661 | 0.770 |

These are inflated by ~0.18–0.30 across the board because the
candidate set was 92–99% positive — any classifier admitting most pairs
wins. The relative ordering is preserved (heuristic > gLM > gLM2) but
the absolute magnitudes don't reflect operon-detection skill.

## Results — gLM-shipped E. coli annot (corroborates above)

| Caller | TP | FP | FN | Precision | Recall | F1 |
|---|---:|---:|---:|---:|---:|---:|
| heuristic | 1710 | 1061 |   25 | 0.617 | 0.986 | **0.759** |
| gLM       | 1043 |  872 |  692 | 0.545 | 0.601 | 0.572 |
| gLM2      |  631 |  259 | 1104 | 0.709 | 0.364 | 0.481 |

The earlier eval against the gLM-shipped annot used a *different*
candidate set (every adjacent same-strand pair, not just ODB4-surveyed
ones), which is why the precision numbers are lower across the board
(more "fake negatives" inflate FP). The relative ordering is the same.

## Implications

1. **The phase-1 F-max NO-GO verdict is corroborated, but with more
   modest margins than the botched eval suggested.** On the trustable
   genome subset (ecoli + hpylori), heuristic F1 ≈ 0.79 vs gLM 0.63
   vs gLM2 0.60. The FMs are real operon callers — they're just not
   beating the trivial adjacency rule on these benchmarks.

2. **Both FMs trade recall for precision but the precision gain is
   small.** On ecoli, gLM2 has the highest precision (0.69) but the
   lowest recall (0.35). The heuristic at 0.63 P / 0.99 R wins on F1.

3. **gLM has a home-field advantage and still loses.** Its shipped
   logreg was trained on the very ground truth used for evaluation
   (gLM-shipped E. coli annot, RegulonDB-lineage). It still scores
   below the heuristic at F1 ≈ 0.58. gLM2's hand-tuned cosine
   threshold scores 0.47.

4. **The "redo with gLM2" verdict is consistent across two metrics
   and the trustable genomes.** The line stays closed.

## What would change the picture

- **Train a proper gLM2 operon predictor on E. coli ground truth** with
  a 190-style attention-contact feature vector (per-pair, like gLM
  does). The cosine-only sigmoid used here is the weakest possible
  operon classifier — a fair calibration could move gLM2's F1 closer
  to gLM's ~0.58, but unlikely past the heuristic's 0.79.

- **Get a denser multi-genome operon ground truth.** ODB4 coverage is
  the limiting factor on bsubtilis and paeruginosa. RegulonDB is
  E. coli-only; DOOR3 was unreachable; ProOpDB might fill the gap.

- **Score on the dark-matter slice.** The phase-1 F-max ablation tested
  on Swiss-Prot-rich genomes where the homology stack saturates.
  The operon-quality eval here is on the operon decision itself, but
  whether better-than-trivial operon decisions matter *downstream*
  is still an open question on the EQ-MAG-style sparse genomes.

## Reproducibility

```bash
# ODB4 known operons (538 KB)
curl -L https://operondb.jp/download/known_operon.download.txt \
     -o benchmark/glm/operon_eval/odb4_known.txt

# Multi-genome eval (uses local copies of operons + GFFs)
python3 benchmark/glm/operon_eval/eval_odb4.py

# gLM E. coli eval (uses local ecoli_operon_data/operon.annot)
python3 benchmark/glm/operon_eval/eval_operons.py
```
