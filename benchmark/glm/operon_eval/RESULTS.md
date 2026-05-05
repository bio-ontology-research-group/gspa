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
each caller predicts "co-operonic" or "not". The candidate set is
restricted to pairs where both genes appear in *some* ODB4 / RegulonDB
known operon — this avoids penalising a caller for predicting
co-operonic on a pair the database simply hasn't surveyed.

Locus-tag bridging:

- `ecoli`: ODB4 uses `b0001`, my GFF `locus_tag` matches directly.
- `bsubtilis`: ODB4 uses `BSU34960` (no underscore), my GFF
  `old_locus_tag=BSU34960` matches; `locus_tag=BSU_34960` is the
  modern form.
- `hpylori`: ODB4 uses `HP1072`, my GFF `old_locus_tag=HP0001%2CHP_0001`
  bridges (URL-encoded comma).
- `paeruginosa`: both use `PA####` directly.

## Results — ODB4 multi-genome

```
genome        caller          TP    FP    FN    TN   prec    rec     F1
------------------------------------------------------------------------
ecoli         heuristic     1714   113    20    35  0.938  0.988  0.963
ecoli         gLM           1051    66   683    82  0.941  0.606  0.737
ecoli         gLM2           614    28  1120   120  0.956  0.354  0.517

bsubtilis     heuristic      640    35     8     5  0.948  0.988  0.967
bsubtilis     gLM            459     6   189    34  0.987  0.708  0.825
bsubtilis     gLM2           514    26   134    14  0.952  0.793  0.865

hpylori       heuristic      719    41     4     8  0.946  0.994  0.970
hpylori       gLM            490    14   233    35  0.972  0.678  0.799
hpylori       gLM2           562    22   161    27  0.962  0.777  0.860

paeruginosa   heuristic       68     1     0     0  0.986  1.000  0.993
paeruginosa   gLM             54     0    14     1  1.000  0.794  0.885
paeruginosa   gLM2            49     0    19     1  1.000  0.721  0.838
```

### Macro mean across 4 genomes

| caller | <P> | <R> | <F1> |
|---|---:|---:|---:|
| **heuristic** | 0.954 | **0.993** | **0.973** |
| gLM           | 0.975 | 0.697     | 0.812 |
| gLM2          | 0.968 | 0.661     | 0.770 |

Both foundation models are slightly **more precise** than the
heuristic (lower FP rate) but pay an enormous recall cost (~30% fewer
gold pairs found). The heuristic's near-perfect recall (99.3%) reflects
a basic property of bacterial genome organisation: co-operonic genes
are almost universally same-strand and ≤ 300 bp apart, so the trivial
adjacency rule captures essentially all of them.

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

1. **The phase-1 F-max NO-GO verdict is now even stronger.** The FMs
   aren't producing better operons in the first place — so even a
   downstream test more sensitive to operon quality wouldn't have
   helped them.

2. **Both FMs trade recall for precision.** Macro F1: heuristic 0.973
   vs gLM 0.812 vs gLM2 0.770. The FMs are on the *correct* end of
   the precision axis but their recall (~70% / ~66%) is the limiting
   factor in F1.

3. **gLM has a home-field advantage and still loses.** Its shipped
   logreg was trained on the very ground truth used for evaluation
   (gLM-shipped E. coli annot). It still scores below the heuristic.
   gLM2 wasn't trained on this data and scores lower still.

4. **The "redo with gLM2" verdict is consistent across two metrics
   and four genomes.** The line stays closed.

## What would change the picture

- **Drop the same-contig same-strand constraint from the heuristic.**
  Most of the heuristic's recall comes from this trivial rule. If you
  forced the comparison on a noisier candidate set (e.g. anti-strand
  pairs, longer intergenic distances), the FMs might pull ahead — but
  that's not how operons are defined biologically.

- **Score under-300-bp same-strand pairs only.** Restrict the candidate
  set further to make the FMs' precision-recall tradeoff matter more.
  Likely doesn't move the verdict.

- **Train a proper gLM2 operon predictor on E. coli ground truth** with
  a 190-dim attention-contact feature vector (gLM-style). The
  cosine-only logreg used here is the weakest possible operon
  classifier — a fair calibration could move gLM2's recall closer to
  gLM's, but unlikely past the heuristic.

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
