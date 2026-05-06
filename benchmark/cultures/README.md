# Culture dark-matter prediction pipeline

Prospective functional-annotation pipeline for user cultures: given a
culture's assembly FASTA, produce ranked candidate assignments for
each reaction that gapsmith cannot fill, using genomic context
(anchor-density around pathway-neighbour catalysts) rather than
sequence similarity to the target.

See [`predictions/README.md`](predictions/README.md) for the actual
output format, dark-matter definition, and validation-ready shortlist.

## Scripts

| file | purpose |
|---|---|
| `manifest.tsv` | `tag<TAB>fasta_path` — one row per culture to process |
| `annotate_culture.sh` | SLURM array script: prodigal → DIAMOND + Pfam → claims → gapsmith → integrate |
| `predict_dark_matter.py` | Takes annotations for one culture, emits per-gap top-K candidates ranked by density (with anchor + gapsmith-assigned exclusion) |
| `embed_cultures.sh` | SLURM GPU array: ESM-2 t30 embeddings for culture proteins (for novelty scoring) |
| `score_novelty.py` | Cross-references predictions against culture ESM2 embeddings vs panel EC centroids (legacy absolute thresholds) |
| `score_novelty_v2.py` | Extended scoring: max cos over all ECs, bins predictions by ESM2 signal strength |

## Running a new culture

```bash
# 1. Add to manifest
echo -e "NEW-TAG\t/path/to/assembly.fa" >> manifest.tsv

# 2. Annotate (prodigal → DIAMOND + Pfam → gapsmith → integrate)
sbatch --array=<row-idx> annotate_culture.sh manifest.tsv

# 3. Dark-matter predictor
python3 predict_dark_matter.py \
    --tag NEW-TAG \
    --reactions-tbl .../NEW-TAG-all-Reactions.tbl \
    --layout .../NEW-TAG_layout.tsv \
    --integrated .../NEW-TAG_integrated.tsv \
    --reactions-tsv $GAPSMITH_DATA/seed_reactions.tsv \
    --diffusion-tsv $GAPSMITH_DATA/diffusion_mets.tsv \
    --ec-aliases-tsv $GAPSMITH_DATA/seed_Enzyme_Class_Reactions_Aliases_unique.tsv \
    --ec2go $REFERENCE/ec2go.txt \
    --tau 0.3 --top-k 5 --min-anchors 3 \
    --out NEW-TAG_dark_matter_predictions.tsv

# 4. (optional) ESM2 novelty scoring — needs GPU
sbatch --array=<row-idx> embed_cultures.sh
python3 score_novelty_v2.py \
    --predictions NEW-TAG_dark_matter_predictions.tsv \
    --plm-dir .../plm --tag NEW-TAG \
    --centroids-dir $BENCH/plm_centroids \
    --top-k 5 --out NEW-TAG_dark_matter_scored.tsv
```

All paths assume the `bench_gtdb30` panel has been assembled (see
`benchmark/panel/` and `benchmark/ml/`).
