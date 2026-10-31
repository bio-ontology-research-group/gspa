# Observed pathway search and adapter repair

An actual gapsmith 0.1.0 / DIAMOND 2.0.15 search on the acquired training
proteome GCA_000012525.1 exposed an unusable pathway adapter. The old production
parser returns **zero annotations** from the complete 13,147-row reaction table.
It expected a synthetic six-column format and an `active` status; the tool emits
28 named columns with alignment statuses and target identifiers in `stitle`.

The corrected adapter parses the complete observed table into **3,469 candidate
claims for 802 proteins**: 2,207 EC candidates and 1,262 reaction/subunit
candidates. This is evidence preservation, not an accuracy result.

## External execution and provenance

Slurm job **10111** completed with exit 0 on node004. The observer copied and
checked 41,274 reference files on node-local scratch, ran the frozen binary on
the exact acquired 2,661-protein FASTA, then rehashed the source and staged
references. The resulting tables contain 1,923 pathways, 111 predicted present,
and the following reaction-row statuses:

| Status | Rows |
|---|---:|
| good_blast | 2,239 |
| bad_blast | 5,135 |
| no_blast | 4,467 |
| no_seq_data | 1,002 |
| spontaneous | 304 |

No held-out functional labels were supplied. Tool/reference training overlap
and biological pathway accuracy have not been established. Input, binary and
reference hashes, complete tables, resources and terminal Slurm state are in
[the evidence archive](../../benchmark/v3/results/gapsmith-observed-v2/README.md).
The initial staging failure, job 10110, is preserved separately: it omitted the
reference-directory sentinel. The repeat copies the original sentinel table and
preflights table loading before copying the sequence collection.

## Adapter behavior

The default search target is now the complete proteome. Gapsmith's protein
aligner no longer receives nucleotide contigs. Pathway configuration now reaches
the adapter, including tool, executable, reference directory, taxonomy, target
and threads. `annotate --threads` overrides its thread budget. Gapseq taxonomy
uses `-t`; the old command incorrectly passed it to the bit-score option `-b`.
The installed gapseq shell source verifies that command syntax, but a full
gapseq search has not been certified by this experiment.

Named-column parsing validates table shape, statuses, booleans, identifiers,
statistics and hit coordinates. Unknown target proteins fail instead of being
silently discarded. Missing protein sequences produce partial coverage. A
valid header-only table is successful with zero results; malformed output is
failed. Nucleotide gapseq hits cannot become protein claims without coordinate
attribution; this missing capability is reported as failure when hits occur.

Strong and weak alignment hits remain candidates. Their raw score is percent
identity divided by 100, with bit score, E-value, coverage, reference accession
and every original supporting row preserved. This value is **not a calibrated
probability**. The evidence is typed as sequence similarity, so it does not
become independent metabolic support in the combiner.

Duplicate hits/pathway entries for a protein/function produce one candidate,
retaining all supporting rows in `hitsJson`. Complex hits produce `CUSTOM`
reaction-subunit candidates, rather than asserting the complete EC activity of
each subunit. Hits without a usable EC produce reaction candidates. These IDs
include both pathway and reaction because custom pathways reuse IDs such as
`rxn1`. Partial ECs remain partial; reference-local identifiers such as
`1.5.1.M4` remain in raw records and reaction candidates, not the EC namespace.

The pipeline writes `<genome>_<tool>_observations.json`, preserving every reaction
row, including absent-hit and unknown-reference rows, the table hash, search
target and outcome. Reports clear between genomes and failed repeated calls. Reusing an output
directory removes an obsolete report when the current run has none.
Subunit metadata is retained, but a trained role-aware assignment layer is still
required; these records do not complete item 7.

## Validation through production routes

Nine input/configuration regressions and fourteen initial output regressions
fail before repair. The final change adds 34 JVM cases. The complete gate runs
**1,147 tests**, zero failures/errors and two optional skips; all four live
HMMER/AntiFam tests execute. All nine live Nextflow checks pass on the final JAR.
The reordered synthetic fixture also passes its 97 nearest adapter/parser
checks.

Seven actual CLI runs replay the verified tables through the repaired engine.
The replay executable checks that submitted protein IDs/sequences exactly match
the observed proteome and verifies the forwarded thread/taxonomy options. It
does not rerun the external search. An independent table-based comparison checks
every candidate identity, score and complete supporting-row multiset. Annotation
with and without integration retains identical raw claims and full observations;
annotation and direct-integration prediction TSVs are byte-identical.

Both current learned models receive the same claims. All 3,469 candidates per
model are explicitly unsupported for the unseen `gapsmith` source; none receives
a learned probability or changes the heuristic prediction TSV. Injected empty
and malformed outputs yield `SUCCESS_ZERO` and `FAILED`, respectively.

The first repaired build passed its sampled fixtures but rejected local EC-like
identifiers in the complete real table. That failure, intermediate engine hash
and replay are retained. Added regressions now cover those identifiers and
partial ECs; the complete-table replay passes on the final build:
`bc9487caf1d2a4b1b787eff2eb40fbad70350f3dda9545aa87c719c99ef2b6f9`.

## Genuine candidate-gap preparation

The [reaction-level gap converter](v3-gap-table-contract.md) now also passes a
complete-table audit on this observation. It produces **50 records covering 46
eligible reactions**, versus 234 records before repair. The old records include
161 with a strong hit elsewhere for the reaction, 17 involving complexes,
17 spontaneous and seven with unknown reference coverage; these reason counts
overlap. Every new record belongs to a predicted pathway, has suitable negative
search evidence, and has an order-independent weak-hit flag. All eligible
reaction groups are represented. Forty-eight records have an EC-to-GO mapping.

These are search-supported hypotheses, not experimentally verified missing
enzymes. Automatic gap/context wiring, explicit pathway namespaces, subunit
inference, calibration, and matched held-out context ablations remain open.

## Reproduction

Use new output directories and disk-backed `TMPDIR` and Java scratch paths.
With the observed files and proteome at their recorded locations:

```bash
python benchmark/v3/check_gapsmith_route.py --jar "$JAR" \
  --fixture "$OBSERVATION" --proteome "$PROTEOME" --out "$REPLAY" \
  --models benchmark/v3/results/cohort-evidence-dedup-v1/models
python benchmark/v3/check_gapsmith_gaps.py --observation "$OBSERVATION" \
  --ec2go benchmark/whole_genome_cafa/data/reference/ec2go.txt \
  --before benchmark/v3/results/gap-table-contract-v1/before-parser.py \
  --after benchmark/parse_gapseq_gaps.py --out "$GAP_AUDIT"
```
