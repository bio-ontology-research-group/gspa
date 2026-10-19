# P02 — Stage 1/2A ablation grid and cross-clade perturbation (both designs)

Manuscript elements: Table 2, Fig. 2a–c (ablation); the nine-genome
cross-clade perturbation experiment (reported under both injection designs).
Compute host: unimatrix01 (SLURM, partition `debug`; never the login node).

## Contents

| File | Role |
| --- | --- |
| `run.sh` | sbatch wrapper: `submit` (three perturbation runs), `perturb` (one run), `collect`, `ablation` |
| `run_perturbation.py` | Perturbation experiment, `--design {v1,v2}`, `--jar`, `--runner {auto,batch,perjob}` |
| `ablation/01_run_ablation_grid.py` | Ablation grid A0–A9 over nine genomes (copied from `benchmark/unified/pillar1_2/`) |
| `ablation/02_run_cafa_eval.py` | CAFA micro / protein-centric F-max per ablation output |
| `ablation/03_run_gaef_eval.py` | Genome-aware (GAEF) metrics per ablation output |
| `ablation/05_compile_results.py` | Compiles `pillar1_2_ablation_*.tsv` and `pillar1_2_perturbation_summary.tsv` |
| `ablation/ablation_idfix_check.py` | Re-integrates the ablation grid on ID-normalised claims and compares with the published outputs (result: identical) |
| `collect_results.py` | Builds `results/pillar1_2_perturbation_{designs,summary}.tsv` from the per-run summaries |

The ablation scripts are byte-identical (md5) to the ones on unimatrix01 in
`/data/hohndor/gspa/proteomes/bench9_unified/` that produced the published
ablation (job script `run_pillar1_2_slurm.sh`, 2026-10-17). The ablation was
**not** re-run for this package; `run.sh ablation` documents how to.

## Perturbation designs

Both designs inject the same 100 GO terms into the same 100 proteins per genome
(one term per protein):

* Terms: 5 mandatory anatomical terms (GO:0007507, GO:0009908, GO:0001755,
  GO:0030282, GO:0048437) + 45 sampled from `go-anatomical-constraints.tsv`
  rows with an UBERON/PO/CL anatomy, and 50 sampled from `go-taxon-constraints.tsv`
  `only_in_taxon` rows for eukaryotic clades (NCBITaxon 2759, 33208, 33090,
  4751, 7711, 40674, 89593, 9606). Sampling: `random.Random(42)` shuffles of
  the sorted unique candidate lists (anatomical first, then taxon, same RNG).
* Proteins: `random.Random(md5(tag) mod 2**32).sample(sorted unique protein_ids, 100)`,
  with protein_ids taken after claims preprocessing. For the six genomes without raw
  FASTA IDs the sample is the same as before the fix.
* Scores: drawn from the same per-genome RNG immediately after the protein
  sample, one draw per term.

| | v1 (commit 326c112) | v2 (commit 3310b97, current) |
| --- | --- | --- |
| raw_score | U(0.55, 0.85) | U(0.62, 0.90) |
| primary source | `deepgo` | `diamond` (even term index) / `deepgo` (odd) |
| corroborating claim | none | 60% of pairs (index mod 5 in {0,1,2}): `pfam` if primary is diamond else `diamond`, score max(0.52, raw − 0.12) |
| claims per genome | 100 | 160 |
| configs in original | A1, A6, A8, A9 | A1, A6, A7, A8, A9 |

This package runs A1, A6, A7, A8, A9 for both designs. Note that `run.sh` excludes node007, which killed every batch job submitted to it on 2026-10-19 (signal 53). A survivor is an
injected (protein, term) pair with posterior ≥ 0.5 in the integrated output.

## Claims preprocessing (since 2026-10-19)

The merged claims of hpylori, mgenitalium and mjannaschii carry raw FASTA IDs
(`sp|ACC|NAME`) on all pfam claims, so before the fix 47 / 42 / 55 of their 100
injected pairs landed on phantom protein IDs that the survivor count never saw.
`run_perturbation.py` now preprocesses all nine claims files before protein
sampling (`--claims-ids normalised`, the default):

1. normalise `protein_id` matching `^(sp|tr)|ACC|NAME` to `ACC`;
2. drop exact duplicate claim records (identical JSON after step 1, keys sorted;
   the first occurrence is kept).

Use `--claims-ids raw` to reproduce the earlier runs. The ablation grid is not
affected; see `OUTPUTS.md`.

## Jars

| Jar (unimatrix01) | Date | sha256 | Notes |
| --- | --- | --- | --- |
| `/data/hohndor/gspa/bin/gspa-unified.jar` | 2026-10-17 12:48 | `e415d9c6e13fbf96382ba7d798074b2ca4ec2991499de6714a5e4da2b19f72a7` | v2.1.0, alpha_cons = alpha_anat_penalty = 0.68; supports `--batch-manifest` ("current jar") |
| `/data/hohndor/gspa/proteomes/bench9_unified/gspa-unified-bench.jar` | 2026-10-14 10:25 | `0553a2e445fae782bd05f1a8ba5f3c010aa84c5665c82a91570d89127c96c0eb` | alpha_cons = 3.0, alpha_anat_penalty = 3.5; no `--batch-manifest`, run one JVM per job ("old jar") |
| `/data/hohndor/gspa/proteomes/bench9_unified/.pillar1_2_private.jar` | 2026-10-17 11:09 | `0ca96abfc8d4f3419f1bbc87b8dc3e07c4f3de8f2ff1d4d5bc13b5a93348177d` | named in the 326c112 script; not used here |

## Inputs (unimatrix01)

| Input | Path |
| --- | --- |
| Merged predictor claims | `/data/hohndor/gspa/proteomes/bench9/with_iprscan/{tag}_claims_merged.jsonl` |
| gapseq gaps | `/data/hohndor/gspa/proteomes/bench9/gapseq/{tag}_gaps.jsonl` (modified 2026-10-17 09:25) |
| Operons | `/data/hohndor/gspa/proteomes/operons/{tag}_operons.tsv` |
| Truth (exp / all) | `/data/hohndor/gspa/proteomes/truth_dual/{tag}_truth_{exp,all}.tsv` |
| GO / KEGG references | `/data/hohndor/gspa/reference/{go.owl,go_aspect_map.tsv,kegg_pathways.tsv,kegg_modules.tsv,ec2go.txt}` |
| Constraint TSVs (term selection) | `/data/hohndor/gspa/proteomes/bench9_unified/go-{anatomical,taxon}-constraints.tsv` |
| Evaluator | `/data/hohndor/gspa/benchmark-py/benchmark_pgap_v2.py` (`--n-bootstrap 50`) |

`{tag}` ∈ bsubtilis, ecoli, ecolo157, hpylori, mgenitalium, mjannaschii, mtb,
paeruginosa, synechocystis. Sizes and sha256 are in
`paper-reproducibility/manifest/data_manifest.tsv` (`used_by_pipeline = P02`).

## Running

```bash
rsync -a paper-reproducibility/pipelines/P02_ablation_perturbation/ \
  unimatrix01:/data/hohndor/gspa/proteomes/bench9_unified/repro_2026-10-19/P02_ablation_perturbation/
ssh unimatrix01 'bash /data/hohndor/gspa/proteomes/bench9_unified/repro_2026-10-19/P02_ablation_perturbation/run.sh submit_idfix'
# optional: historical raw-ID runs
ssh unimatrix01 'bash /data/hohndor/gspa/proteomes/bench9_unified/repro_2026-10-19/P02_ablation_perturbation/run.sh submit'
```

See `OUTPUTS.md` for the result files and how they are copied into `results/`.
