# Paper reproducibility package

Everything the manuscript (`paper/main.tex`, `paper/supplement.tex`) reports is
regenerated from this folder. Nothing in the manuscript is typed by hand: tables
are generated LaTeX files under `paper/tables/`, figures are produced by scripts
under `figures/`, and every number quoted in prose is checked against a result
file by `checks/check_numbers.py`.

Layout

| Folder | Content |
| --- | --- |
| `manifest/data_manifest.tsv` | Every raw input: logical name, host and path (ibex / orix / unimatrix01 / local), size, sha256, producing job or script, date |
| `pipelines/<id>_<name>/` | One folder per benchmark: `run.sh` (or sbatch), the Python/Groovy entry points it calls, and `OUTPUTS.md` listing the `results/*.tsv` it writes |
| `tables/make_tables.py` | Reads `results/*.tsv`, writes `paper/tables/tab_*.tex` (one file per manuscript table) |
| `figures/` | `generate_figures.py` and friends; writes `paper/figures/*.pdf` |
| `checks/check_numbers.py` | Extracts every numeric value in the manuscript prose and verifies it against a results file or a documented derivation (`checks/number_ledger.tsv`) |
| `env/` | Environment specs (pixi/conda/requirements) for local, orix, ibex, unimatrix01 runs |

Large inputs excluded from git (Stage 2B ModelSEED bundle; whole-genome CAFA caches, embeddings, retrained-baseline predictions, GSPA-Seq checkpoint) are published as checksummed archives at https://bio2vec.net/data/gspa/ ; extract them at the repository root.

Regeneration order: `make results` (runs or re-aggregates each pipeline),
`make tables`, `make figures`, `make check`, then build the paper.

Pipelines

| id | Manuscript element | Compute host | Status |
| --- | --- | --- | --- |
| P01 | Table 1, Suppl. S4: whole-genome CAFA splits | orix (GPU baselines) + local CPU (GSPA) | done |
| P02 | Table 2, Fig. 2a-c: Stage 1/2A ablation and perturbation (both designs) | unimatrix01 | done |
| P03 | Fig. 2d, Suppl. S2: 13/21-genome panels | unimatrix01 | done |
| P04 | Table 3, Fig. 4a-b: leave-subunit-out | unimatrix01 | done |
| P05 | Table 4, Fig. 4c-d: eukaryote anatomy | unimatrix01 | done |
| P06 | Table 5, Fig. 3, Fig. 6c: leave-reaction-out | local | done |
| P07 | Fig. 6a, Suppl. S12: 108-genome stoichiometric audit | local | done |
| P08 | Table 6, Fig. 6b: Price-149 | local (+ibex raw predictor outputs) | done |
| P09 | Table 7, Fig. 6d: BacDive | local | done |
| P10 | Stage 3 benchmarks | unimatrix01 (vLLM) | done |
| P11 | Section 2.8, Suppl. S16-S17: AllTheBacteria | ibex | done for 7,009 species representatives; Stage D over all batches still running |
| P12 | Suppl. S8-S9: 1,000-genome and 88-proteome runs | unimatrix01 / ibex | done (six external tools); GSPA rows of S8 pending the ibex array |
