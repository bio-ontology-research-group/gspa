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
| P01 | Table 1, Suppl. S4: whole-genome CAFA splits | orix (GPU baselines) + local CPU (GSPA) | being regenerated (calibration fix) |
| P02 | Table 2, Fig. 2a-c: Stage 1/2A ablation and perturbation (both designs) | unimatrix01 | pending |
| P03 | Fig. 2d, Suppl. S2: 13/21-genome panels | unimatrix01 | pending (raw outputs being located) |
| P04 | Table 3, Fig. 4a-b: leave-subunit-out | unimatrix01 | pending |
| P05 | Table 4, Fig. 4c-d: eukaryote anatomy | unimatrix01 | pending |
| P06 | Table 5, Fig. 3, Fig. 6c: leave-reaction-out | local | statistics being added |
| P07 | Fig. 6a, Suppl. S12: 108-genome stoichiometric audit | local | feasibility being computed |
| P08 | Table 6, Fig. 6b: Price-149 | local (+ibex raw predictor outputs) | pending |
| P09 | Table 7, Fig. 6d: BacDive | local | cohort sizes being corrected |
| P10 | Stage 3 benchmarks | unimatrix01 (vLLM) | pending regeneration |
| P11 | Section 2.8, Suppl. S7: AllTheBacteria | ibex | pending |
| P12 | Suppl. S8-S9: 1,000-genome and 88-proteome runs | unimatrix01 | pending |
