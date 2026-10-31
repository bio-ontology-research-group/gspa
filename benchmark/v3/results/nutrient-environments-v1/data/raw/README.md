# Nutrient environments shape amino acid auxotrophy and cross-feeding

This repository contains the data and analysis workflow associated with the manuscript:

> **Nutrient environments shape amino acid auxotrophy and cross-feeding**  
> Amichai Baichman-Kass, Lianet Noda-García, and Jonathan Friedman

The study examines amino acid auxotrophy and pairwise cross-feeding in *Bacillus subtilis* and *Escherichia coli* across combinations of carbon and nitrogen sources. The repository can be used to inspect the source data, reproduce the processed tables and statistical analyses, rerun the random-forest model, and regenerate the manuscript figures.

## Contents

- [Analysis scope](#analysis-scope)
- [Repository structure](#repository-structure)
- [Quick start](#quick-start)
- [Dependencies](#dependencies)
- [Analysis workflow](#analysis-workflow)
- [Input-data key](#input-data-key)
- [Plate and sample-map key](#plate-and-sample-map-key)
- [Generated-table documentation](#generated-table-documentation)
- [Abbreviations](#abbreviations)
- [Reproducibility notes](#reproducibility-notes)

## Analysis scope

The primary analysis includes **40 nutrient environments**, generated from ten supplied nitrogen sources and four carbon-source conditions.

`NoN` measurements are no-nitrogen controls and are **excluded from `filtered_df.csv` and from all downstream analyses and models**. They may remain in the raw plate-reader files, `Sample Names.xlsx`, and the reconstructed pre-filtering table `big_norm_df.csv` so that the original plate layout and processing steps remain traceable.

## Repository structure

```text
.
├── README.md
├── Consolidated Analysis and Figures.ipynb
├── Data/
│   ├── Crossfeeding Experiments/
│   │   ├── Bs_data_Ex1.xlsx
│   │   ├── Bs_data_Ex2.xlsx
│   │   ├── Bs_data_Ex3.xlsx
│   │   ├── Ec_data_Ex1.xlsx
│   │   ├── Ec_data_Ex2.xlsx
│   │   └── Ec_data_Ex3.xlsx
│   ├── AA Requirments Experiments/
│   │   ├── AA_per_OD_df_Bsub_ex1.csv
│   │   ├── AA_per_OD_df_Bsub_ex2.csv
│   │   ├── AA_per_OD_df_Ecoli_ex1.csv
│   │   └── AA_per_OD_df_Ecoli_ex2.csv
│   ├── AA_synthesis_pathways.csv
│   ├── Sample Names.xlsx
│   ├── Koo et al Data.xlsx
│   ├── Tong et al Data.xlsx
│   ├── Wetmore et al Data/
│   │   ├── fit_logratios_good.tab
│   │   ├── fit_genes.tab
│   │   └── fit_quality.tab
│   ├── Nichols et al Data.xlsx
│   ├── Monk et al 2017.xlsx
│   ├── Orth et al 2011.xls
│   └── Ito et al 2005.xlsx
├── Generated Tables/
│   ├── AA_requirments_dataset.csv
│   ├── big_norm_df.csv
│   ├── filtered_df.csv
│   ├── Supplemented_AA_data.csv
│   ├── Koo_auxs_df_bs.csv
│   ├── Tong_aux_df_ec.csv
│   ├── Wetmore_aux_df_ec.csv
│   └── Random Forest Model/
│       ├── perm_test_df.csv
│       ├── perm_cv_df.csv
│       ├── block_cv_fast.csv
│       ├── block_cv_folds_fast.csv
│       └── block_cv_perms.csv
│       └── best_model.joblib

└── images/
```

## Quick start

1. Download or clone the repository.
2. Open a terminal in the repository root.
3. Create and activate a Python environment.
4. Install the dependencies listed below.
5. Start Jupyter and run `Consolidated Analysis and Figures.ipynb` from top to bottom.

For example:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install jupyterlab numpy pandas openpyxl xlrd scipy matplotlib seaborn scikit-learn statsmodels statannotations joblib cld4py requests
jupyter lab
```

On Windows, activate the environment with `.venv\Scripts\activate`.

The notebook metadata records Python 3.12.7. Exact package versions from the final analysis environment should also be recorded in a `requirements.txt` or `environment.yml` file before publication.

## Dependencies

| Dependency | Installation name | Purpose in this repository |
|---|---|---|
| Python | `python` | Notebook runtime; the supplied notebook records Python 3.12.7 |
| JupyterLab or Jupyter Notebook | `jupyterlab` | Opens and runs the analysis notebook |
| NumPy | `numpy` | Numerical operations, arrays, random permutations, and error propagation |
| pandas | `pandas` | Imports, reshapes, merges, filters, and exports tabular data |
| openpyxl | `openpyxl` | Allows pandas to read `.xlsx` workbooks |
| xlrd | `xlrd` | Allows pandas to read `Orth et al 2011.xls` |
| SciPy | `scipy` | Regression, correlations, statistical tests, and geometric means |
| Matplotlib | `matplotlib` | Base figure generation and export |
| seaborn | `seaborn` | Heatmaps, scatterplots, boxplots, and other statistical plots |
| scikit-learn | `scikit-learn` | Random forest, preprocessing, grouped data splits, model evaluation, and permutation importance |
| statsmodels | `statsmodels` | Multiple-testing correction |
| statannotations | `statannotations` | Statistical annotations on plotted comparisons |
| joblib | `joblib` | Parallel permutation analyses and optional local model serialization |
| cld4py | `cld4py` | Compact-letter displays for multiple comparisons |
| Requests | `requests` | Retrieves the BiGG gene-name mapping used by the Monk analysis; replacing this with a local mapping would remove the live web dependency |

Python standard-library modules such as `pathlib`, `itertools`, `statistics`, and `warnings` do not require installation.

## Analysis workflow

`Consolidated Analysis and Figures.ipynb` contains four main stages:

1. **Import source data** — load the current-study plate-reader data, sample map, amino-acid requirement summaries, pathway annotations, and published knockout or fitness datasets.
2. **Generate analysis tables** — blank-correct endpoint OD measurements, reconstruct sample identities, calculate means and SEMs, normalize coculture growth, identify conditionally non-obligate auxotrophs, remove `NoN` controls, and export the processed tables.
3. **Fit and evaluate the random-forest model** — split data by species–environment groups, select model parameters on the training data, evaluate the held-out groups, and calculate feature and feature-block permutation importance.
4. **Generate figures** — create the main and supplementary figure panels in `images/`.

The cached `best_model.joblib` file is not required as a permanent repository file. The model should be reproducible from the processed data and recorded settings. If a local cache is used to save time, the notebook should train the model when the cache is unavailable.

## Input-data key

### Current-study experimental data

| File or folder | Contents | Data used by the notebook |
|---|---|---|
| `Data/Crossfeeding Experiments/Bs_data_Ex1.xlsx`–`Bs_data_Ex3.xlsx` | Endpoint OD600 plate-reader exports for *B. subtilis*; one worksheet per strain pair or control | The 96-well numeric matrix from each worksheet. Experiment 1 is read after 32 header rows from columns A:M; the other workbooks are read after 23 header rows from columns B:N. |
| `Data/Crossfeeding Experiments/Ec_data_Ex1.xlsx`–`Ec_data_Ex3.xlsx` | Endpoint OD600 plate-reader exports for *E. coli*; one worksheet per strain pair or control | The 96-well numeric matrix from each worksheet, read after 23 header rows from columns B:N. `Ec_data_Ex2.xlsx` has no `F+S` worksheet; available measurements are retained. |
| `Data/Sample Names.xlsx` | Plate map connecting each well to its nitrogen and carbon environment | The complete A1:H12 grid. Environment labels are divided into `N` and `C` at the `+` separator. |
| `Data/AA Requirments Experiments/AA_per_OD_df_Bsub_ex1.csv` and `...ex2.csv` | Processed amino-acid titration regression summaries for *B. subtilis* | `Aux`, `Environment`, `Num_of_points`, `R2`, `Slope`, `C`, and `N`; `Unnamed: 0` is treated as a saved index and is not analyzed. |
| `Data/AA Requirments Experiments/AA_per_OD_df_Ecoli_ex1.csv` | Processed amino-acid titration regression summaries for *E. coli* | `Aux`, `Environment`, `Num_of_points`, `R2`, `Slope`, `C`, and `N`; `Unnamed: 0` is treated as a saved index. |
| `Data/AA Requirments Experiments/AA_per_OD_df_Ecoli_ex2.csv` | Second processed *E. coli* amino-acid titration experiment | `Aux`, `R2`, `Slope`, and `Num_of_points`; `Time` is removed and the environment is assigned as `Gluc_Glcn`. |
| `Data/AA_synthesis_pathways.csv` | Curated amino-acid biosynthesis pathway annotations | `Gene` is used to define the knockout genes retained from the external datasets. Entries containing `+` or `/` are split into individual gene names. `Unnamed: 19` is an empty export column and is not analyzed. |

The amino-acid requirement CSVs are processed inputs rather than raw titration measurements. `aa_per_od` is calculated as `1 / Slope` in the notebook.

### Published datasets used for auxotrophy comparisons

| File | Source section | Columns used and processing |
|---|---|---|
| `Data/Koo et al Data.xlsx` | Sheet `B`, beginning after row 24 | `locus_tag`, `gene`, all columns beginning with `nRF_ErmR_`, and all columns beginning with `nRF_KanR_`. The final analysis retains 13 specified environments and records the resistance cassette as `Erm` or `Kan`. |
| `Data/Tong et al Data.xlsx` | Sheet `S1A - Endpoint Biomass` | `B-numbers`, `Gene`, and all 30 environmental endpoint-biomass columns. Values are normalized environment-by-environment to the `kdsC` reference before genes with relative growth below 0.1 are selected. |
| `Data/Wetmore et al Data/fit_logratios_good.tab` | Quality-filtered gene-fitness matrix | `locusId`, `sysName`, `desc`, and experiment columns matched to approved records in `fit_quality.tab`. Replicate values sharing a short environment name are averaged. |
| `Data/Wetmore et al Data/fit_genes.tab` | Gene annotation table | `locusId` and `name`; `name` is renamed `gene` and joined to the fitness matrix. |
| `Data/Wetmore et al Data/fit_quality.tab` | Experiment metadata and quality table | `u`, `name`, `short`, `num`, and `t0set`. Only `u == True` experiments are retained, and time-zero records are excluded. |

### Published datasets used for supplementary heatmaps

| File | Source section | Columns used and processing |
|---|---|---|
| `Data/Nichols et al Data.xlsx` | Sheet `C_N_sources` | `Gene` and the eight condition-score columns: acetate, N-acetylglucosamine, glucose, glycerol, maltose, succinate, glucosamine, and NH4Cl. |
| `Data/Monk et al 2017.xlsx` | Sheet `Table 7` | `c_source`, `gene`, and `upper asymptote`. The `gene` field contains BiGG/b-number identifiers, which are mapped to gene names before analysis. |
| `Data/Orth et al 2011.xls` | Sheet `Table 1` | Excel columns B and I:P: gene identifier and eight endpoint-OD measurements. Duplicate plate columns are averaged. |
| `Data/Ito et al 2005.xlsx` | Sheet `Supplemental Table 2` | `Gene name deleted` and all binary phenotype columns. `Cluster` and `Appearance in Supplemental Figure` are not used. |

`Data/Aida et al 2024.xlsx` is not used by the supplied notebook. If it remains unused after the final notebook review, it should be removed from the repository to avoid suggesting that it contributes to the published analysis.

## Plate and sample-map key

Each cross-feeding workbook contains one worksheet for each strain pair or control. Worksheet names encode the pair being measured.

| Worksheet pattern | Meaning |
|---|---|
| `WT` | Prototroph control; converted to `WT+WT` during processing |
| `M+M`, `F+F`, ..., `L+L` | Unsupplemented auxotroph monocultures |
| `M+Met`, `F+Phe`, ..., `L+Leu` | Amino-acid-supplemented monocultures |
| `M+F`, `M+R`, ..., `C+L` | Pairwise auxotroph cocultures |
| `Blank` | Blank-control worksheet |

The six auxotroph labels are `M`, `F`, `R`, `C`, `S`, and `L`. See [Abbreviations](#abbreviations) for their full names.

### Well layout

| Wells | Carbon condition | Replicate structure |
|---|---|---|
| Rows A–B | Glucose (`Gluc`) | Duplicate wells for each nitrogen condition |
| Rows C–D | Glycerol (`Gly`) | Duplicate wells for each nitrogen condition |
| Rows E–F | Fumarate (`Fum`) | Duplicate wells for each nitrogen condition |
| Rows G–H | No added carbon (`NoC`) | Duplicate wells for each nitrogen condition |

| Columns | Nitrogen condition |
|---|---|
| 1 | Glutamate (`Glut`) |
| 2 | Ammonium sulfate (`AmSu`) |
| 3 | Urea |
| 4 | Aspartate (`Asp`) |
| 5 | Asparagine (`Asn`) |
| 6 | Glutamine (`Gln`) |
| 7 | N-acetylglucosamine (`GNAc`) |
| 8 | Glucosamine (`Glcn`) |
| 9 | Alanine (`Ala`) |
| 10 | Proline (`Pro`) |
| 11–12 | No added nitrogen (`NoN`) controls |

For every worksheet, the mean of wells G11, G12, H11, and H12 is used for blank correction. Blank-corrected values are rounded to three decimal places and clipped to a lower bound of 0.005 OD600. Because every `NoN` condition is a control, all rows for which `N == "NoN"` are removed from the analysis-ready dataset.

## Generated-table documentation

### `Generated Tables/big_norm_df.csv`

This is the reconstructed pre-filtering table. It may include `NoN` control rows so that the plate processing remains traceable.

| Column | Definition |
|---|---|
| `Environment` | Combined nitrogen and carbon label from the sample map |
| `Species` | `Bs` for *B. subtilis* or `Ec` for *E. coli* |
| `Pair` | Strain-pair or control identity |
| `N` | Nitrogen-source label |
| `C` | Carbon-source label |
| `Aux1`, `Aux2` | First and second auxotroph/control labels parsed from `Pair` |
| `Mean` | Mean blank-corrected endpoint OD600 for the species, environment, and pair |
| `SEM` | Standard error of `Mean` across the available measurements |
| `mono_aux1`, `mono_aux2` | Auxotroph labels used to join each partner to its unsupplemented monoculture measurement |
| `Aux_OD_for_norm_mean1`, `Aux_OD_for_norm_mean2` | Mean unsupplemented monoculture OD600 for partners 1 and 2 in the same species and environment |
| `Aux_OD_for_norm_sem1`, `Aux_OD_for_norm_sem2` | SEM of the corresponding unsupplemented monoculture measurements |
| `Auxs_growth_wo_AA` | Sum of `Aux_OD_for_norm_mean1` and `Aux_OD_for_norm_mean2` |
| `Norm_mean` | `Mean - Auxs_growth_wo_AA`, clipped to a minimum of 0.005; when matched monoculture values are unavailable, `Mean` is retained |
| `Norm_sem` | Propagated SEM: `sqrt(SEM^2 + Aux_OD_for_norm_sem1^2 + Aux_OD_for_norm_sem2^2)`; unavailable values are set to 0 |
| `NS`, `CS` | Nitrogen and carbon labels copied from the supplemented-control table during merging |
| `WT_mean`, `WT_sem` | Mean and SEM for the prototroph control in the same species and environment |
| `M_mean`, `M_sem` | Mean and SEM for the methionine-supplemented `M` monoculture |
| `F_mean`, `F_sem` | Mean and SEM for the phenylalanine-supplemented `F` monoculture |
| `R_mean`, `R_sem` | Mean and SEM for the arginine-supplemented `R` monoculture |
| `C_mean`, `C_sem` | Mean and SEM for the cysteine-supplemented `C` monoculture |
| `S_mean`, `S_sem` | Mean and SEM for the serine-supplemented `S` monoculture |
| `L_mean`, `L_sem` | Mean and SEM for the leucine-supplemented `L` monoculture |
| `max_aux_mean`, `max_aux_sem` | Mean and SEM of the better-growing supplemented control corresponding to the two members of `Pair` |
| `min_aux_mean`, `min_aux_sem` | Mean and SEM of the lower-growing supplemented control corresponding to the two members of `Pair` |

### `Generated Tables/filtered_df.csv`

This is the main analysis-ready table. It contains the columns in `big_norm_df.csv` plus the following filter fields.

| Column | Definition |
|---|---|
| `Aux1_bad` | `True` when partner 1 grows above 0.02 OD600 without amino-acid supplementation in that species–environment combination |
| `Aux2_bad` | `True` when partner 2 grows above 0.02 OD600 without amino-acid supplementation in that species–environment combination |
| `filtered_out` | Logical OR of `Aux1_bad` and `Aux2_bad` |

Only rows with `filtered_out == False` and `N != "NoN"` are exported to `filtered_df.csv`. Consequently, the retained rows normally have `False` in all three filter columns. The table is used for cross-feeding analyses and as the source of the random-forest model dataset.

### `Generated Tables/Supplemented_AA_data.csv`

| Column | Definition |
|---|---|
| `Environment` | Combined nitrogen and carbon label |
| `Pair` | Prototroph or amino-acid-supplemented monoculture identity |
| `Species` | `Bs` or `Ec` |
| `Mean` | Mean blank-corrected endpoint OD600 |
| `SEM` | Standard error of the endpoint OD600 |
| `Aux1`, `Aux2` | Components of `Pair` |
| `N`, `C` | Nitrogen and carbon labels parsed from `Environment` |

### `Generated Tables/AA_requirments_dataset.csv`

| Column | Definition |
|---|---|
| `Aux` | Auxotroph label |
| `Environment` | Nutrient environment used for the titration |
| `Num_of_points` | Number of data points included in the fitted titration relationship |
| `R2` | Coefficient of determination for the fitted relationship |
| `Slope` | Fitted slope supplied by the processed amino-acid requirement input |
| `C`, `N` | Carbon and nitrogen labels |
| `Species` | `Bs` or `Ec` |
| `Exp` | Experiment number |
| `aa_per_od` | Estimated extracellular amino-acid requirement per OD unit, calculated as `1 / Slope` |

### `Generated Tables/Koo_auxs_df_bs.csv`

| Column or group | Definition |
|---|---|
| `locus_tag` | *B. subtilis* locus identifier supplied by Koo et al. |
| `gene` | Gene name |
| 13 environment columns | Normalized relative-fitness values for the selected Koo environments |
| `Casette` | Resistance-cassette background, `Erm` or `Kan`; spelling retained from the analysis file |
| `n_grows_01` | Percentage of the 13 environments in which the value is greater than 0.1 |
| `n_grows_05` | Percentage of the 13 environments in which the value is greater than 0.5 |
| `n_grows_1` | Percentage of the 13 environments in which the value is greater than 1.0 |
| `rank_01`, `rank_05`, `rank_1` | Descending ranks corresponding to the three percentages; ties are resolved by first appearance |

Only pathway-associated genes with a value below 0.1 in at least one selected environment are retained.

### `Generated Tables/Tong_aux_df_ec.csv`

| Column or group | Definition |
|---|---|
| `B-numbers` | *E. coli* b-number identifier supplied by Tong et al. |
| `gene` | Gene name |
| 30 environment columns | Endpoint biomass relative to `kdsC` in the corresponding environment |
| `n_grows_01` | Percentage of the 30 environments in which relative growth is greater than 0.1 |
| `n_grows_05` | Percentage of the 30 environments in which relative growth is greater than 0.5 |
| `n_grows_1` | Percentage of the 30 environments in which relative growth is greater than 1.0 |
| `rank_01`, `rank_05`, `rank_1` | Descending ranks corresponding to the three percentages; ties are resolved by first appearance |

Normalization to `kdsC` is performed **before** applying the below-0.1 auxotrophy filter. Only pathway-associated genes with normalized growth below 0.1 in at least one environment are retained.

### `Generated Tables/Wetmore_aux_df_ec.csv`

| Column or group | Definition |
|---|---|
| `locusId` | Wetmore locus identifier |
| `sysName` | Systematic gene identifier |
| `gene` | Gene name joined from `fit_genes.tab` |
| `desc` | Gene description |
| Environment columns | Mean gene-fitness log2 ratio for each quality-approved, non-time-zero condition after experiments sharing the same short name are averaged |

Only genes represented in the curated amino-acid biosynthesis gene list are retained.

### `Generated Tables/Random Forest Model/`

The model directory contains small tabular outputs used to reproduce the model figures. A serialized model may be stored locally for convenience but is not required in the public repository.

| File | Columns and meaning |
|---|---|
| `perm_test_df.csv` | `feature`, `importance_mean`, and `importance_std`: held-out test-set permutation importance summarized across repeats |
| `perm_cv_df.csv` | `feature`, `importance_mean`, and `importance_std`: grouped cross-validation permutation importance summarized across folds |
| `block_cv_fast.csv` | `block`, `rmse_increase_mean`, `rmse_increase_std`, and `n_folds`: summary of feature-block permutation importance |
| `block_cv_folds_fast.csv` | `fold`, `block`, `rmse_increase`, and `n_test`: mean importance for each feature block in each held-out fold |
| `block_cv_perms.csv` | `fold`, `block`, `repeat`, `rmse_increase`, `base_rmse`, `rmse_perm`, and `n_test`: individual block-permutation results |

## Abbreviations

### Species

| Label | Meaning |
|---|---|
| `Bs` | *Bacillus subtilis* |
| `Ec` | *Escherichia coli* |

### Auxotrophs and controls

| Label | Meaning |
|---|---|
| `M` | Methionine auxotroph |
| `F` | Phenylalanine auxotroph |
| `R` | Arginine auxotroph |
| `C` | Cysteine auxotroph |
| `S` | Serine auxotroph |
| `L` | Leucine auxotroph |
| `WT` | Prototroph control |

### Carbon sources

| Label | Meaning |
|---|---|
| `Fum` | Fumarate |
| `Gluc` | Glucose |
| `Gly` | Glycerol |
| `NoC` | No added carbon source |

### Nitrogen sources

| Label | Meaning |
|---|---|
| `Ala` | Alanine |
| `Asn` | Asparagine |
| `Asp` | Aspartate |
| `Glut` | Glutamate |
| `Gln` | Glutamine |
| `Pro` | Proline |
| `GNAc` | N-acetylglucosamine |
| `Glcn` | Glucosamine |
| `AmSu` | Ammonium sulfate |
| `Urea` | Urea |
| `NoN` | No added nitrogen control; excluded from all analyses |

## Reproducibility notes

- Run the notebook from the repository root because input and output paths are relative to that location.
- Keep the current directory names unchanged unless the corresponding paths are also updated in the notebook.
- `NoN` controls must be excluded before every downstream statistical analysis, model fit, and figure described as covering 40 environments.
- The model groups observations by species and environment so that the same species–environment combination does not occur in both training and held-out test sets.
- A cached serialized model is not a substitute for the processed inputs, split settings, selected hyperparameters, model metrics, and permutation-output tables.
- Internet access is currently needed only for the BiGG mapping used in the Monk heatmap. Saving the mapping as a small local source table is recommended for fully offline reproducibility.
- Remove `.ipynb_checkpoints/`, `desktop.ini`, and other operating-system or notebook temporary files before creating the public repository.
- Figure panels are written to `images/` by the notebook.

## Citation

If you use these data or code, please cite the associated manuscript and the original publications corresponding to the external datasets. Publication identifiers can be added here when the manuscript citation is finalized.
