# Measured nutrient dependence and supplementation

This adds a raw-data audit and outcome-independent inputs for the biological
validation required by audit items 8, 9 and 12. It does not report model accuracy
or accept those items as complete. The earlier
[absence comparison](v3-absence-comparison.md) was validated on synthetic
mechanisms; these observations provide real matched environmental interventions.

## Primary sources and scope

The source is version 3 of Baichman-Kass, Noda-Garcia and Friedman,
[Nutrient environments shape amino acid auxotrophy and cross-feeding](https://www.biorxiv.org/content/10.64898/2026.05.20.726445v3),
a preprint posted on 6 August 2026. Its XML and PDF were downloaded and the XML
identity/version checked. The author's repository is pinned to commit
[`ee3a2079a4a7e47c090a7ad9945b053c72723490`](https://github.com/amichaibk/nutrient-dependent-auxotrophy_and_CF/tree/ee3a2079a4a7e47c090a7ad9945b053c72723490).
Thirteen input files retain Git blob hashes, SHA-256 hashes and byte sizes.
The preparation code reads the notebook as source evidence; it never executes it.

The experiment uses six auxotrophs in each of two species:

| Label | Required amino acid | *E. coli* deletion | *B. subtilis* deletion |
|---|---|---|---|
| M | Methionine | metA | metA |
| F | Phenylalanine | pheA | pheA |
| R | Arginine | argA | argH |
| C | Cysteine | cysE | cysE |
| S | Serine | serA | serA |
| L | Leucine | leuB | leuB |

The *E. coli* background is EcNR1, derived from MG1655; the *B. subtilis*
background is 168. Those identities come from the paper's Methods, not from
auxotroph letters or presumed gene orthology. Inputs keep each coculture organism
separate; a pair is never encoded as one strain with two deletions.

The primary design has ten named nitrogen sources and four carbon conditions
(glucose, glycerol, fumarate and no added carbon), giving 40 environments. Four
no-added-nitrogen environments remain as controls. “No added carbon” does not
mean chemically carbon-free: the base medium contains tryptophan and some named
nitrogen sources also contain carbon.

The exported protocol preserves the printed MSgg ingredients, concentrations,
30°C incubation, 96-hour endpoint, 200-µl volume and initial OD600 of 0.005.
Supplementation is 200 µM for five amino acids and 1 mM for serine. Concentrations
are not converted into invented uptake fluxes. Oxygen transfer, maintenance,
biomass and model-specific transport remain explicit unresolved modeling inputs.

## Raw measurements and independent reconstruction

`prepare_nutrient_environments.py` locates the labelled 8×12 grid in each workbook
sheet, checks the complete sample map, and records the file, sheet and cell for
each observation. It preserves raw OD600, the four-well plate background,
rounded corrected OD and the author's 0.005 floor. Missing values and sheets
stay unavailable; they are never converted to no growth.

The full acquisition contains 16,608 measured wells, including blank sheets.
The case catalogue is fixed by assay identities, independently of outcomes:

| Assay | All environments | Primary environments |
|---|---:|---:|
| Unsupplemented monoculture | 528 | 480 |
| Supplemented monoculture | 528 | 480 |
| Coculture | 1,320 | 1,200 |
| Wild type | 88 | 80 |
| Total | 2,464 | 2,240 |

The F+S sheet is absent from *E. coli* experiment 2. Its cases retain the other
two experiments, explicitly recording the missing experiment and wells. A
wholly missing assay would still remain in the catalogue with unknown outcome.
No case is excluded because its measured phenotype disagrees with a model or
because an auxotroph grew without supplementation. No-nitrogen controls are
tagged separately, not discarded.

The output reports both the authors' pooled-well SEM and an SEM across experiment
means. This keeps technical wells distinguishable from experiment-level repeats.
The author's OD600>0.02 monoculture threshold is retained as an observation flag;
it is not a calibrated function score or proof of gene activity.

## Export discrepancy

Independent reconstruction agrees with all 1,936 means and pooled SEMs in
`big_norm_df.csv`, and all 616 means in `Supplemented_AA_data.csv`.
Every one of the latter file's 616 SEM entries instead equals its mean and
disagrees with the raw-replicate SEM. Notebook cell 15 constructs that SEM column
from `Mean_Bs` and `Mean_Ec` again after discarding the actual SEM column.

The original files are preserved unchanged. The audit lists every discrepancy,
and the new observation output uses the uncertainty recomputed from the raw
wells. It does not copy those erroneous exported uncertainties or silently edit
the author's source.

## Model preflight

`register_nutrient_models.py` reads only the outcome-free case input file and
supplied SBML models. It records exact gene-symbol-to-locus mappings, canonical
gene identities, model GPR effects, objectives, native maintenance and hashes.
Model GPRs are used for known-genotype phenotype evaluation, never as independent
gene-function truth.

- iML1515: the previously acquired author-repository model,
  SHA-256 `9c772d44ca43350e40dc7ee86c7aa148796856be1eea45e5406c6df8f7dcde28`.
- [iYO844](https://bigg.ucsd.edu/models/iYO844): downloaded from the primary BiGG
  distribution for *B. subtilis* 168,
  SHA-256 `28e4ca32d412102481ac29f84b0267b1449108bbf1fa6efc604794b4a93b9740`.

All twelve study gene symbols map exactly. In iML1515, deleting serA blocks three
reactions, while deleting leuB blocks none because its reaction has an alternative
GPR route. The other four *E. coli* deletions and all six *B. subtilis* deletions
each block one reaction. The preflight retains those differences rather than
forcing all genes into single-reaction absence scenarios. The EcNR1/model-strain
difference and these model assumptions remain part of the eventual interpretation.

## Reproduction and acceptance status

```bash
python benchmark/v3/prepare_nutrient_environments.py acquire --root "$DATA"
python benchmark/v3/prepare_nutrient_environments.py prepare \
  --root "$DATA" --out "$DATA/derived-new"
python benchmark/v3/register_nutrient_models.py \
  --inputs "$DATA/derived-new/inputs.jsonl" \
  --ec-model "$EC_MODEL" --bs-model "$BS_MODEL" --out "$DATA/model-preflight-new.json"
```

`DATA`, `EC_MODEL` and `BS_MODEL` name local acquisition/model paths. Existing
derived outputs are not overwritten. `measurements.jsonl` retains raw provenance;
`observations.jsonl` contains measured summaries; `inputs.jsonl` contains only
genotypes and conditions. The report binds every source and output. Tests prove
that changing OD values or measurement coverage cannot change simulation inputs.

Thirteen focused checks pass, covering both workbook layouts, malformed/duplicate
grids, blank-well identities, explicit missingness, experiment versus pooled SEM,
the author export defect, source/Git identity, complete input catalogues,
species-specific genotypes, separated coculture organisms and canonical GPR
mapping with an isozyme alternative. Initial acquisition metadata errors and
passing final logs are retained.

The complete benchmark gate passes 626 cases, with three optional skips and ten
executed live Nextflow cases. JVM production code is unchanged in this slice;
the preceding 1,251-case JVM gate remains the current engine verification.
The raw source, every derived version, model files, exact commands, tests and
hash receipts are retained in `benchmark/v3/results/nutrient-environments-v1/`.

This acquisition audit makes no biological accuracy claim. The subsequent
[matched phenotype benchmark](v3-nutrient-phenotypes.md) freezes six medium/oxygen
arms, executes all 1,144 single-organism assays and joins these measurements only
after prediction. Its measured rescue accuracy remains below a constant-rescue
comparator. Biological absence discrimination, separate gene attribution,
community and actual sequence-loss evaluation, transporter knockout validation
and the other original acceptance gates remain open.
