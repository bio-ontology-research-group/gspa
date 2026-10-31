# Secondary metabolism: Sco-GEM source intake

This intake identifies usable primary material and reasons it is not yet an
accepted benchmark. **No biological cases or function labels have been accepted.**

Source: Sulheim et al. (2020), *Enzyme-Constrained Models and Omics Analysis of
Streptomyces coelicolor Reveal Metabolic Changes that Enhance Heterologous
Production*, [DOI 10.1016/j.isci.2020.101525](https://doi.org/10.1016/j.isci.2020.101525),
PMID 32942174, PMCID PMC7501462. The full article XML was retrieved from Europe
PMC and its identifiers verified. The authors' recommended paper reproduction
release is [Sco-GEM v1.2.1](https://github.com/SysBioChalmers/Sco-GEM/tree/v1.2.1),
commit `63003f8b01a692dc2b73f2c9e0bbde871a3fe6f0`.

All **26** acquired repository files match the pinned Git tree blob hashes and
their SHA-256 receipts. They include the base SBML model, media, growth and
production measurements, estimated rates and fitting scripts. The original
CC BY 4.0 license and authors' README are retained. Downloading and hashing these
files establishes provenance; it does not resolve experimental identities.

The two offline tables contain **52 rows each**. Their product columns preserve
the following original observation types; blanks are not converted to zero:

| Directory / product column | Numeric | Below detection limit | Missing |
|---|---:|---:|---:|
| M145 / Undecylprodigiosin 2 | 6 | 4 | 42 |
| M145 / TBP | 10 | 0 | 42 |
| M145 / Germicidin-A | 8 | 2 | 42 |
| M145 / Germicidin-B | 9 | 1 | 42 |
| M1152 / Undecylprodigiosin 2 | 0 | 3 | 49 |
| M1152 / TBP | 4 | 0 | 48 |
| M1152 / Germicidin-A | 2 | 1 | 49 |
| M1152 / Germicidin-B | 3 | 0 | 49 |

Both files carry the same M145/F134 header, including the file under M1152.
The directory names therefore cannot by themselves establish strain identity.
TBP uses `spec units`; it is not silently mapped to a molar actinorhodin flux.
Exact medium composition, phosphate-dependent sampling phases, replicate
handling and assay/strain identity still need reconciliation with the methods.
The article describes M1152 as a four-BGC deletion derivative with an `rpoB`
mutation; its phenotype is not a single-enzyme deletion experiment.

The source fitting script
[`fixFluxes.m`](https://github.com/SysBioChalmers/Sco-GEM/blob/v1.2.1/ComplementaryScripts/ecModel/simulation/fixFluxes.m)
uses estimated RED and germicidin rates to set narrow flux bounds. Testing those
fitted models against the same production measurements would be circular. A
future benchmark must freeze the base model and permitted nutrient inputs,
withhold product observations, and document any model-development overlap.
Dynamic concentration measurements also require an explicit conversion or
comparison protocol before use with steady-state FBA.

The base model has external actinorhodin and gamma-actinorhodin exchanges, but
RED and germicidin A/B are represented by **intracellular demand reactions**.
Those sinks do not establish extracellular secretion or transport. The current
GSPA product menu correctly requires an external exchange; relabelling sinks as
exchanges to make this dataset run would change the biological question. Any
extension to intracellular synthesis feasibility needs its own explicit output
semantics and validation.

Next benchmark work is to resolve the source identities and medium, define
withheld product outcomes with detection limits, and compare condition-specific
predictions against those outcomes without using them as bounds. Function
recovery additionally needs independently characterized target genes and
masking of their assignments. Until then, this source collection is a research
intake, not evidence for an improvement in EC prediction or secondary metabolism.

Archive: `benchmark/v3/results/secondary-source-intake-v1/`. Its verifier checks
every retained byte and reproduces the source audit with
`benchmark/v3/audit_secondary_metabolism_sources.py`. This replay validates the
intake counts and model representation; it does not validate the unresolved
biological interpretations. Item 12 and all other whole items remain open.
