# Pathway coherence: evidence-gated KEGG module completeness and module-based refinement

## Problem with the v2.1.0 metric

`PathwayDatabase.computePathwayCoherenceDetailed` treats a pathway as *triggered* when a single one of its
enzymes is annotated and then scores the fraction of its enzymes present. Shared enzymes of central
metabolism trigger dozens of pathways in every genome, KEGG maps are enzyme bags without step logic, and
an enzyme whose EC number has no GO mapping counts as missing. The number reports breadth of enzyme
coverage, not whether the organism can run a pathway. The coherence prior then boosted candidates for any
missing enzyme of any triggered pathway.

## What this adds

* **Module definitions with step logic** (`fetch_kegg_module_definitions.py` → `data/kegg_module_defs.json`):
  582 KEGG modules with their DEFINITION expressions parsed into required steps, alternative orthologs,
  complexes with optional subunits and KEGG's own gaps, plus KO → EC and the organisms KEGG lists as
  complete (`complete_organisms`). `data/kegg_organism_modules.json` holds gene-to-module links for the
  nine validation organisms. 2,549 steps; 69 % map to GO through ec2go, the rest are *unmappable* and are
  reported as such instead of as missing functions.
* **Assessment** (`gspa.ontology.ModuleDefinitions`, `gspa.ontology.ModuleAssessor`): per step a status
  (present, partial complex, sub-threshold evidence, missing, unmappable, forbidden for the taxon), module
  completeness over assessable required steps, and an activation rule weighted by enzyme specificity
  (1 / number of modules requiring the enzyme) so shared enzymes activate nothing on their own. Signature
  modules are assessed but never scored. `gspa evaluate --module-defs` adds `module_coherence`,
  `modules_activated`, `modules_complete` and `module_detail` to the quality JSON and keeps the old number
  under `pathway_coherence` and `pathway_map_coverage` as a breadth measure.
* **Refinement** (`CoherencePrior`, mode `module`): only activated modules at least half complete are
  refined, one candidate per open step and at most two steps per module; a candidate is a sub-threshold
  claim on a step enzyme scored by Stage 0 log-odds, corroboration, an unannotated-protein bonus, a penalty
  for a protein already confidently annotated elsewhere in the same aspect, an operon-partner bonus and the
  module's completeness. The chosen claim is pinned to 0.60. Steps without a candidate stay open.
* **Provenance**: `gspa integrate --module-defs` writes `<out>.modules.json` with every module before and
  after refinement, the proteins behind each present step, the annotations that changed with their
  evidence, and each unresolved gap with its best sub-threshold candidate and the reason it was left
  alone; rescues also go to `<out>.rescues.tsv`.

## Validation (`run_bench9.py`, unimatrix01 `/data/hohndor/gspa/pathway_coherence/`)

Nine reference genomes, three variants (`legacy`, `assess` = definitions loaded without module
refinement, `module`), judged by F-max against curated annotations, coherence-rescue precision against the
curated annotation of the rescued protein, and module-level agreement with the organisms KEGG lists as
complete. Results in `results/`.

## Round 1 findings (`results/` from job 9692, before the assessability rules)

Refinement was neutral on F-max (micro all 0.7734 to 0.7744, CAFA all 0.8057 to 0.8064, experimental
truth within 0.001) and made 331 coherence rescues across the nine genomes, 131 of the 324 with a curated
protein correct at term level (40 %). Module completeness judged against KEGG's complete-organism lists
reached recall 0.95 in E. coli but precision 0.50: 105 modules called complete that KEGG does not list.
Three causes, two of them fixable:

* 29 of the 105 had a single GO-mappable step and 26 had fewer than half their steps mappable, so one
  enzyme made them "complete". Fix: a module is *assessable* only with at least two mappable required
  steps covering at least half of its steps; the rest are reported as not assessable.
* Partial EC numbers (1.1.1.-) map through ec2go to enzyme-class terms (hydro-lyase, acid-thiol ligase,
  methyltransferase activity) that every genome carries, so coarse steps counted as present. Fix:
  partial ECs are treated as no mapping (337 of 4,107 orthologs).
* KEGG names organism-specific variants of one chemistry (isoprenoid biosynthesis for archaea, plants and
  non-plant eukaryotes; C1-unit interconversion, eukaryotes) whose orthologs differ below EC resolution,
  so a GO-based annotation satisfies all variants. Fix: the module name's domain hint against the genome's
  domain marks such modules as not expected; variants without a name hint remain indistinguishable and
  are a floor on the precision measured against KEGG.

Round 2 (`bench9_v2`) adds the `nocoh` baseline (all priors except coherence) so the annotations each
pathway refinement adds, and their precision against the curated truth, are measured directly.
