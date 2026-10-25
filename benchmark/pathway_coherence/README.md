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
