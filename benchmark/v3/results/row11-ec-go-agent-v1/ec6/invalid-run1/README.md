First scoring pass (2026-11-03), archived as INVALID and not used for any result. Three defects found after this pass
(see the plan's Deviations, 2026-11-03, second entry): (1) the shared guardrail failed closed on every PROMOTE whose
EC maps to a GO term because the Price organisms' NCBI taxids are not in the bundled taxon hierarchy
(UNKNOWN_TAXON), whereas the production engine does not apply taxon constraints in that state; (2) the locus agent
never forced a resolution on its last turn, so every case abstained with "did not emit a resolution"; (3) policy
files reported 0 known calibration decisions whenever no threshold qualified. The public bundles, PUBLIC_SHA256 and
truth files are unchanged; only methods, policies and evaluations were rerun.
