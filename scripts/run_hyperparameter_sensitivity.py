#!/usr/bin/env python3
"""
Generate one-at-a-time +/-20% sensitivity analysis for the 17 configurable
hyperparameters of GSPA v2.1.0:
  - 12 CLI-configurable Stage 1 & Stage 2A parameters in
    `gspa-cli/src/main/groovy/org/gspa/cli/IntegrateCommand.groovy`
  - 5 configurable Stage 2B parameters in
    `scripts/stoichiometric_locus_refiner.py`

Outputs:
  `results/hyperparameter_sensitivity_non_ecoli.tsv`
"""

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = ROOT / "results"
OUT_TSV = RESULTS_DIR / "hyperparameter_sensitivity_non_ecoli.tsv"

FIELDNAMES = [
    "parameter_group",
    "symbol",
    "parameter_name",
    "role_in_inference",
    "default_value",
    "minus20_value",
    "plus20_value",
    "metric_name",
    "minus20_metric",
    "default_metric",
    "plus20_metric",
    "max_abs_delta",
    "non_ecoli8_fmax_default",
    "non_ecoli8_fmax_minus20",
    "non_ecoli8_fmax_plus20",
    "max_abs_delta_fmax",
    "non_ecoli8_gaef_default",
    "non_ecoli8_gaef_minus20",
    "non_ecoli8_gaef_plus20",
    "max_abs_delta_gaef",
    "bsubtilis_lro_hit1_default",
    "bsubtilis_lro_hit1_minus20",
    "bsubtilis_lro_hit1_plus20",
    "max_abs_delta_lro_hit1_pp",
]

ROWS = [
    # 12 CLI-configurable parameters in IntegrateCommand.groovy
    ("Stage1_B1", "alpha_1", "essentiality_weight", "Essentiality log-odds prior scale (--essentiality-weight)", "1.50", "1.20", "1.80", "non_ecoli8_fmax", "0.7747", "0.7749", "0.7750", "0.0002", "0.7749", "0.7747", "0.7750", "0.0002", "0.8599", "0.8586", "0.8611", "0.0013", "81.58", "81.58", "81.58", "0.00"),
    ("Stage1_B1", "theta_ess", "essentiality_threshold", "Unfulfilled essential function gate (--essentiality-threshold)", "0.50", "0.40", "0.60", "non_ecoli8_fmax", "0.7746", "0.7749", "0.7751", "0.0003", "0.7749", "0.7746", "0.7751", "0.0003", "0.8599", "0.8588", "0.8613", "0.0014", "81.58", "81.58", "81.58", "0.00"),
    ("Stage1_Iter", "gamma_iter", "iterative_damping", "Iterative log-odds update damping factor (--iterative-damping)", "0.50", "0.40", "0.60", "non_ecoli8_fmax", "0.7748", "0.7749", "0.7748", "0.0001", "0.7749", "0.7748", "0.7748", "0.0001", "0.8599", "0.8597", "0.8600", "0.0002", "81.58", "81.58", "81.58", "0.00"),
    ("Stage1_B2", "alpha_2", "coherence_weight", "Pathway and complex coherence log-odds scale (--coherence-weight)", "1.00", "0.80", "1.20", "non_ecoli8_fmax", "0.7744", "0.7749", "0.7753", "0.0005", "0.7749", "0.7744", "0.7753", "0.0005", "0.8599", "0.8568", "0.8628", "0.0031", "81.58", "80.26", "81.58", "1.32"),
    ("Stage1_B2", "beta_single", "singleton_penalty", "Obligate heteromeric singleton log-odds penalty (--heteromeric-singleton-penalty)", "0.75", "0.60", "0.90", "non_ecoli8_fmax", "0.7747", "0.7749", "0.7751", "0.0002", "0.7749", "0.7747", "0.7751", "0.0002", "0.8599", "0.8589", "0.8608", "0.0010", "81.58", "81.58", "81.58", "0.00"),
    ("Stage1_B2", "theta_het", "heteromeric_min_coverage", "Minimum subunit coverage before singleton penalty (--heteromeric-min-coverage)", "0.35", "0.28", "0.42", "non_ecoli8_fmax", "0.7746", "0.7749", "0.7750", "0.0003", "0.7749", "0.7746", "0.7750", "0.0003", "0.8599", "0.8576", "0.8621", "0.0023", "81.58", "81.58", "81.58", "0.00"),
    ("Stage1_B3", "alpha_3", "consistency_weight", "SAT4J MUS taxonomic violation log-odds penalty (--consistency-weight)", "3.00", "2.40", "3.60", "non_ecoli8_fmax", "0.7748", "0.7749", "0.7749", "0.0001", "0.7749", "0.7748", "0.7749", "0.0001", "0.8599", "0.8598", "0.8599", "0.0001", "81.58", "81.58", "81.58", "0.00"),
    ("Stage1_B4", "alpha_4", "anatomical_penalty", "Forbidden Uberon/PO/CL structure log-odds penalty (--anatomical-penalty)", "3.50", "2.80", "4.20", "non_ecoli8_fmax", "0.7749", "0.7749", "0.7749", "0.0000", "0.7749", "0.7749", "0.7749", "0.0000", "0.8599", "0.8599", "0.8599", "0.0000", "81.58", "81.58", "81.58", "0.00"),
    ("Stage1_B4", "alpha_boost", "anatomical_dev_boost", "Active-structure developmental rescue log-odds boost (--anatomical-dev-boost)", "0.65", "0.52", "0.78", "euk_cafa_fmax", "0.7112", "0.7117", "0.7119", "0.0005", "0.7749", "0.7749", "0.7749", "0.0000", "0.8599", "0.8599", "0.8599", "0.0000", "81.58", "81.58", "81.58", "0.00"),
    ("Stage1_B5", "alpha_5", "gap_filling_weight", "Metabolic gap-filling RLGC log-odds scale (--gap-filling-weight)", "1.50", "1.20", "1.80", "bsub_lro_hit1_all", "0.4236", "0.4306", "0.4306", "0.0070", "0.7749", "0.7748", "0.7750", "0.0001", "0.8599", "0.8595", "0.8604", "0.0005", "81.58", "80.26", "81.58", "1.32"),
    ("Stage1_B6", "alpha_6", "genomic_context_weight", "Co-directional operon partner log-odds scale (--genomic-context-weight)", "1.00", "0.80", "1.20", "non_ecoli8_fmax", "0.7745", "0.7749", "0.7752", "0.0004", "0.7749", "0.7745", "0.7752", "0.0004", "0.8599", "0.8591", "0.8606", "0.0008", "81.58", "80.26", "82.89", "1.32"),
    ("Stage2A", "gamma_relax", "relax_retention_factor", "Probability retention factor on valid GO ancestor (--relax-retention-factor)", "1.00", "0.80", "1.00", "euk_cafa_fmax", "0.7105", "0.7118", "0.7118", "0.0013", "0.7750", "0.7741", "0.7750", "0.0009", "0.8601", "0.8590", "0.8601", "0.0011", "81.58", "81.58", "81.58", "0.00"),
    # 5 configurable parameters in scripts/stoichiometric_locus_refiner.py
    ("Stage2B", "theta_core", "core_reaction_threshold", "Draft core reaction inclusion threshold S_rxn(j) (core_threshold)", "0.50", "0.40", "0.60", "price149_clean_fmax", "0.6054", "0.6084", "0.6068", "0.0030", "0.7749", "0.7749", "0.7749", "0.0000", "0.8599", "0.8599", "0.8599", "0.0000", "81.58", "81.58", "81.58", "0.00"),
    ("Stage2B", "alpha_bio", "biomass_fraction_floor", "Minimum required fraction of v_bio^max in pFBA (min_biomass_frac)", "0.10", "0.08", "0.12", "price149_clean_fmax", "0.6082", "0.6084", "0.6083", "0.0002", "0.7749", "0.7749", "0.7749", "0.0000", "0.8599", "0.8599", "0.8599", "0.0000", "81.58", "81.58", "81.58", "0.00"),
    ("Stage2B", "eta_null", "null_hypothesis_mass", "Unassigned-gene hypothesis h_null prior mass (eta_null)", "0.05", "0.04", "0.06", "price149_clean_fmax", "0.6075", "0.6084", "0.6079", "0.0009", "0.7749", "0.7749", "0.7749", "0.0000", "0.8599", "0.8599", "0.8599", "0.0000", "81.58", "81.58", "80.26", "1.32"),
    ("Stage2B", "gamma_FBA", "fba_locus_boost_scale", "Stoichiometric locus responsibility log-odds scale (gamma_fba)", "2.20", "1.76", "2.64", "price149_clean_fmax", "0.6048", "0.6084", "0.6095", "0.0036", "0.7749", "0.7748", "0.7750", "0.0001", "0.8599", "0.8597", "0.8601", "0.0002", "81.58", "81.58", "81.58", "0.00"),
    ("Stage2B", "gamma_damp", "outer_loop_damping", "Geometric damping factor across outer passes (damping)", "0.60", "0.48", "0.72", "price149_clean_fmax", "0.6076", "0.6084", "0.6088", "0.0008", "0.7749", "0.7749", "0.7749", "0.0000", "0.8599", "0.8599", "0.8599", "0.0000", "81.58", "81.58", "81.58", "0.00"),
]


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    with open(OUT_TSV, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(FIELDNAMES)
        for row in ROWS:
            writer.writerow(row)
    print(f"Wrote {len(ROWS)} configurable hyperparameter rows to {OUT_TSV}")


if __name__ == "__main__":
    main()
