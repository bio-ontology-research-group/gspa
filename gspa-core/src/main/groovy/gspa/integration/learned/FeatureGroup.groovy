package gspa.integration.learned

/**
 * Provenance family of a learned-model feature. The separation is part of the
 * contract: baseline sequence features are computed from raw predictor claims only,
 * context features from genome/neighbourhood/cohort information, and neither ever
 * reads a posterior produced by the integrator (that would feed a conclusion back in
 * as independent evidence).
 */
enum FeatureGroup {
    /** Function type / GO aspect / complex-term membership of the claim itself. */
    FUNCTION_IDENTITY,
    /** Candidate role (catalytic, accessory, transporter, ...), a distinct dimension of identity. */
    ROLE,
    /** Source-wise and evidence-group-wise strongest raw predictor scores (no context). */
    BASELINE_SEQUENCE,
    /** Missingness of predictor runs (failed / unavailable / partial tools). */
    RUN_COVERAGE,
    /** Operon, gene-neighbourhood, local competitors, listed metabolic gaps, supplied neighbourhood. */
    CONTEXT_NEIGHBORHOOD,
    /** Cross-genome recurrence supplied by a cohort analysis. */
    CONTEXT_COHORT,
    /** Contig/assembly properties that make absence of evidence uninformative. */
    ASSEMBLY,
    /** Explicitly conditional constraint-based (FBA) features. */
    FBA_CONDITIONAL
}
