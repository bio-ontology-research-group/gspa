package gspa.integration.cohort

/**
 * Functional role of a protein within a metabolic neighbourhood. Roles are
 * resolved from protein-intrinsic evidence only; genomic or cohort context
 * never assigns {@link #CATALYTIC}.
 */
enum Role {
    CATALYTIC,
    ACCESSORY,
    TRANSPORTER,
    REGULATOR,
    COMPLEX_SUBUNIT,
    UNKNOWN
}
