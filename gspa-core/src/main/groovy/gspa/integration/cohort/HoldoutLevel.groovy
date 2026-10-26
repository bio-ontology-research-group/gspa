package gspa.integration.cohort

/** Which genomes are removed from the peer set when scoring a target genome. */
enum HoldoutLevel {
    /** Only the target genome itself. */
    GENOME,
    /** The target genome and every genome of its strain cluster (default). */
    STRAIN_CLUSTER,
    /** The target genome and every genome of its clade. */
    CLADE,
    /** Nobody: descriptive cohort aggregates only, never for scoring a member. */
    NONE
}
