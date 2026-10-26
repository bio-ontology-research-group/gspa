package gspa.integration.cohort

/**
 * Whether a genome assembly may be treated as containing every gene of the
 * organism. {@link #UNKNOWN} is the explicit value for missing metadata;
 * completeness is never inferred.
 */
enum AssemblyStatus {
    COMPLETE,
    DRAFT,
    UNKNOWN
}
