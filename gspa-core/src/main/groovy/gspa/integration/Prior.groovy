package gspa.integration

/**
 * A prior that adds a log-odds contribution to a candidate (protein,
 * function) assignment during iterative refinement.
 *
 * <p>Implemented in Phase 7.3. Phase 7.1/7.2 ship {@link PriorEngine}
 * with an empty {@code priors} list so refinement is likelihood-only.</p>
 */
interface Prior {
    String name()

    /**
     * Log-odds contribution for (proteinId, functionKey). Zero means no
     * effect; negative values downweight the claim; positive values boost.
     */
    double logOddsBoost(String proteinId, String functionKey, IntegrationState state)

    /** Declares which parts of state this prior reads, for re-run triggers. */
    Set<String> inputs()
}
