package gspa.model

import groovy.transform.Canonical

/**
 * Represents a consistency violation found during quality assessment.
 * Typically a taxonomic constraint violation identified by the SAT solver.
 */
@Canonical
class ConsistencyViolation {

    /** Type of violation */
    ViolationType type

    /** Severity level */
    Severity severity

    /** Human-readable description */
    String description

    /** Proteins involved in the violation */
    List<String> involvedProteins = []

    /** GO terms involved in the violation */
    List<String> involvedGoTerms = []

    /** Suggested corrective action */
    String suggestedAction

    /** SAT solver UNSAT core or OWL justification */
    String justification

    enum ViolationType {
        TAXON_CONFLICT,         // never_in_taxon / only_in_taxon violation
        MISSING_ESSENTIAL,      // Essential function not covered
        INCOHERENT_PATHWAY,     // Incomplete pathway
        INCOHERENT_COMPLEX,     // Missing complex subunits
        INCOHERENT_PROCESS,     // Missing has_part dependency
        CONTRADICTORY_FUNCTION  // Contradictory annotations on same protein
    }

    enum Severity {
        ERROR,   // Definitive violation (e.g., taxon constraint)
        WARNING, // Likely issue (e.g., incomplete pathway)
        INFO     // Informational (e.g., low IC annotation)
    }
}
