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

    /** Claims (GO or EC) directly attributed to the violation */
    List<String> involvedClaims = []

    /** For EC violations, maps the EC claim ID to the mapped GO term(s) that caused the violation */
    Map<String, List<String>> ecToGoCauses = [:]

    /** Suggested corrective action */
    String suggestedAction

    /** SAT solver UNSAT core or OWL justification */
    String justification

    enum ViolationType {
        TAXON_CONFLICT,         // never_in_taxon / only_in_taxon violation
        ANATOMICAL_CONFLICT,    // GO-Uberon-PO anatomical/developmental structure violation
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
