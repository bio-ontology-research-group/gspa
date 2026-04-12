package gspa.ontology

import gspa.model.ConsistencyViolation
import gspa.model.ConsistencyViolation.Severity
import gspa.model.ConsistencyViolation.ViolationType
import org.sat4j.core.VecInt
import org.sat4j.minisat.SolverFactory
import org.sat4j.specs.ContradictionException
import org.sat4j.specs.ISolver
import org.sat4j.specs.IVecInt
import org.sat4j.tools.xplain.Xplain
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * SAT-based consistency checker for taxon constraints.
 *
 * Encodes the taxon constraint satisfaction problem as a Boolean SAT instance:
 * - Variables: one per taxon node in the relevant portion of the NCBI taxonomy
 * - Clauses:
 *   - Taxonomy hierarchy: child → parent (if taxon A is child of B, then A implies B)
 *   - only_in_taxon T: unit clause requiring T to be true
 *   - never_in_taxon T: unit clause requiring T to be false
 *
 * If UNSAT, extracts minimal unsatisfiable core to identify conflicting GO annotations.
 */
class SatConsistencyChecker {

    private static final Logger log = LoggerFactory.getLogger(SatConsistencyChecker)

    TaxonConstraints taxonConstraints

    /** NCBI taxonomy parent map: child taxon ID -> parent taxon ID */
    Map<String, String> taxonomyParentMap = [:]

    /** Reverse map: parent taxon ID -> set of children */
    private Map<String, Set<String>> taxonomyChildrenMap = [:].withDefault { [] as Set }

    SatConsistencyChecker(TaxonConstraints taxonConstraints) {
        this.taxonConstraints = taxonConstraints
    }

    /**
     * Load the taxonomy hierarchy from a simple TSV file.
     * Format: child_taxon_id\tparent_taxon_id
     */
    void loadTaxonomyHierarchy(File taxonomyFile) {
        log.info("Loading taxonomy hierarchy from: ${taxonomyFile}")
        taxonomyFile.eachLine { line ->
            if (line.startsWith('#') || line.trim().isEmpty()) return
            def fields = line.split('\t')
            if (fields.length >= 2) {
                taxonomyParentMap[fields[0].trim()] = fields[1].trim()
            }
        }
        log.info("Loaded ${taxonomyParentMap.size()} taxonomy parent-child relationships")
        buildChildrenMap()
    }

    /**
     * Load taxonomy hierarchy from a map (for programmatic use / testing).
     */
    void loadTaxonomyHierarchy(Map<String, String> parentMap) {
        this.taxonomyParentMap = new LinkedHashMap<>(parentMap)
        buildChildrenMap()
    }

    private void buildChildrenMap() {
        taxonomyChildrenMap.clear()
        taxonomyParentMap.each { child, parent ->
            taxonomyChildrenMap[parent] << child
        }
    }

    /**
     * Check consistency of a set of GO annotations against taxon constraints.
     *
     * @param goTerms The set of GO term IDs annotated to proteins in the genome
     * @return A ConsistencyResult with satisfiability status and any violations
     */
    ConsistencyResult check(Set<String> goTerms) {
        def requirements = taxonConstraints.getRequirements(goTerms)

        if (requirements.positiveTaxa.isEmpty() && requirements.negativeTaxa.isEmpty()) {
            return new ConsistencyResult(consistent: true)
        }

        // Collect all relevant taxa (positive, negative, and their ancestors)
        Set<String> relevantTaxa = new HashSet<>()
        relevantTaxa.addAll(requirements.positiveTaxa)
        relevantTaxa.addAll(requirements.negativeTaxa)

        // Add ancestors of all relevant taxa
        Set<String> toExpand = new HashSet<>(relevantTaxa)
        while (!toExpand.isEmpty()) {
            Set<String> newTaxa = []
            toExpand.each { taxon ->
                String parent = taxonomyParentMap[taxon]
                if (parent && !relevantTaxa.contains(parent)) {
                    relevantTaxa.add(parent)
                    newTaxa.add(parent)
                }
            }
            toExpand = newTaxa
        }

        // Assign SAT variables
        Map<String, Integer> taxonToVar = [:]
        Map<Integer, String> varToTaxon = [:]
        int varCounter = 1
        relevantTaxa.each { taxon ->
            taxonToVar[taxon] = varCounter
            varToTaxon[varCounter] = taxon
            varCounter++
        }

        // Build SAT problem with explanation support
        Xplain<ISolver> solver = new Xplain<>(SolverFactory.newDefault())
        solver.newVar(varCounter)
        solver.setTimeout(60) // 60 second timeout

        // Track which clauses correspond to which GO annotations
        Map<Integer, String> clauseToSource = [:]

        try {
            // 1. Taxonomy hierarchy clauses: child → parent  (¬child ∨ parent)
            relevantTaxa.each { child ->
                String parent = taxonomyParentMap[child]
                if (parent && taxonToVar.containsKey(parent)) {
                    int childVar = taxonToVar[child]
                    int parentVar = taxonToVar[parent]
                    solver.addClause(new VecInt([-childVar, parentVar] as int[]))
                }
            }

            // 1b. Sibling disjointness: an organism can be in at most one child of a given parent.
            // For each parent with multiple children in our relevant set,
            // add pairwise mutual exclusion: ¬A ∨ ¬B for each sibling pair (A, B).
            taxonomyChildrenMap.each { parent, children ->
                def relevantChildren = children.findAll { taxonToVar.containsKey(it) }
                if (relevantChildren.size() >= 2) {
                    def childList = relevantChildren.toList()
                    for (int i = 0; i < childList.size(); i++) {
                        for (int j = i + 1; j < childList.size(); j++) {
                            int varA = taxonToVar[childList[i]]
                            int varB = taxonToVar[childList[j]]
                            solver.addClause(new VecInt([-varA, -varB] as int[]))
                        }
                    }
                }
            }

            // 2. Positive constraints: only_in_taxon T → must be in T (unit clause: T)
            requirements.positiveTaxa.each { taxon ->
                int var = taxonToVar[taxon]
                IVecInt clause = new VecInt([var] as int[])
                solver.addClause(clause)

                String goTerms_str = requirements.positiveSource[taxon]?.join(', ')
                clauseToSource[var] = "only_in_taxon(${taxon}) required by: ${goTerms_str}"
            }

            // 3. Negative constraints: never_in_taxon T → must NOT be in T (unit clause: ¬T)
            requirements.negativeTaxa.each { taxon ->
                int var = taxonToVar[taxon]
                IVecInt clause = new VecInt([-var] as int[])
                solver.addClause(clause)

                String goTerms_str = requirements.negativeSource[taxon]?.join(', ')
                clauseToSource[var] = "never_in_taxon(${taxon}) required by: ${goTerms_str}"
            }

        } catch (ContradictionException e) {
            // Contradiction detected during clause addition — definitely UNSAT
            log.debug("Contradiction detected during SAT encoding: ${e.message}")
            return buildUnsatResult(requirements, "Direct contradiction in taxon constraints")
        }

        // Solve
        if (solver.isSatisfiable()) {
            return new ConsistencyResult(consistent: true)
        } else {
            // Extract minimal unsatisfiable core
            try {
                int[] unsatCore = solver.minimalExplanation()
                List<String> coreExplanation = []
                for (int idx = 0; idx < unsatCore.length; idx++) {
                    int var = Math.abs(unsatCore[idx])
                    String taxon = varToTaxon[var]
                    if (taxon) {
                        coreExplanation << (clauseToSource[var] ?: "taxon constraint on ${taxon}")
                    }
                }

                return buildUnsatResult(requirements, coreExplanation.join('; '))
            } catch (Exception e) {
                log.warn("Could not extract UNSAT core: ${e.message}")
                return buildUnsatResult(requirements, "Taxon constraints are unsatisfiable")
            }
        }
    }

    private ConsistencyResult buildUnsatResult(TaxonRequirements requirements, String justification) {
        List<ConsistencyViolation> violations = []

        // Build violation from the conflicting GO terms
        Set<String> allGoTerms = []
        requirements.positiveSource.values().each { allGoTerms.addAll(it) }
        requirements.negativeSource.values().each { allGoTerms.addAll(it) }

        violations << new ConsistencyViolation(
            type: ViolationType.TAXON_CONFLICT,
            severity: Severity.ERROR,
            description: "Taxon constraints from GO annotations are contradictory",
            involvedGoTerms: allGoTerms.toList(),
            suggestedAction: "Remove GO annotations that violate taxon constraints for this organism",
            justification: justification
        )

        new ConsistencyResult(consistent: false, violations: violations)
    }
}

/**
 * Result of a consistency check.
 */
class ConsistencyResult {
    boolean consistent
    List<ConsistencyViolation> violations = []
}
