package gspa.metrics

import gspa.model.ConsistencyViolation
import gspa.model.Genome
import gspa.ontology.GoOntology
import gspa.ontology.SatConsistencyChecker
import gspa.ontology.TaxonConstraints
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Evaluates consistency of genome-scale function annotations.
 * Paper Definition 3: checks that annotations satisfy taxon constraints.
 *
 * Uses SAT4J solver to encode taxon constraints as a satisfiability problem.
 */
class Consistency {

    private static final Logger log = LoggerFactory.getLogger(Consistency)

    SatConsistencyChecker satChecker

    Consistency(SatConsistencyChecker satChecker) {
        this.satChecker = satChecker
    }

    /**
     * Evaluate consistency of a genome's annotations.
     */
    ConsistencyResult evaluate(Genome genome) {
        Set<String> allGoTerms = genome.allGoTerms()

        if (allGoTerms.isEmpty()) {
            return new ConsistencyResult(consistent: true)
        }

        def satResult = satChecker.check(allGoTerms)

        new ConsistencyResult(
            consistent: satResult.consistent,
            violations: satResult.violations,
            totalAnnotatedTerms: allGoTerms.size(),
            constrainedTermCount: satChecker.taxonConstraints.constrainedTermCount()
        )
    }

    /**
     * Evaluate consistency with per-protein violation attribution.
     * For each violation, identifies which proteins contributed the conflicting annotations.
     */
    ConsistencyResult evaluateWithAttribution(Genome genome) {
        def result = evaluate(genome)

        if (!result.consistent) {
            // For each violation, find which proteins have the conflicting GO terms
            result.violations.each { violation ->
                Set<String> conflictingTerms = violation.involvedGoTerms as Set
                genome.proteins.each { protein ->
                    def proteinTerms = protein.annotations.goTermIds()
                    if (proteinTerms.intersect(conflictingTerms)) {
                        violation.involvedProteins << protein.id
                    }
                }
            }
        }

        result
    }
}

class ConsistencyResult {
    boolean consistent
    List<ConsistencyViolation> violations = []
    int totalAnnotatedTerms = 0
    int constrainedTermCount = 0
}
