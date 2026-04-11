package gspa.metrics

import gspa.config.EssentialFunctions
import gspa.model.Genome
import gspa.ontology.GoOntology
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Evaluates completeness of genome-scale function annotations.
 * Paper Definition 1: fraction of essential functions present in the annotation set.
 *
 * For each essential GO term F, checks whether any protein in the genome
 * has an annotation that is equal to or a subclass of F (using GO hierarchy).
 */
class Completeness {

    private static final Logger log = LoggerFactory.getLogger(Completeness)

    GoOntology goOntology
    EssentialFunctions essentialFunctions

    Completeness(GoOntology goOntology, EssentialFunctions essentialFunctions) {
        this.goOntology = goOntology
        this.essentialFunctions = essentialFunctions
    }

    /**
     * Evaluate completeness of a genome's annotations.
     *
     * @param genome The genome to evaluate
     * @return CompletenessResult with score and details
     */
    CompletenessResult evaluate(Genome genome) {
        Set<String> annotatedTerms = genome.allGoTerms()

        // Propagate annotations upward through the GO hierarchy
        Set<String> propagated = goOntology.propagateAnnotations(annotatedTerms)

        Set<String> presentFunctions = []
        Set<String> missingFunctions = []

        essentialFunctions.goTerms.each { essentialTerm ->
            if (propagated.contains(essentialTerm)) {
                presentFunctions << essentialTerm
            } else {
                missingFunctions << essentialTerm
            }
        }

        double score = essentialFunctions.goTerms.isEmpty() ? 1.0 :
            presentFunctions.size() / (double) essentialFunctions.goTerms.size()

        // Per-category breakdown
        Map<String, Double> categoryScores = [:]
        essentialFunctions.categories.each { category ->
            def categoryTerms = essentialFunctions.getTermsByCategory(category)
            int present = categoryTerms.count { propagated.contains(it) }
            categoryScores[category] = categoryTerms.isEmpty() ? 1.0 :
                present / (double) categoryTerms.size()
        }

        new CompletenessResult(
            score: score,
            presentFunctions: presentFunctions,
            missingFunctions: missingFunctions,
            categoryScores: categoryScores,
            totalEssential: essentialFunctions.goTerms.size(),
            profileName: essentialFunctions.profileName
        )
    }

    /**
     * Evaluate completeness with MAG adjustment.
     * Adjusts expected completeness based on CheckM2 estimates.
     */
    CompletenessResult evaluateMAG(Genome genome) {
        def result = evaluate(genome)

        if (genome.assemblyInfo?.completeness != null) {
            double checkMCompleteness = genome.assemblyInfo.completeness / 100.0
            result.adjustedScore = result.score / checkMCompleteness
            result.adjustedScore = Math.min(result.adjustedScore, 1.0) // cap at 1.0
            result.magCompleteness = genome.assemblyInfo.completeness
        }

        result
    }
}

class CompletenessResult {
    double score
    Set<String> presentFunctions = []
    Set<String> missingFunctions = []
    Map<String, Double> categoryScores = [:]
    int totalEssential
    String profileName

    // MAG-specific
    Double adjustedScore
    Double magCompleteness
}
