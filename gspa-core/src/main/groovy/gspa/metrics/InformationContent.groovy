package gspa.metrics

import gspa.model.Genome
import gspa.ontology.GoOntology
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Computes Information Content (IC) for GO annotations.
 * IC measures how specific a GO term is based on its frequency in a reference corpus.
 *
 * IC(t) = -log2(p(t)) where p(t) is the probability of term t being annotated.
 * Higher IC = more specific term.
 */
class InformationContent {

    private static final Logger log = LoggerFactory.getLogger(InformationContent)

    GoOntology goOntology

    /** Term frequencies from a reference corpus (e.g., UniProt-GOA) */
    Map<String, Integer> termFrequencies = [:]
    int totalAnnotations = 0

    InformationContent(GoOntology goOntology) {
        this.goOntology = goOntology
    }

    /**
     * Load term frequencies from a reference annotation file.
     * Builds a frequency table by propagating each annotation upward.
     */
    void loadFrequencies(Map<String, Set<String>> proteinToGoTerms) {
        termFrequencies.clear()
        totalAnnotations = proteinToGoTerms.size()

        proteinToGoTerms.values().each { terms ->
            // Propagate each term's ancestors
            Set<String> allTerms = goOntology.propagateAnnotations(terms)
            allTerms.each { term ->
                termFrequencies[term] = (termFrequencies[term] ?: 0) + 1
            }
        }

        log.info("Loaded IC frequencies for ${termFrequencies.size()} terms from ${totalAnnotations} proteins")
    }

    /**
     * Get the IC of a specific GO term.
     */
    double getIC(String goTerm) {
        if (totalAnnotations == 0) return 0.0
        int freq = termFrequencies[goTerm] ?: 0
        if (freq == 0) return 0.0
        -Math.log(freq / (double) totalAnnotations) / Math.log(2)
    }

    /**
     * Compute IC-based metrics for a genome's annotations.
     */
    ICResult evaluate(Genome genome) {
        List<Double> ics = []
        Map<String, Double> termICs = [:]

        genome.allGoTerms().each { term ->
            double ic = getIC(term)
            ics << ic
            termICs[term] = ic
        }

        double meanIC = ics.isEmpty() ? 0.0 : ics.sum() / ics.size()
        double maxIC = ics.isEmpty() ? 0.0 : ics.max()

        // IC breadth: average number of proteins per annotated GO class
        Map<String, Integer> proteinCountPerTerm = [:].withDefault { 0 }
        genome.proteins.each { protein ->
            protein.annotations.goTermIds().each { term ->
                proteinCountPerTerm[term]++
            }
        }
        double icBreadth = proteinCountPerTerm.isEmpty() ? 0.0 :
            proteinCountPerTerm.values().sum() / (double) proteinCountPerTerm.size()

        new ICResult(
            meanIC: meanIC,
            maxIC: maxIC,
            icBreadth: icBreadth,
            termICs: termICs,
            annotatedTermCount: ics.size()
        )
    }
}

class ICResult {
    double meanIC
    double maxIC
    double icBreadth
    Map<String, Double> termICs = [:]
    int annotatedTermCount
}
