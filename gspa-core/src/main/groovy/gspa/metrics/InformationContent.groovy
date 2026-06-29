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

    /** Lazily-built intrinsic IC (from the GO DAG), used when no external corpus
     * frequencies have been loaded — so mean IC is meaningful given just the GO. */
    private Map<String, Double> intrinsicIcCache = null

    InformationContent(GoOntology goOntology) {
        this.goOntology = goOntology
    }

    /**
     * Intrinsic information content from the ontology structure alone
     * (Seco et al. 2004): a term subsuming many descendants is general (low IC),
     * a leaf is specific (high IC). IC(t) = -log2((|descendants(t)| + 1) / N).
     * Computed once and cached.
     */
    private Map<String, Double> intrinsicIc() {
        if (intrinsicIcCache != null) return intrinsicIcCache
        Map<String, Integer> descendantCount = [:].withDefault { 0 }
        Set<String> allTerms = goOntology.getAllGoTerms()
        int n = Math.max(allTerms.size(), 1)
        allTerms.each { t ->
            descendantCount[t]                                  // ensure leaves register (0)
            goOntology.getAncestors(t).each { a -> descendantCount[a] = descendantCount[a] + 1 }
        }
        Map<String, Double> ic = [:]
        double log2 = Math.log(2.0d)
        descendantCount.each { term, d ->
            ic[term] = -Math.log((d + 1) / (double) n) / log2
        }
        intrinsicIcCache = ic
        log.info("Computed intrinsic IC for ${ic.size()} GO terms (N=${n})")
        ic
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
        if (totalAnnotations > 0) {
            int freq = termFrequencies[goTerm] ?: 0
            if (freq == 0) return 0.0
            return -Math.log(freq / (double) totalAnnotations) / Math.log(2)
        }
        // No reference corpus loaded: fall back to intrinsic (DAG-structure) IC.
        goOntology == null ? 0.0 : (intrinsicIc()[goTerm] ?: 0.0)
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
