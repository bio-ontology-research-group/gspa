package gspa.metrics

import gspa.model.Genome
import gspa.ontology.GoOntology
import gspa.ontology.GoReasoner
import gspa.ontology.PathwayDatabase
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Evaluates coherence of genome-scale function annotations.
 * Paper Definition 2: three sub-metrics measuring functional dependencies.
 */
class Coherence {

    private static final Logger log = LoggerFactory.getLogger(Coherence)

    GoOntology goOntology
    GoReasoner goReasoner
    PathwayDatabase pathwayDatabase

    /** Complex terms classification (from complex_terms.tsv) */
    Map<String, String> complexClassification = [:] // GO term -> h/n/a

    Coherence(GoOntology goOntology, GoReasoner goReasoner) {
        this.goOntology = goOntology
        this.goReasoner = goReasoner
    }

    /**
     * Load protein complex classification from TSV.
     * Format: GO_term\tclassification\tdescription
     * classification: h=homomeric, n=heteromeric, a=ambiguous
     */
    void loadComplexTerms(File complexTermsFile) {
        complexTermsFile.eachLine { line ->
            if (line.startsWith('GO_term') || line.startsWith('#')) return
            def fields = line.split('\t')
            if (fields.length >= 2) {
                complexClassification[fields[0].trim()] = fields[1].trim()
            }
        }
        log.info("Loaded ${complexClassification.size()} complex term classifications")
    }

    /**
     * Evaluate all three coherence sub-metrics.
     */
    CoherenceResult evaluate(Genome genome) {
        Set<String> annotatedTerms = genome.allGoTerms()

        def processResult = evaluateProcessCoherence(annotatedTerms)
        def pathwayResult = evaluatePathwayCoherence(annotatedTerms)
        def complexResult = evaluateComplexCoherence(genome)

        new CoherenceResult(
            processCoherence: processResult.score,
            processTriggered: processResult.triggered,
            processSatisfied: processResult.satisfied,
            pathwayCoherence: pathwayResult,
            complexCoherence: complexResult.score,
            complexTotal: complexResult.total,
            complexCoherent: complexResult.coherent
        )
    }

    /**
     * Process coherence: fraction of triggered has_part dependencies that are satisfied.
     */
    ProcessCoherenceResult evaluateProcessCoherence(Set<String> annotatedTerms) {
        def result = goReasoner.checkProcessCoherence(annotatedTerms)
        int triggered = result.triggered as int
        int satisfied = result.satisfied as int

        double score = triggered == 0 ? 1.0 : satisfied / (double) triggered

        new ProcessCoherenceResult(
            score: score,
            triggered: triggered,
            satisfied: satisfied,
            unsatisfied: result.unsatisfied as List
        )
    }

    /**
     * Pathway coherence: average max completeness across triggered pathways.
     */
    double evaluatePathwayCoherence(Set<String> annotatedTerms) {
        if (pathwayDatabase == null) {
            log.warn("No pathway database loaded, skipping pathway coherence")
            return 1.0
        }
        pathwayDatabase.computePathwayCoherence(annotatedTerms)
    }

    /**
     * Complex coherence: fraction of triggered complexes that are coherent.
     * Heteromeric complexes require >= 2 distinct proteins annotated.
     * Homomeric complexes are always coherent if at least one protein is annotated.
     */
    ComplexCoherenceResult evaluateComplexCoherence(Genome genome) {
        if (complexClassification.isEmpty()) {
            return new ComplexCoherenceResult(score: 1.0)
        }

        Map<String, Set<String>> proteinsByComplex = [:].withDefault { [] as Set }

        // For each protein, check if it's annotated with any complex term
        genome.proteins.each { protein ->
            protein.annotations.goTermIds().each { goTerm ->
                if (complexClassification.containsKey(goTerm)) {
                    proteinsByComplex[goTerm] << protein.id
                }
            }
        }

        int total = 0
        int coherent = 0

        proteinsByComplex.each { complexTerm, proteinIds ->
            String classification = complexClassification[complexTerm]
            total++

            switch (classification) {
                case 'h': // homomeric: always coherent if present
                    coherent++
                    break
                case 'n': // heteromeric: need >= 2 distinct proteins
                    if (proteinIds.size() >= 2) coherent++
                    break
                case 'a': // ambiguous: treat as coherent (conservative)
                    coherent++
                    break
            }
        }

        double score = total == 0 ? 1.0 : coherent / (double) total

        new ComplexCoherenceResult(score: score, total: total, coherent: coherent)
    }
}

class CoherenceResult {
    double processCoherence
    int processTriggered
    int processSatisfied
    double pathwayCoherence
    double complexCoherence
    int complexTotal
    int complexCoherent
}

class ProcessCoherenceResult {
    double score
    int triggered
    int satisfied
    List unsatisfied = []
}

class ComplexCoherenceResult {
    double score = 1.0
    int total = 0
    int coherent = 0
}
