package gspa.metrics

import gspa.config.EssentialFunctions
import gspa.config.GspaConfig
import gspa.model.Genome
import gspa.model.QualityReport
import gspa.ontology.GoOntology
import gspa.ontology.GoReasoner
import gspa.ontology.SatConsistencyChecker
import gspa.ontology.TaxonConstraints
import gspa.ontology.PathwayDatabase
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Orchestrates all quality metrics evaluation for a genome.
 * Combines completeness, coherence, and consistency into a single QualityReport.
 */
class QualityScorer {

    private static final Logger log = LoggerFactory.getLogger(QualityScorer)

    GoOntology goOntology
    GoReasoner goReasoner
    EssentialFunctions essentialFunctions
    TaxonConstraints taxonConstraints
    SatConsistencyChecker satChecker
    PathwayDatabase pathwayDatabase
    File complexTermsFile

    private Completeness completenessChecker
    private Coherence coherenceChecker
    private Consistency consistencyChecker
    private InformationContent icCalculator

    /**
     * Initialize all metric evaluators.
     * Must be called before evaluate().
     */
    void initialize() {
        log.info("Initializing quality scorer...")

        completenessChecker = new Completeness(goOntology, essentialFunctions)

        coherenceChecker = new Coherence(goOntology, goReasoner)
        if (pathwayDatabase != null) {
            coherenceChecker.pathwayDatabase = pathwayDatabase
        }
        if (complexTermsFile?.exists()) {
            coherenceChecker.loadComplexTerms(complexTermsFile)
        }

        consistencyChecker = new Consistency(satChecker)
        icCalculator = new InformationContent(goOntology)

        log.info("Quality scorer initialized")
    }

    /**
     * Evaluate all quality metrics for a genome.
     */
    QualityReport evaluate(Genome genome) {
        log.info("Evaluating quality for genome: ${genome.id}")

        // Completeness
        def completenessResult = genome.mag ?
            completenessChecker.evaluateMAG(genome) :
            completenessChecker.evaluate(genome)

        // Coherence
        def coherenceResult = coherenceChecker.evaluate(genome)

        // Consistency
        def consistencyResult = consistencyChecker.evaluateWithAttribution(genome)

        // Information content
        def icResult = icCalculator.evaluate(genome)

        // Build report
        def report = new QualityReport(
            genomeId: genome.id,
            assessmentDate: new Date().format('yyyy-MM-dd'),
            completeness: completenessResult.score,
            presentEssentialFunctions: completenessResult.presentFunctions,
            missingEssentialFunctions: completenessResult.missingFunctions,
            processCoherence: coherenceResult.processCoherence,
            pathwayCoherence: coherenceResult.pathwayCoherence,
            complexCoherence: coherenceResult.complexCoherence,
            consistent: consistencyResult.consistent,
            violations: consistencyResult.violations,
            meanIC: icResult.meanIC,
            annotatedProteinCount: genome.proteins.count { !it.annotations.isEmpty() },
            totalProteinCount: genome.proteinCount
        )

        // Count annotations by source
        genome.proteins.each { protein ->
            protein.annotations.sources().each { source ->
                report.annotationCountBySource[source] =
                    (report.annotationCountBySource[source] ?: 0) +
                    protein.annotations.bySource(source).size()
            }
        }

        genome.qualityReport = report
        log.info("Quality evaluation complete for ${genome.id}: " +
            "completeness=${String.format('%.1f%%', report.completeness * 100)}, " +
            "coherence=(proc=${String.format('%.1f%%', report.processCoherence * 100)}, " +
            "path=${String.format('%.1f%%', report.pathwayCoherence * 100)}, " +
            "cplx=${String.format('%.1f%%', report.complexCoherence * 100)}), " +
            "consistent=${report.consistent}")

        report
    }
}
