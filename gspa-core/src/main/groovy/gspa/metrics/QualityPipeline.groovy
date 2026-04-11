package gspa.metrics

import gspa.config.EssentialFunctions
import gspa.config.GspaConfig
import gspa.model.Genome
import gspa.model.OrganismDomain
import gspa.model.QualityReport
import gspa.ontology.*
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Full quality evaluation pipeline.
 * Loads all required resources (GO ontology, taxon constraints, pathways, complex terms)
 * and orchestrates completeness, coherence, and consistency evaluation.
 *
 * This is the main entry point for quality assessment.
 */
class QualityPipeline {

    private static final Logger log = LoggerFactory.getLogger(QualityPipeline)

    // Resources
    private GoOntology goOntology
    private GoReasoner goReasoner
    private TaxonConstraints taxonConstraints
    private SatConsistencyChecker satChecker
    private PathwayDatabase pathwayDatabase
    private EssentialFunctions essentialFunctions
    private File complexTermsFile

    // Evaluators
    private Completeness completenessEvaluator
    private Coherence coherenceEvaluator
    private Consistency consistencyEvaluator
    private InformationContent icEvaluator

    // Config
    private GspaConfig config
    private boolean initialized = false
    private File reasonerCacheDir

    QualityPipeline() {
        this.config = new GspaConfig()
    }

    QualityPipeline(GspaConfig config) {
        this.config = config
    }

    // --- Builder-style configuration ---

    QualityPipeline goOwlFile(String path) {
        goOntology = new GoOntology()
        goOntology.loadOwl(path)
        this
    }

    QualityPipeline reasonerCacheDir(File dir) {
        this.reasonerCacheDir = dir
        this
    }

    QualityPipeline goOwlFile(File file) {
        goOntology = new GoOntology()
        goOntology.loadOwl(file)
        this
    }

    QualityPipeline goOntology(GoOntology ontology) {
        this.goOntology = ontology
        this
    }

    QualityPipeline taxonomyFile(File file) {
        ensureTaxonConstraintsInitialized()
        satChecker.loadTaxonomyHierarchy(file)
        this
    }

    QualityPipeline taxonomyHierarchy(Map<String, String> parentMap) {
        ensureTaxonConstraintsInitialized()
        satChecker.loadTaxonomyHierarchy(parentMap)
        this
    }

    QualityPipeline taxonConstraintsFromOntology() {
        if (goOntology == null) throw new IllegalStateException("GO ontology must be loaded first")
        taxonConstraints = new TaxonConstraints()
        taxonConstraints.loadFromGoOntology(goOntology)
        satChecker = new SatConsistencyChecker(taxonConstraints)
        this
    }

    QualityPipeline taxonConstraintsFile(File file) {
        taxonConstraints = new TaxonConstraints()
        taxonConstraints.loadFromTsv(file)
        satChecker = new SatConsistencyChecker(taxonConstraints)
        this
    }

    QualityPipeline ec2goFile(File file) {
        def ec2go = PathwayLoader.loadEc2Go(file)
        pathwayDatabase = PathwayLoader.createFromEc2Go(ec2go)
        this
    }

    QualityPipeline pathwayFile(File pathwayFile, File ec2goFile = null) {
        def ec2go = ec2goFile ? PathwayLoader.loadEc2Go(ec2goFile) : [:]
        pathwayDatabase = PathwayLoader.loadPathways(pathwayFile, ec2go)
        this
    }

    QualityPipeline pathwayDatabase(PathwayDatabase db) {
        this.pathwayDatabase = db
        this
    }

    QualityPipeline complexTermsFile(File file) {
        this.complexTermsFile = file
        this
    }

    QualityPipeline essentialFunctions(EssentialFunctions ef) {
        this.essentialFunctions = ef
        this
    }

    QualityPipeline essentialFunctionsForDomain(OrganismDomain domain) {
        this.essentialFunctions = EssentialFunctions.getDefault(domain)
        this
    }

    // --- Initialization ---

    /**
     * Initialize all evaluators. Must be called after configuration.
     * Will initialize ELK reasoner for coherence evaluation.
     */
    QualityPipeline initialize() {
        if (goOntology == null) {
            throw new IllegalStateException("GO ontology must be loaded before initialization")
        }

        // Initialize ELK reasoner
        goReasoner = new GoReasoner(goOntology)
        if (reasonerCacheDir != null) {
            goReasoner.cacheDir = reasonerCacheDir
        }
        goReasoner.initialize()

        // Default essential functions if not set
        if (essentialFunctions == null) {
            def domain = config.resolveOrganismDomain()
            essentialFunctions = EssentialFunctions.getDefault(domain)

            // Apply runtime modifications from config
            if (config.quality.completeness.addTerms || config.quality.completeness.removeTerms) {
                essentialFunctions = essentialFunctions.withModifications(
                    config.quality.completeness.addTerms,
                    config.quality.completeness.removeTerms
                )
            }
        }

        // Default taxon constraints if not set
        if (taxonConstraints == null) {
            taxonConstraints = new TaxonConstraints()
            if (config.quality.consistency.taxonConstraints) {
                taxonConstraints.loadFromGoOntology(goOntology)
            }
        }
        if (satChecker == null) {
            satChecker = new SatConsistencyChecker(taxonConstraints)
        }

        // Build evaluators
        completenessEvaluator = new Completeness(goOntology, essentialFunctions)
        coherenceEvaluator = new Coherence(goOntology, goReasoner)
        if (pathwayDatabase != null) {
            coherenceEvaluator.pathwayDatabase = pathwayDatabase
        }
        if (complexTermsFile?.exists()) {
            coherenceEvaluator.loadComplexTerms(complexTermsFile)
        }
        consistencyEvaluator = new Consistency(satChecker)
        icEvaluator = new InformationContent(goOntology)

        initialized = true
        log.info("Quality pipeline initialized")
        this
    }

    /**
     * Initialize without ELK (skip process coherence).
     * Useful for quick evaluation when ELK initialization is too slow.
     */
    QualityPipeline initializeLite() {
        if (goOntology == null) {
            throw new IllegalStateException("GO ontology must be loaded before initialization")
        }

        if (essentialFunctions == null) {
            essentialFunctions = EssentialFunctions.getDefault(config.resolveOrganismDomain())
        }
        if (taxonConstraints == null) {
            taxonConstraints = new TaxonConstraints()
        }
        if (satChecker == null) {
            satChecker = new SatConsistencyChecker(taxonConstraints)
        }

        completenessEvaluator = new Completeness(goOntology, essentialFunctions)
        consistencyEvaluator = new Consistency(satChecker)
        icEvaluator = new InformationContent(goOntology)
        // coherenceEvaluator stays null — process coherence will be skipped

        initialized = true
        log.info("Quality pipeline initialized (lite mode, no ELK)")
        this
    }

    // --- Evaluation ---

    /**
     * Evaluate all quality metrics for a genome.
     */
    QualityReport evaluate(Genome genome) {
        if (!initialized) throw new IllegalStateException("Pipeline not initialized. Call initialize() first.")

        log.info("Evaluating quality for: ${genome.id} (${genome.proteinCount} proteins, ${genome.allGoTerms().size()} GO terms)")

        def report = new QualityReport(
            genomeId: genome.id,
            assessmentDate: new java.text.SimpleDateFormat('yyyy-MM-dd').format(new Date()),
            totalProteinCount: genome.proteinCount,
            annotatedProteinCount: genome.proteins.count { !it.annotations.isEmpty() }
        )

        // Completeness
        if (config.quality.completeness.profile != 'skip') {
            def cr = genome.mag ?
                completenessEvaluator.evaluateMAG(genome) :
                completenessEvaluator.evaluate(genome)
            report.completeness = cr.score
            report.presentEssentialFunctions = cr.presentFunctions
            report.missingEssentialFunctions = cr.missingFunctions
            log.info("  Completeness: ${String.format('%.1f%%', cr.score * 100)} " +
                "(${cr.presentFunctions.size()}/${cr.totalEssential} essential functions)")
        }

        // Coherence
        if (coherenceEvaluator != null) {
            def cohr = coherenceEvaluator.evaluate(genome)
            report.processCoherence = cohr.processCoherence
            report.pathwayCoherence = cohr.pathwayCoherence
            report.complexCoherence = cohr.complexCoherence
            log.info("  Process coherence: ${String.format('%.1f%%', cohr.processCoherence * 100)} " +
                "(${cohr.processSatisfied}/${cohr.processTriggered} dependencies)")
            log.info("  Pathway coherence: ${String.format('%.1f%%', cohr.pathwayCoherence * 100)}")
            log.info("  Complex coherence: ${String.format('%.1f%%', cohr.complexCoherence * 100)}")
        } else {
            report.processCoherence = -1.0
            report.pathwayCoherence = -1.0
            report.complexCoherence = -1.0
        }

        // Consistency
        if (config.quality.consistency.taxonConstraints) {
            def consr = consistencyEvaluator.evaluateWithAttribution(genome)
            report.consistent = consr.consistent
            report.violations = consr.violations
            log.info("  Consistency: ${consr.consistent ? 'PASS' : 'FAIL'} " +
                "(${consr.violations.size()} violations)")
        }

        // Information content
        def icr = icEvaluator.evaluate(genome)
        report.meanIC = icr.meanIC
        log.info("  Mean IC: ${String.format('%.2f', icr.meanIC)}")

        // Per-source counts
        genome.proteins.each { protein ->
            protein.annotations.sources().each { source ->
                report.annotationCountBySource[source] =
                    (report.annotationCountBySource[source] ?: 0) +
                    protein.annotations.bySource(source).size()
            }
        }

        log.info("  Composite score: ${String.format('%.3f', report.compositeScore)}")

        genome.qualityReport = report
        report
    }

    /**
     * Evaluate multiple genomes (e.g., MAG bins).
     */
    List<QualityReport> evaluateAll(List<Genome> genomes) {
        genomes.collect { evaluate(it) }
    }

    /**
     * Clean up resources (dispose ELK reasoner).
     */
    void dispose() {
        goReasoner?.dispose()
    }

    private void ensureTaxonConstraintsInitialized() {
        if (taxonConstraints == null) {
            taxonConstraints = new TaxonConstraints()
        }
        if (satChecker == null) {
            satChecker = new SatConsistencyChecker(taxonConstraints)
        }
    }
}
