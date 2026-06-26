package gspa.metrics

import gspa.config.EssentialFunctions
import gspa.config.GspaConfig
import gspa.model.Contig
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
            loadTaxonConstraints()
        }
        if (satChecker == null) {
            satChecker = new SatConsistencyChecker(taxonConstraints)
            applyTaxonomyHierarchy()
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
            loadTaxonConstraints()
        }
        if (satChecker == null) {
            satChecker = new SatConsistencyChecker(taxonConstraints)
            applyTaxonomyHierarchy()
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
            // Per-pair / per-pathway detail for the GAEF report.
            report.incoherentProcessPairs = cohr.processUnsatisfiedPairs ?: []
            report.incoherentPathways = (cohr.pathwayPerPathway ?: []).collect { p ->
                [
                    id           : p.pathwayId,
                    name         : p.pathwayName,
                    pathway_term : p.pathwayGoTerm,
                    completeness : p.completeness,
                    required     : p.requiredCount,
                    present      : p.presentTerms.toList(),
                    missing      : p.missingTerms.toList(),
                ] as Map<String, Object>
            }
            log.info("  Process coherence: ${String.format('%.1f%%', cohr.processCoherence * 100)} " +
                "(${cohr.processSatisfied}/${cohr.processTriggered} dependencies; " +
                "${report.incoherentProcessPairs.size()} unsatisfied pair(s))")
            log.info("  Pathway coherence: ${String.format('%.1f%%', cohr.pathwayCoherence * 100)} " +
                "(${report.incoherentPathways.count { (it.completeness as Double) < 1.0 }}/${report.incoherentPathways.size()} triggered pathways incoherent)")
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

        // Populate the GO ID → name lookup the writer uses to render
        // human-readable GAEF detail. Covers every term referenced by the
        // report (essentials, incoherent process pairs, incoherent pathways,
        // consistency violations).
        if (goOntology != null) {
            Set<String> ids = new LinkedHashSet<>()
            ids.addAll(report.presentEssentialFunctions)
            ids.addAll(report.missingEssentialFunctions)
            report.incoherentProcessPairs.each { p ->
                if (p.key) ids << p.key
                if (p.value) ids << p.value
            }
            report.incoherentPathways.each { pw ->
                if (pw.pathway_term) ids << (pw.pathway_term as String)
                (pw.present as List)?.each { ids << (it as String) }
                (pw.missing as List)?.each { ids << (it as String) }
            }
            report.violations.each { v -> v.involvedGoTerms?.each { ids << it } }
            ids.each { id ->
                try {
                    String label = goOntology.getLabel(id)
                    if (label) report.goLabels[id] = label
                } catch (ignored) { /* unknown id */ }
            }
        }

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
     * Evaluate genome-scale metrics <b>per contig</b> rather than pooled
     * across the whole assembly. Each contig is wrapped as a standalone
     * single-contig {@link Genome} (id = {@code <genomeId>:<contigId>}) so
     * completeness, coherence, consistency and IC are computed over only that
     * contig's proteins. This is the right unit for a metagenome or a
     * multi-replicon assembly, where pooling would mix organisms.
     *
     * Contigs with no annotated proteins are skipped. The shared Protein
     * objects are referenced (not reparented), so the source genome is left
     * untouched.
     */
    List<QualityReport> evaluatePerContig(Genome genome) {
        genome.contigs.findAll { it.proteins }.collect { contig ->
            def sub = new Genome(id: "${genome.id}:${contig.id}", mag: genome.mag,
                taxonomy: genome.taxonomy)
            def wrap = new Contig(id: contig.id, sequence: contig.sequence)
            wrap.genome = sub
            wrap.proteins.addAll(contig.proteins)   // reference, do not reparent
            sub.contigs << wrap
            evaluate(sub)
        }
    }

    /**
     * Run all enabled enforcement passes on the genome's predicted GO
     * annotations, sharing one provenance {@link EnforcementReport}. Order:
     * consistency (demote impossible terms) &rarr; coherence (fix complex
     * singletons, promote missing has_part partners) &rarr; completeness
     * (promote missing essentials). Promotions are consistency-gated. Returns
     * the report (empty when nothing is enabled).
     */
    EnforcementReport enforceAll(Genome genome) {
        if (!initialized) throw new IllegalStateException("Pipeline not initialized. Call initialize() first.")
        def report = new EnforcementReport(provenance: config.quality.provenance)
        def cc = config.quality.consistency
        def cmp = config.quality.completeness
        def coh = config.quality.coherence

        if (cc.enforce) {
            def r = new ConsistencyEnforcer(goOntology: goOntology, satChecker: satChecker,
                mode: cc.enforceMode, downrankFactor: cc.downrankFactor, report: report).enforce(genome)
            log.info("  consistency: ${r}")
        }

        if (coh.enforce) {
            List<Map.Entry<String, String>> pairs = []
            boolean canProcess = coh.enforceProcess && coherenceEvaluator != null && goReasoner != null
            if (canProcess) {
                def pr = coherenceEvaluator.evaluateProcessCoherence(
                    goOntology.propagateAnnotations(genome.allGoTerms()))
                pairs = (pr.unsatisfied ?: []) as List<Map.Entry<String, String>>
            }
            def r = new CoherenceEnforcer(goOntology: goOntology, satChecker: satChecker, report: report,
                complexClassification: coherenceEvaluator?.complexClassification ?: [:],
                promotePartner: coh.promotePartner, enforceProcess: canProcess,
                promoteEvidence: coh.promoteEvidence, promoteMinScore: coh.promoteMinScore).enforce(genome, pairs)
            log.info("  coherence: ${r}")
        }

        if (cmp.enforce) {
            def compRes = genome.mag ?
                completenessEvaluator.evaluateMAG(genome) : completenessEvaluator.evaluate(genome)
            def r = new CompletenessEnforcer(goOntology: goOntology, satChecker: satChecker, report: report,
                promoteEvidence: cmp.promoteEvidence, promoteMinScore: cmp.promoteMinScore)
                .enforce(genome, compRes.missingFunctions)
            log.info("  completeness: ${r}")
        }
        report
    }

    /** True if any enforcement pass is enabled in config. */
    boolean anyEnforcementEnabled() {
        config.quality.consistency.enforce || config.quality.completeness.enforce ||
            config.quality.coherence.enforce
    }

    /**
     * Clean up resources (dispose ELK reasoner).
     */
    void dispose() {
        goReasoner?.dispose()
    }

    /** Bundled (vendored) taxon-constraint resources — see resources/taxon-constraints/PROVENANCE.md. */
    private static final String BUNDLED_CONSTRAINTS = '/taxon-constraints/go-taxon-constraints.tsv'
    private static final String BUNDLED_HIERARCHY = '/taxon-constraints/ncbi-taxon-hierarchy.tsv'

    /**
     * Populate {@link #taxonConstraints} for consistency checking / enforcement.
     * Priority: explicit constraints file (.obo / .tsv) &rarr; the bundled GO
     * taxon constraints (Asaad et al.) &rarr; only_in/never_in axioms from the
     * loaded GO ontology (go-basic carries none, so the bundled file is the
     * default that makes consistency actually act).
     */
    private void loadTaxonConstraints() {
        boolean want = config.quality.consistency.taxonConstraints || config.quality.consistency.enforce
        if (!want) return
        String cf = config.quality.consistency.constraintsFile
        if (cf) {
            File f = new File(cf)
            if (!f.exists()) {
                log.warn("Taxon constraints file not found: ${cf}; falling back to bundled constraints")
            } else if (cf.toLowerCase().endsWith('.obo')) {
                taxonConstraints.loadFromObo(f); return
            } else {
                taxonConstraints.loadFromTsv(f); return
            }
        }
        File bundled = resourceToTempFile(BUNDLED_CONSTRAINTS, 'go-taxon-constraints', '.tsv')
        if (bundled) {
            taxonConstraints.loadFromTsv(bundled)
        } else if (goOntology != null) {
            taxonConstraints.loadFromGoOntology(goOntology)
        }
    }

    /**
     * Load the taxon hierarchy + explicit disjointness into the SAT checker
     * (configured file, else the bundled NCBI disjointness backbone), and
     * assert the organism's taxon if one was configured.
     */
    private void applyTaxonomyHierarchy() {
        String tf = config.quality.consistency.taxonomyFile
        File hierarchy = tf ? new File(tf) : resourceToTempFile(BUNDLED_HIERARCHY, 'ncbi-taxon-hierarchy', '.tsv')
        if (hierarchy?.exists()) {
            satChecker.loadTaxonomyTsv(hierarchy)
        } else if (tf) {
            log.warn("Taxonomy hierarchy file not found: ${tf}")
        }
        String organism = resolveOrganismTaxon()
        if (organism) {
            satChecker.organismTaxon = organism
            log.info("Asserting organism taxon for consistency: ${organism}")
        }
    }

    /** Resolve the configured organism taxon to an {@code NCBITaxon_<id>} string. */
    private String resolveOrganismTaxon() {
        String raw = config.quality.consistency.organismTaxon
        if (!raw) return null
        raw = raw.trim()
        switch (raw.toLowerCase()) {
            case ['bacteria', 'bacterium', 'bacterial']: return 'NCBITaxon_2'
            case ['archaea', 'archaeal']:                return 'NCBITaxon_2157'
            case ['eukaryote', 'eukaryota', 'eukaryotic']: return 'NCBITaxon_2759'
            case ['virus', 'viruses', 'viral']:          return 'NCBITaxon_10239'
        }
        if (raw.startsWith('NCBITaxon')) return raw.replace(':', '_')
        if (raw ==~ /\d+/) return "NCBITaxon_${raw}"
        log.warn("Unrecognized organism taxon '${raw}'; ignoring")
        null
    }

    /** Copy a classpath resource to a temp file so File-based loaders can read it. */
    private File resourceToTempFile(String resourcePath, String prefix, String suffix) {
        def stream = getClass().getResourceAsStream(resourcePath)
        if (stream == null) {
            log.warn("Bundled resource not found on classpath: ${resourcePath}")
            return null
        }
        try {
            File tmp = File.createTempFile(prefix, suffix)
            tmp.deleteOnExit()
            tmp.withOutputStream { out -> stream.transferTo(out) }
            return tmp
        } finally {
            stream.close()
        }
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
