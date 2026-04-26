package gspa.predictor

import gspa.config.GspaConfig
import gspa.integration.ClaimExtractor
import gspa.integration.EvidenceClaim
import gspa.integration.EvidenceCombiner
import gspa.integration.IntegratedAnnotationSet
import gspa.integration.IntegrationState
import gspa.integration.IntegrationWriter
import gspa.integration.IterativeRefiner
import gspa.integration.suggester.DarkMatterSuggester
import gspa.integration.suggester.DarkMatterWriter
import gspa.io.FastaReader
import gspa.io.GffReader
import gspa.io.GffWriter
import gspa.model.*
import gspa.predictor.cluster.ClusterAnnotationPropagator
import gspa.predictor.cluster.IntragenomeClusterer
import gspa.predictor.context.OperonPredictor
import gspa.predictor.genecalling.GeneCaller
import gspa.predictor.genecalling.ProdigalCaller
import gspa.predictor.genecalling.PyrodigalCaller
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Main annotation pipeline.
 * Orchestrates: gene calling -> predictor execution -> annotation merging -> output.
 *
 * Steps:
 * 1. Load input genome (FASTA or GFF3)
 * 2. Call genes (if needed)
 * 3. Run configured predictors in parallel
 * 4. Merge annotations from all sources
 * 5. Run operon/context predictors
 * 6. Output annotated genome (GFF3, GenBank, etc.)
 */
class AnnotationPipeline {

    private static final Logger log = LoggerFactory.getLogger(AnnotationPipeline)

    GspaConfig config
    PredictorRegistry registry = new PredictorRegistry()
    GeneCaller geneCaller
    List<GenomePredictor> genomePredictors = []

    /** Run protein-level predictors in parallel threads */
    boolean parallel = true

    AnnotationPipeline(GspaConfig config) {
        this.config = config
    }

    /**
     * Auto-configure the pipeline based on config settings.
     * Discovers available tools and registers them.
     */
    void configure() {
        log.info("Configuring annotation pipeline...")

        // Configure gene caller
        configureGeneCaller()

        // Register available predictors
        registerPredictors()

        log.info("Pipeline configured: geneCaller=${geneCaller?.name ?: 'none'}, " +
            "${registry.getAll().size()} predictors registered")
    }

    private void configureGeneCaller() {
        def domain = config.resolveOrganismDomain()

        boolean isMeta = config.input.type in ['mag', 'metagenome']

        // For eukaryotes, Augustus would be preferred but we fall back to Prodigal
        // (works for intron-poor eukaryotes like microsporidia, fungi)
        // Try pyrodigal first (faster), fall back to prodigal
        def pyrodigal = new PyrodigalCaller(metaMode: isMeta)
        if (pyrodigal.isAvailable()) {
            geneCaller = pyrodigal
        } else {
            def prodigal = new ProdigalCaller(metaMode: isMeta)
            if (prodigal.isAvailable()) {
                geneCaller = prodigal
            }
        }

        if (domain.isEukaryote() && geneCaller != null) {
            log.warn("Using ${geneCaller.name} for eukaryote gene calling — " +
                "results may miss intron-containing genes. Install Augustus for better eukaryote support.")
        }

        if (geneCaller != null) {
            log.info("Gene caller: ${geneCaller.name}")
        } else {
            log.warn("No gene caller available; input must include protein sequences")
        }
    }

    private void registerPredictors() {
        // Register all known predictor types and check availability
        def allPredictors = createAllPredictors()

        allPredictors.each { predictor ->
            String name = predictor.name
            boolean disabled = config.predictors.disable.contains(name)
            if (disabled) {
                log.debug("Predictor ${name} explicitly disabled")
                return
            }
            registry.register(predictor)
        }
    }

    private List<Predictor> createAllPredictors() {
        List<Predictor> predictors = []

        // Similarity
        if (config.predictors.similarity.enabled) {
            switch (config.predictors.similarity.tool) {
                case 'diamond':
                    def p = new gspa.predictor.similarity.DiamondPredictor(
                        database: config.predictors.similarity.database,
                        evalue: config.predictors.similarity.evalue,
                        maxTargetSeqs: config.predictors.similarity.maxTargetSeqs,
                    )
                    predictors << p
                    break
                case 'mmseqs2':
                    def p = new gspa.predictor.similarity.MMseqs2Predictor(
                        database: config.predictors.similarity.database,
                        evalue: config.predictors.similarity.evalue,
                    )
                    predictors << p
                    break
            }
        }

        // Structure (homology-transfer by default; centroid mode when configured)
        if (config.predictors.structure.enabled) {
            def structureCfg = config.predictors.structure
            String modeStr = (structureCfg.centroidMode ?: 'none').toUpperCase(Locale.ROOT)
            def mode = gspa.predictor.structure.FoldSeekPredictor.CentroidMode.valueOf(modeStr)
            String dbPath = mode != gspa.predictor.structure.FoldSeekPredictor.CentroidMode.NONE &&
                    structureCfg.centroidDb ?
                    structureCfg.centroidDb :
                    structureCfg.database
            def p = new gspa.predictor.structure.FoldSeekPredictor(
                database: dbPath,
                prostt5Model: structureCfg.prostt5Model,
                centroidMode: mode,
            )
            predictors << p
        }

        // Domains
        if (config.predictors.domains.enabled) {
            switch (config.predictors.domains.tool) {
                case 'interproscan':
                    def p = new gspa.predictor.domain.InterProScanPredictor(
                        applications: config.predictors.domains.applications,
                    )
                    predictors << p
                    break
                case 'hmmer':
                    predictors << new gspa.predictor.domain.HmmerPredictor()
                    break
            }
        }

        // Operon predictor (genome-level)
        if (config.predictors.operons.enabled) {
            def op = new OperonPredictor(
                maxIntergenicDistance: config.predictors.operons.maxIntergenicDistance,
                requireSameStrand: config.predictors.operons.requireSameStrand,
            )
            predictors << op
            genomePredictors << op
        }

        // Pathway
        if (config.predictors.pathway.enabled) {
            def gp = new gspa.predictor.pathway.GapseqPredictor()
            predictors << gp
            genomePredictors << gp
        }

        // Localization
        if (config.predictors.localization.enabled) {
            if (config.predictors.localization.signalP) {
                predictors << new gspa.predictor.localization.SignalPPredictor(
                    organismType: config.resolveOrganismDomain().isProkaryote() ? 'gram-' : 'eukarya'
                )
            }
            if (config.predictors.localization.tmhmm) {
                predictors << new gspa.predictor.localization.DeepTmhmmPredictor()
            }
        }

        // eggNOG-mapper
        if (config.predictors.enable.contains('eggnog-mapper')) {
            predictors << new gspa.predictor.function.EggNogMapperPredictor()
        }

        // Disorder (Metapredict)
        if (config.predictors.disorder.enabled) {
            predictors << new gspa.predictor.disorder.MetapredictPredictor(
                minRegionLen: config.predictors.disorder.minRegionLen,
                minScore: config.predictors.disorder.minScore,
            )
        }

        // Neural sidecar predictors
        def neural = config.predictors.neural
        if (neural.esm2DeepGoPlus.enabled) {
            predictors << new gspa.predictor.neural.DeepGoPlusEsm2Predictor(
                sidecarScript: neural.sidecarScript,
                pythonExecutable: neural.pythonExecutable,
                checkpoint: neural.esm2DeepGoPlus.checkpoint,
                terms: neural.esm2DeepGoPlus.terms,
                esm2Model: neural.esm2DeepGoPlus.model,
                batchSize: neural.esm2DeepGoPlus.batchSize,
                minScore: neural.esm2DeepGoPlus.minScore,
            )
        }
        if (neural.proteinfer.enabled) {
            predictors << new gspa.predictor.neural.ProteInferPredictor(
                sidecarScript: neural.sidecarScript,
                pythonExecutable: neural.pythonExecutable,
                modelDir: neural.proteinfer.modelDir,
                batchSize: neural.proteinfer.batchSize,
                minScore: neural.proteinfer.minScore,
            )
        }
        if (neural.clean.enabled) {
            predictors << new gspa.predictor.neural.CleanPredictor(
                sidecarScript: neural.sidecarScript,
                pythonExecutable: neural.pythonExecutable,
                modelDir: neural.clean.modelDir,
                batchSize: neural.clean.batchSize,
                minScore: neural.clean.minScore,
            )
        }
        if (neural.esm2Centroid.enabled) {
            predictors << new gspa.predictor.neural.Esm2CentroidPredictor(
                sidecarScript: neural.sidecarScript,
                pythonExecutable: neural.pythonExecutable,
                centroidDb: neural.esm2Centroid.db,
                esm2Model: neural.esm2Centroid.model,
                topK: neural.esm2Centroid.topK,
                batchSize: neural.esm2Centroid.batchSize,
                minScore: neural.esm2Centroid.minScore,
            )
        }

        // FOSS region predictors (use run_region_predictors.py sidecar)
        def loc = config.predictors.localization
        def foss = config.predictors.foss
        if (loc.deepSig) {
            predictors << new gspa.predictor.localization.DeepSigPredictor(
                sidecarScript: foss.regionSidecar,
                pythonExecutable: foss.pythonExecutable,
                kingdom: loc.deepSigKingdom,
            )
        }
        if (loc.tmbed) {
            predictors << new gspa.predictor.localization.TmbedPredictor(
                sidecarScript: foss.regionSidecar,
                pythonExecutable: foss.pythonExecutable,
            )
        }
        if (loc.tppred3) {
            predictors << new gspa.predictor.localization.TPpred3Predictor(
                sidecarScript: foss.regionSidecar,
                pythonExecutable: foss.pythonExecutable,
                kingdom: loc.tppred3Kingdom,
            )
        }
        // FOSS term predictors (use run_term_predictors.py sidecar)
        if (loc.psortb) {
            predictors << new gspa.predictor.localization.PSORTbPredictor(
                sidecarScript: foss.termSidecar,
                pythonExecutable: foss.pythonExecutable,
                gram: loc.psortbGram,
            )
        }
        if (foss.deepFri.enabled) {
            predictors << new gspa.predictor.structure.DeepFriPredictor(
                sidecarScript: foss.termSidecar,
                pythonExecutable: foss.pythonExecutable,
                modelDir: foss.deepFri.modelDir,
                structureMode: foss.deepFri.mode,
                minScore: foss.deepFri.minScore,
            )
        }
        if (foss.deepEc.enabled) {
            if (!foss.deepEc.acknowledgeAgpl) {
                log.warn("DeepEC is AGPL-3.0; if running as a service you must " +
                         "publish source. Set predictors.foss.deepEc.acknowledgeAgpl=true to silence this warning.")
            }
            predictors << new gspa.predictor.neural.DeepEcPredictor(
                sidecarScript: foss.termSidecar,
                pythonExecutable: foss.pythonExecutable,
                modelDir: foss.deepEc.modelDir,
                minScore: foss.deepEc.minScore,
            )
        }
        if (foss.deepArg.enabled) {
            predictors << new gspa.predictor.specialized.DeepArgPredictor(
                sidecarScript: foss.termSidecar,
                pythonExecutable: foss.pythonExecutable,
                modelDir: foss.deepArg.modelDir,
                type: foss.deepArg.type,
                minScore: foss.deepArg.minScore,
            )
        }
        // FOSS site predictors (use run_site_predictors.py sidecar)
        if (foss.musiteDeep.enabled) {
            predictors << new gspa.predictor.sites.MusiteDeepPredictor(
                sidecarScript: foss.siteSidecar,
                pythonExecutable: foss.pythonExecutable,
                modelDir: foss.musiteDeep.modelDir,
                residueTypes: foss.musiteDeep.residueTypes,
                minScore: foss.musiteDeep.minScore,
            )
        }
        if (foss.scanNet.enabled) {
            predictors << new gspa.predictor.sites.ScanNetPredictor(
                sidecarScript: foss.siteSidecar,
                pythonExecutable: foss.pythonExecutable,
                modelDir: foss.scanNet.modelDir,
                structureDir: foss.scanNet.structureDir,
                minScore: foss.scanNet.minScore,
            )
        }

        predictors
    }

    /**
     * Run the full annotation pipeline on a genome.
     */
    Genome annotate(File inputFile) {
        log.info("=== GSPA Annotation Pipeline ===")
        log.info("Input: ${inputFile}")

        // 1. Load genome
        Genome genome = loadGenome(inputFile)
        log.info("Loaded: ${genome}")

        // 2. Call genes if needed
        if (genome.proteinCount == 0 && geneCaller != null) {
            log.info("Calling genes with ${geneCaller.name}...")
            geneCaller.callGenes(genome)
            log.info("Called ${genome.proteinCount} genes")
        }

        if (genome.proteinCount == 0) {
            log.warn("No proteins found. Check input format or install a gene caller.")
            return genome
        }

        // 3. Run protein-level predictors (in parallel if multiple available)
        def availablePredictors = registry.getAvailable()
            .findAll { !(it instanceof GenomePredictor) }
        log.info("Running ${availablePredictors.size()} protein-level predictors...")

        // Phase 10: optional intra-genome clustering — predictors see only
        // cluster representatives; annotations are propagated to members
        // afterwards. When disabled (default) clusterSet is null and
        // predictorInput == all proteins.
        def clusterCfg = config.integration?.intragenomeCluster
        def clusterSet = null
        List<Protein> predictorInput = genome.proteins
        if (clusterCfg?.enabled) {
            log.info("Clustering proteome at identity=${clusterCfg.identity}, coverage=${clusterCfg.coverage}...")
            clusterSet = new IntragenomeClusterer().cluster([genome], clusterCfg.identity, clusterCfg.coverage)
            Set<String> repIds = clusterSet.representatives()*.proteinId as Set
            predictorInput = genome.proteins.findAll { it.id in repIds }
            log.info("  clustered ${genome.proteinCount} proteins → ${clusterSet.clusterCount()} reps " +
                "(${String.format(Locale.ROOT, '%.1f', 100.0d * clusterSet.clusterCount() / Math.max(1, genome.proteinCount))}%)")
        }

        if (availablePredictors.size() > 1 && parallel) {
            log.info("  Running predictors in parallel...")
            def results = Collections.synchronizedMap(new LinkedHashMap<String, Map<String, List>>())
            def threads = availablePredictors.collect { predictor ->
                Thread.start {
                    log.info("  [parallel] Running ${predictor.name}...")
                    try {
                        results[predictor.name] = predictor.predictBatch(predictorInput)
                    } catch (Exception e) {
                        log.error("  ${predictor.name} failed: ${e.message}")
                    }
                }
            }
            threads.each { it.join() }
            results.each { name, annotations ->
                applyAnnotations(genome, annotations as Map<String, List>, name)
            }
        } else {
            availablePredictors.each { predictor ->
                log.info("  Running ${predictor.name}...")
                try {
                    def annotations = predictor.predictBatch(predictorInput)
                    applyAnnotations(genome, annotations, predictor.name)
                } catch (Exception e) {
                    log.error("  ${predictor.name} failed: ${e.message}")
                }
            }
        }

        // Phase 10: propagate representative annotations to cluster members.
        if (clusterSet != null) {
            int copied = new ClusterAnnotationPropagator().propagate(clusterSet, [genome])
            log.info("  propagated ${copied} annotations from ${clusterSet.clusterCount()} reps to members")
        }

        // 4. Run genome-level predictors (sequential — they may depend on protein annotations)
        def availableGenomePredictors = genomePredictors.findAll { it.isAvailable() }
        if (availableGenomePredictors) {
            log.info("Running ${availableGenomePredictors.size()} genome-level predictors...")
            availableGenomePredictors.each { predictor ->
                log.info("  Running ${predictor.name}...")
                try {
                    def annotations = predictor.predictGenome(genome)
                    applyAnnotations(genome, annotations, predictor.name)
                } catch (Exception e) {
                    log.error("  ${predictor.name} failed: ${e.message}")
                }
            }
        }

        // 5. Evidence integration (Phase 7+), gated by config.
        if (config.integration?.enabled) {
            runIntegration(genome)
        }

        // Summary
        int totalAnnotations = genome.proteins.sum { it.annotations.size() } ?: 0
        int annotatedProteins = genome.proteins.count { !it.annotations.isEmpty() }
        log.info("Annotation complete: ${annotatedProteins}/${genome.proteinCount} proteins annotated, " +
            "${totalAnnotations} total annotations")

        genome
    }

    /**
     * Run the Phase 7 evidence integration layer. Collects raw claims from
     * every protein's AnnotationSet, runs the iterative refiner, and
     * replaces each protein's annotations with the integrated posteriors.
     * Optionally follows with Phase 8's DarkMatterSuggester.
     */
    private void runIntegration(Genome genome) {
        log.info("=== Evidence integration (Phase 7) ===")
        def extractor = new ClaimExtractor()
        List<EvidenceClaim> claims = extractor.fromGenome(genome)
        log.info("  extracted ${claims.size()} claims from ${genome.proteinCount} proteins")

        if (claims.isEmpty()) {
            log.info("  no claims to integrate; skipping")
            return
        }

        def combiner = new EvidenceCombiner()
        def refiner = new IterativeRefiner(combiner)
        refiner.maxIter = config.integration.maxIter
        refiner.epsilon = config.integration.epsilon
        refiner.damping = config.integration.damping

        def state = new IntegrationState(genome)
        IntegratedAnnotationSet integrated = refiner.refine(claims, state)
        log.info("  ${integrated}")

        IntegrationWriter.applyIntegratedAnnotations(genome, integrated)

        // Phase 8: dark-matter / contextual gap suggester.
        if (config.integration.darkMatter?.enabled) {
            log.info("=== Dark-matter suggester (Phase 8) ===")
            def suggester = new DarkMatterSuggester()
            suggester.bfMin             = config.integration.darkMatter.bfMin
            suggester.gammaInP          = config.integration.darkMatter.gammaInP
            suggester.coverageThreshold = config.integration.darkMatter.coverageThreshold
            suggester.suggest(state, integrated)
            int touched = DarkMatterWriter.apply(genome, integrated, config.integration.darkMatter.mode)
            log.info("  ${integrated.suggestions.size()} suggestions emitted; ${touched} proteins updated")
        }
    }

    /**
     * Apply predicted annotations to genome proteins.
     */
    private void applyAnnotations(Genome genome, Map<String, List<Annotation>> annotations, String source) {
        int count = 0
        annotations.each { proteinId, anns ->
            def protein = genome.findProtein(proteinId)
            if (protein) {
                protein.annotations.addAll(anns)
                count += anns.size()
            }
        }
        log.info("  ${source}: ${count} annotations assigned")
    }

    /**
     * Load genome from various input formats.
     */
    private Genome loadGenome(File inputFile) {
        String name = inputFile.name.toLowerCase()

        if (name.endsWith('.gff3') || name.endsWith('.gff')) {
            return GffReader.readGff3(inputFile)
        }

        // Default: treat as nucleotide FASTA
        def genome = FastaReader.readGenome(inputFile)
        genome.mag = config.input.type in ['mag', 'metagenome']

        // If protein FASTA provided separately
        if (config.input.proteinFasta) {
            def proteins = FastaReader.readProteins(new File(config.input.proteinFasta))
            if (genome.contigs.isEmpty()) {
                genome.addContig(new Contig(id: 'unknown'))
            }
            proteins.each { protein ->
                genome.contigs.first().addProtein(protein)
            }
        }

        genome
    }

    /**
     * Write annotated genome to output files.
     */
    void writeOutput(Genome genome, File outputDir) {
        outputDir.mkdirs()

        if (config.output.formats.contains('gff3')) {
            def gffFile = new File(outputDir, "${genome.id}.gff3")
            GffWriter.writeGff3(genome, gffFile, true)
            log.info("GFF3 written: ${gffFile}")
        }

        if (config.output.formats.contains('tsv')) {
            def tsvFile = new File(outputDir, "${genome.id}_annotations.tsv")
            writeAnnotationsTsv(genome, tsvFile)
            log.info("TSV written: ${tsvFile}")
        }
    }

    private void writeAnnotationsTsv(Genome genome, File output) {
        output.withWriter { writer ->
            writer.writeLine(['protein_id', 'type', 'value', 'score', 'source', 'evidence'].join('\t'))
            genome.proteins.each { protein ->
                protein.annotations.annotations.each { ann ->
                    writer.writeLine([
                        protein.id,
                        ann.type,
                        ann.value,
                        String.format('%.4f', ann.score),
                        ann.source ?: '',
                        ann.evidence ?: '',
                    ].join('\t'))
                }
            }
        }
    }
}
