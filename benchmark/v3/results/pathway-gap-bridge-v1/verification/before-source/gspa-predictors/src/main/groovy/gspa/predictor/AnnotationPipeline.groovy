package gspa.predictor

import gspa.config.GspaConfig
import gspa.integration.ClaimExtractor
import gspa.integration.EvidenceClaim
import gspa.integration.IntegrationRunAudit
import gspa.integration.IntegrationState
import gspa.integration.workflow.IntegrationWorkflow
import gspa.integration.workflow.IntegrationWorkflowConfig
import gspa.io.CdsTranslator
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

    /**
     * Typed completion status of every registered predictor from the last
     * {@link #annotate} call (unavailable / failed / successful-zero / successful).
     * Predictors that were not requested are deliberately absent: those switched off in the
     * config (never created) and those named in {@code predictors.disable} (never registered)
     * are "not requested", not "unavailable". Everything that was requested and registered
     * appears exactly once, plus {@code gene-caller:<name>} when genes were called.
     */
    final List<PredictorRunStatus> runStatuses = Collections.synchronizedList([])

    /** Directory for integration artifacts when the workflow config names no output file (null: none written). */
    File integrationOutputDir

    /** Claims of the last {@link #annotate} call, captured before integration or filtering. */
    final Map<String, Map> pathwaySearchReports = [:]

    List<EvidenceClaim> rawClaims

    /** Outcome of the last {@link #integrate} call (null when integration is off or failed). */
    IntegrationWorkflow.Result lastIntegration

    /** Stage ledger of the last integration attempt, kept even when a requested stage failed. */
    IntegrationRunAudit lastIntegrationAudit

    /**
     * Heavyweight integration references (ontologies, pathways, reaction graph) loaded once by
     * {@link IntegrationWorkflow#loadSharedReferences} and reused across genomes; null loads per run.
     */
    IntegrationState sharedIntegrationReferences

    /** Explicit Pfam-A HMM file for {@code domains.tool: hmmer}; set before {@link #configure}. */
    String hmmerDatabase

    /** Environment variable naming the Pfam-A HMM file for the direct HMMER domain tool. */
    static final String PFAM_HMM_ENV = 'GSPA_PFAM_HMM'

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
        if (config.input.requireExistingGenes || config.input.proteinsOnly || config.input.proteinFasta || config.input.gff3) return
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
                genomePredictors.removeAll { it.name == name }
                return
            }
            registry.register(predictor)
        }
    }

    /**
     * Pfam-A HMM path for the direct HMMER tool. Precedence: the explicit
     * {@link #hmmerDatabase} override, then {@code predictors.domains.database}
     * from the config, then {@link #PFAM_HMM_ENV}.
     */
    private String resolveHmmerDatabase() {
        hmmerDatabase ?: config.predictors.domains.database ?: System.getenv(PFAM_HMM_ENV)
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
                        assignProductNames: config.predictors.similarity.assignProductNames,
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
                        threads: config.predictors.domains.threads,
                        disablePrecalc: config.predictors.domains.disablePrecalc,
                    )
                    predictors << p
                    break
                case 'hmmer':
                    predictors << new gspa.predictor.domain.HmmerPredictor(
                        hmmDatabase: resolveHmmerDatabase(),
                        threads: config.predictors.domains.threads,
                    )
                    break
            }
        }

        // AntiFam — spurious / pseudogene HMM scan (same step Bakta uses)
        if (config.predictors.antifam.enabled) {
            def antifamCfg = config.predictors.antifam
            predictors << new gspa.predictor.domain.AntiFamPredictor(
                hmmDatabase : antifamCfg.database,
                threads     : antifamCfg.threads,
                useCutGa    : antifamCfg.useCutGa,
                domainEvalue: antifamCfg.evalue,
            )
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
            def pathway = config.predictors.pathway
            def gp = new gspa.predictor.pathway.GapseqPredictor(
                tool: pathway.tool, target: pathway.target, dataDir: pathway.dataDir,
                executablePath: pathway.executablePath, taxonomy: pathway.taxonomy,
                threads: pathway.threads,
            )
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
        // Resolve the DeepGO-PlusPlus base-predictor selector (full | light | none)
        // into the two enabled flags before registering them.
        neural.resolveBasePredictor()
        if (neural.deepGoPlusPlus.enabled && neural.deepGoPlusPlusLight.enabled) {
            log.warn('Both DeepGO-PlusPlus (full) and DeepGO-PlusPlus-Light are enabled; ' +
                     'set neural.basePredictor to full or light to run a single base predictor.')
        }
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
        if (neural.deepGoPlusPlus.enabled) {
            predictors << new gspa.predictor.neural.DeepGoPlusPlusPredictor(
                sidecarScript: neural.sidecarScript,
                pythonExecutable: neural.pythonExecutable,
                integrator: neural.deepGoPlusPlus.integrator,
                componentsDir: neural.deepGoPlusPlus.componentsDir,
                dag: neural.deepGoPlusPlus.dag,
                batchSize: neural.deepGoPlusPlus.batchSize,
                minScore: neural.deepGoPlusPlus.minScore,
            )
        }
        if (neural.deepGoPlusPlusLight.enabled) {
            predictors << new gspa.predictor.neural.DeepGoPlusPlusLightPredictor(
                sidecarScript: neural.sidecarScript,
                pythonExecutable: neural.pythonExecutable,
                assets: neural.deepGoPlusPlusLight.assets,
                modelsDir: neural.deepGoPlusPlusLight.modelsDir,
                diamond: neural.deepGoPlusPlusLight.diamond,
                interpro: neural.deepGoPlusPlusLight.interpro,
                cnn: neural.deepGoPlusPlusLight.cnn,
                interproscan: neural.deepGoPlusPlusLight.interproscan,
                cnnModel: neural.deepGoPlusPlusLight.cnnModel,
                threads: neural.deepGoPlusPlusLight.threads,
                topK: neural.deepGoPlusPlusLight.topK,
                batchSize: neural.deepGoPlusPlusLight.batchSize,
                minScore: neural.deepGoPlusPlusLight.minScore,
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
        if (foss.mdf.enabled) {
            predictors << new gspa.predictor.structure.MdFPredictor(
                sidecarScript: foss.termSidecar,
                pythonExecutable: foss.pythonExecutable,
                weightsDir: foss.mdf.weightsDir,
                mdfExecutable: foss.mdf.mdfExecutable,
                skipPdb: foss.mdf.skipPdb,
                foldcompDb: foss.mdf.foldcompDb,
                threads: foss.mdf.threads,
                minScore: foss.mdf.minScore,
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

        // v1.3 viral / prophage predictors (genome-level)
        def viral = config.predictors.viral
        // These consume contig nucleotide sequence + coordinates, so they are
        // GenomePredictors and run in the genome stage, never on protein FASTA.
        List<AbstractGenomicRegionSidecarPredictor> regionPredictors = []
        if (viral.genomad.enabled) {
            regionPredictors << new gspa.predictor.viral.GenomadPredictor(
                sidecarScript: viral.genomicSidecar,
                pythonExecutable: viral.pythonExecutable,
                dbPath: viral.genomad.dbPath,
                genomadSif: viral.genomad.sif,
                minScore: viral.genomad.minScore,
            )
        }
        if (viral.checkv.enabled) {
            regionPredictors << new gspa.predictor.viral.CheckVPredictor(
                sidecarScript: viral.genomicSidecar,
                pythonExecutable: viral.pythonExecutable,
                dbPath: viral.checkv.dbPath,
                checkvSif: viral.checkv.sif,
                threads: viral.checkv.threads,
                minScore: viral.checkv.minScore,
            )
        }
        if (viral.phispy.enabled) {
            regionPredictors << new gspa.predictor.viral.PhiSpyPredictor(
                sidecarScript: viral.genomicSidecar,
                pythonExecutable: viral.pythonExecutable,
                phispySif: viral.phispy.sif,
                trainset: viral.phispy.trainset,
                minScore: viral.phispy.minScore,
            )
        }
        regionPredictors.each { predictors << it; genomePredictors << it }

        predictors
    }

    /**
     * Run the full annotation pipeline on a genome.
     */
    Genome annotate(File inputFile) {
        log.info("=== GSPA Annotation Pipeline ===")
        log.info("Input: ${inputFile}")
        runStatuses.clear()
        pathwaySearchReports.clear()

        // 1. Load genome
        Genome genome = loadGenome(inputFile)
        log.info("Loaded: ${genome}")

        // 2. Call genes if needed
        if (genome.proteinCount == 0 && config.input.requireExistingGenes) {
            throw new IllegalArgumentException('This workflow requires existing gene models. Use --proteins proteome.faa or --genome genome.fna --gff3 genes.gff3; GSPA has no splice-aware eukaryotic gene caller.')
        }
        if (genome.proteinCount == 0 && geneCaller != null) {
            log.info("Calling genes with ${geneCaller.name}...")
            geneCaller.callGenes(genome)
            log.info("Called ${genome.proteinCount} genes")
            // GeneCaller reports failure only by logging and returning no genes, so the exit
            // status of the caller is not observable here: zero genes is recorded as FAILED
            // (a gene caller that runs on a genome with sequence but finds no ORF is not a
            // result worth annotating on), with that limitation stated.
            runStatuses << (genome.proteinCount > 0 ?
                PredictorRunStatus.completed("gene-caller:${geneCaller.name}", genome.proteinCount) :
                PredictorRunStatus.failed("gene-caller:${geneCaller.name}",
                    'no genes called; caller failure vs genuinely gene-free input is not distinguishable (unverified)'))
        }

        if (genome.proteinCount == 0) {
            log.warn("No proteins found. Check input format or install a gene caller.")
            return genome
        }

        // 3. Run protein-level predictors (in parallel if multiple available)
        // Probe availability once; tools that are not runnable are recorded as
        // UNAVAILABLE so "never ran" is never confused with "ran and found nothing".
        def (List<Predictor> runnable, List<Predictor> missing) = registry.getAll().split { it.isAvailable() }
        missing.each { Predictor p ->
            String why = (p instanceof StatusReportingPredictor ? p.unavailableReason() : null) ?: 'tool not available'
            log.warn("  ${p.name} unavailable: ${why}")
            runStatuses << PredictorRunStatus.unavailable(p.name, why)
        }
        def availablePredictors = runnable.findAll { !(it instanceof GenomePredictor) }
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
                    def annotations = runTracked(predictor) { predictor.predictBatch(predictorInput) }
                    if (annotations != null) results[predictor.name] = annotations
                }
            }
            threads.each { it.join() }
            results.each { name, annotations ->
                applyAnnotations(genome, annotations as Map<String, List>, name)
            }
        } else {
            availablePredictors.each { predictor ->
                log.info("  Running ${predictor.name}...")
                def annotations = runTracked(predictor) { predictor.predictBatch(predictorInput) }
                if (annotations != null) applyAnnotations(genome, annotations, predictor.name)
            }
        }

        // Phase 10: propagate representative annotations to cluster members.
        if (clusterSet != null) {
            int copied = new ClusterAnnotationPropagator().propagate(clusterSet, [genome])
            log.info("  propagated ${copied} annotations from ${clusterSet.clusterCount()} reps to members")
        }

        // 4. Run genome-level predictors (sequential — they may depend on protein annotations)
        // Availability was probed above; genome-stage tools are the runnable
        // GenomePredictors (registered ones, plus any manually added and not disabled).
        List<GenomePredictor> availableGenomePredictors = runnable.findAll { it instanceof GenomePredictor } as List<GenomePredictor>
        genomePredictors.each { GenomePredictor gp ->
            if (registry.get(gp.name) == null && !config.predictors.disable.contains(gp.name) &&
                    !availableGenomePredictors.contains(gp)) {
                if (gp.isAvailable()) {
                    availableGenomePredictors << gp
                } else {
                    runStatuses << PredictorRunStatus.unavailable(gp.name, 'tool not available')
                }
            }
        }
        if (availableGenomePredictors) {
            log.info("Running ${availableGenomePredictors.size()} genome-level predictors...")
            availableGenomePredictors.each { predictor ->
                log.info("  Running ${predictor.name}...")
                def annotations = runTracked(predictor) { predictor.predictGenome(genome) }
                if (predictor instanceof gspa.predictor.pathway.GapseqPredictor && predictor.reactionTableSha256) {
                    pathwaySearchReports[predictor.name] = [schema:'gspa.pathway-search-observations/1',
                        genome_id:genome.id, source:predictor.name, target:predictor.target,
                        outcome:runStatuses.last().outcome.name(),
                        reaction_table_sha256:predictor.reactionTableSha256,
                        reaction_rows:predictor.reactionObservations.collect { new LinkedHashMap(it) }]
                }
                if (annotations != null) applyAnnotations(genome, annotations, predictor.name)
            }
        }

        // The predictors' claims, exactly as produced: before integration replaces annotations
        // and before any display threshold prunes weak candidates.
        captureRawClaims(genome)

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
     * Execute one predictor and record its typed completion status. Returns the
     * annotations, or null when the predictor threw (status FAILED). A status the
     * predictor reported itself (exit code, tool diagnostics) takes precedence;
     * otherwise it is derived from the result: empty = SUCCESS_ZERO.
     */
    private Map<String, List<Annotation>> runTracked(Predictor predictor, Closure<Map<String, List<Annotation>>> body) {
        PredictorRunStatus before = predictor instanceof StatusReportingPredictor ? predictor.lastRunStatus : null
        long t0 = System.currentTimeMillis()
        try {
            Map<String, List<Annotation>> result = body.call() ?: [:]
            PredictorRunStatus reported = predictor instanceof StatusReportingPredictor ? predictor.lastRunStatus : null
            PredictorRunStatus status = (reported != null && !reported.is(before)) ? reported :
                PredictorRunStatus.completed(predictor.name, result.values().sum { it.size() } as Integer ?: 0)
            if (status.elapsedMillis == 0L) status.elapsedMillis = System.currentTimeMillis() - t0
            runStatuses << status
            if (status.outcome == PredictorOutcome.FAILED) {
                log.error("  ${predictor.name} failed: ${status.detail}")
                return null
            }
            return result
        } catch (Exception e) {
            log.error("  ${predictor.name} failed: ${e.message}")
            runStatuses << PredictorRunStatus.failed(predictor.name, "${e.class.simpleName}: ${e.message}")
                .with { it.elapsedMillis = System.currentTimeMillis() - t0; it }
            return null
        }
    }

    /**
     * Run the shared inference workflow ({@link IntegrationWorkflow}) on a genome:
     * the same priors, references, context stages, enforcement and stage ledger as
     * {@code gspa integrate}. Claims default to the predictions on the genome's
     * proteins; raw {@link EvidenceClaim}s supplied by the caller are used directly
     * instead. State is keyed to the genome itself, and the final posteriors replace
     * its annotations. A requested stage whose data or tool is unavailable fails the
     * call; the outcome of every stage is in {@code lastIntegration.audit}.
     */
    IntegrationWorkflow.Result integrate(Genome genome, List<EvidenceClaim> claims = null, String observedAntiFamStatus = null) {
        log.info("=== Evidence integration (shared workflow) ===")
        List<EvidenceClaim> input = claims != null ? claims : new ClaimExtractor().fromGenome(genome)
        log.info("  ${input.size()} claims from ${genome.proteinCount} proteins")

        IntegrationWorkflowConfig workflowConfig = IntegrationWorkflowConfig.fromGspa(config)
        workflowConfig.genomeId = genome.id
        if (workflowConfig.outFile == null && integrationOutputDir != null) {
            workflowConfig.outFile = new File(integrationOutputDir, "${genome.id}_integrated.tsv")
        }
        def workflow = new IntegrationWorkflow(workflowConfig)
        workflow.antiFamScanStatus = observedAntiFamStatus
        workflow.sharedReferences = sharedIntegrationReferences
        workflow.info = { String m -> log.info(m) } as Closure<Void>
        lastIntegration = null
        lastIntegrationAudit = null
        IntegrationWorkflow.Result result
        try {
            result = workflow.run(input, genome)
        } finally {
            lastIntegrationAudit = workflow.lastAudit    // present even when a stage failed
        }
        workflow.writeArtifacts(result)
        lastIntegration = result
        log.info("  ${result.integrated}")
        result
    }

    private void runIntegration(Genome genome) {
        integrate(genome, rawClaims, runStatuses.find { it.predictor == 'antifam' }?.outcome?.name())
    }

    /**
     * Snapshot the claims carried by the genome's current annotations (source, evidence type, raw and
     * calibrated scores, metadata). Export them with {@link ClaimExtractor#writeClaimsJsonl} to run any
     * other route (web, Nextflow, {@code gspa integrate}) on identical evidence.
     */
    List<EvidenceClaim> captureRawClaims(Genome genome) {
        rawClaims = new ClaimExtractor().fromGenome(genome)
        rawClaims
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
     * Load genome from various input forms. Supported combinations:
     *
     * <ul>
     *   <li><b>nucleotide FASTA</b> (.fna/.fa/...) — contigs only; genes are
     *       called later (Phase 2). Multi-sequence FASTA = one contig each
     *       (metagenome).</li>
     *   <li><b>genome FASTA + GFF3</b> (config.input.gff3) — CDS features from
     *       the GFF3 are translated against the FASTA sequences, mapped to
     *       their contig. No gene calling.</li>
     *   <li><b>single GFF3 with embedded {@code ##FASTA}</b> — same, sequences
     *       come from the GFF3 itself.</li>
     *   <li><b>protein FASTA + GFF3</b> — supplied proteins are joined to
     *       contigs via the CDS attributes (protein_id/ID/locus_tag).</li>
     *   <li><b>protein FASTA only</b> — single synthetic contig (legacy).</li>
     * </ul>
     */
    private Genome loadGenome(File inputFile) {
        // Bare protein FASTA: no genome/GFF3, annotate on a single synthetic contig.
        if (config.input.proteinsOnly && config.input.proteinFasta) {
            File pf = new File(config.input.proteinFasta)
            Genome g = new Genome(id: pf.name.replaceAll(/\.(faa|fasta|fa|fsa)(\.gz)?$/, ''))
            Contig c = new Contig(id: 'proteins')
            g.addContig(c)
            FastaReader.readProteins(pf).each { c.addProtein(it) }
            g.mag = config.input.type in ['mag', 'metagenome']
            log.info("Loaded ${c.proteins.size()} proteins (protein-only input)")
            return g
        }

        String name = inputFile.name.toLowerCase()
        boolean primaryIsGff = name.endsWith('.gff3') || name.endsWith('.gff') ||
            name.endsWith('.gff3.gz') || name.endsWith('.gff.gz')

        File gffFile = config.input.gff3 ? new File(config.input.gff3) : (primaryIsGff ? inputFile : null)
        File seqFile = primaryIsGff ? null : inputFile
        File proteinFile = config.input.proteinFasta ? new File(config.input.proteinFasta) : null
        boolean isMeta = config.input.type in ['mag', 'metagenome']

        // --- GFF3-driven path: reuse supplied CDS coordinates ---
        if (gffFile) {
            Genome genome = GffReader.readGff3(gffFile)
            genome.mag = isMeta
            // Attach contig sequences from a separate genome FASTA if given
            // (a GFF3 with embedded ##FASTA already carries them).
            if (seqFile) {
                attachContigSequences(genome, seqFile)
            }
            if (proteinFile) {
                mapProteinsToContigs(genome, FastaReader.readProteins(proteinFile))
            } else {
                int n = CdsTranslator.populateProteins(genome, config.input.translateAltStart)
                log.info("Translated ${n} proteins from ${genome.features.count { it.type == FeatureType.CDS }} CDS features")
            }
            return genome
        }

        // --- nucleotide FASTA path: gene calling happens in Phase 2 ---
        Genome genome = FastaReader.readGenome(seqFile)
        genome.mag = isMeta

        if (proteinFile) {
            def proteins = FastaReader.readProteins(proteinFile)
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
     * Attach nucleotide sequences from a genome FASTA to the contigs parsed
     * from a GFF3 (matched by contig id == FASTA header first token). FASTA
     * contigs with no GFF3 features are added too, so assembly statistics and
     * per-contig metrics cover the whole assembly.
     */
    private void attachContigSequences(Genome genome, File seqFile) {
        int attached = 0
        FastaReader.parseFasta(seqFile).each { String id, String seq ->
            def contig = genome.findContig(id)
            if (contig) {
                contig.sequence = seq
                attached++
            } else {
                genome.addContig(new Contig(id: id, sequence: seq))
            }
        }
        log.info("Attached sequences to ${attached} contig(s) from ${seqFile.name}")
    }

    /**
     * Join pre-called proteins to contigs using the GFF3 CDS attributes
     * (protein_id / ID / locus_tag). Proteins with no matching feature land on
     * a synthetic 'unplaced' contig so they are still annotated.
     */
    private void mapProteinsToContigs(Genome genome, List<Protein> proteins) {
        Map<String, Feature> featureById = [:]
        genome.contigs.each { Contig c ->
            c.cdsFeatures().each { Feature f ->
                [f.attributes['protein_id']?.first(), f.id, f.locusTag].findAll { it }.each { String key ->
                    featureById.putIfAbsent(key, f)
                }
            }
        }
        Contig unplaced = null
        int placed = 0
        proteins.each { Protein p ->
            Feature f = featureById[p.id]
            if (f) {
                p.sourceFeature = f
                genome.findContig(f.seqId)?.addProtein(p)
                placed++
            } else {
                if (unplaced == null) {
                    unplaced = new Contig(id: 'unplaced')
                    genome.addContig(unplaced)
                }
                unplaced.addProtein(p)
            }
        }
        log.info("Joined ${placed}/${proteins.size()} proteins to contigs via GFF3 attributes" +
            (unplaced ? "; ${unplaced.proteins.size()} unplaced" : ''))
    }

    /**
     * Write annotated genome to output files.
     */
    void writeOutput(Genome genome, File outputDir) {
        outputDir.mkdirs()
        def manifestFile = new File(outputDir, "${genome.id}_predictor_status.json")
        PredictorRunManifest.write(genome.id, runStatuses, manifestFile)
        log.info("Predictor status manifest written: ${manifestFile}")
        // Reusing this genome's output directory must not certify an earlier
        // successful search when the current tool is unavailable or has failed.
        for (String source : ['gapsmith', 'gapseq']) {
            File previous = new File(outputDir, "${genome.id}_${source}_observations.json")
            if (!pathwaySearchReports.containsKey(source) && previous.exists() && !previous.delete()) {
                throw new IOException("Cannot remove stale pathway observations: ${previous}")
            }
        }
        pathwaySearchReports.each { source, report ->
            new File(outputDir, "${genome.id}_${source}_observations.json").text =
                groovy.json.JsonOutput.toJson(report) + '\n'
        }
        gspa.io.FastaWriter.writeProteins(genome.proteins, new File(outputDir, "${genome.id}.faa"))
        new File(outputDir, "${genome.id}_proteins.tsv").withWriter { w ->
            w.writeLine('protein_id\tcontig\tstart\tend\tstrand\tlocus_tag\tgene\tproduct\tproduct_source\thomology_accession\tidentity\tquery_coverage\tsubject_coverage\tlength_aa')
            genome.proteins.each { p ->
                def f = p.sourceFeature
                String product = f?.attributes?.get('product')?.first() ?: ''
                String origin = product ? 'input annotation' : ''
                def hit = p.annotations.annotations.findAll { it.type == gspa.model.AnnotationType.CUSTOM && it.metadata?.annotationType == 'product' }.max { it.metadata.bitscore ?: it.score }
                if ((!product || product.toLowerCase().startsWith('hypothetical') || product.toLowerCase().startsWith('uncharacterized') || product == 'Predicted CDS') && hit) {
                    product = hit.value
                    origin = 'Swiss-Prot homology (putative)'
                    if (f != null) {
                        f.attributes['product'] = [product]
                        f.attributes['product_source'] = [origin]
                        f.attributes['homology_accession'] = [hit.metadata.hitAccession ?: '']
                    }
                }
                def row = [p.id, p.contig?.id ?: '', p.start, p.end, p.strand.symbol,
                    f?.locusTag ?: p.id, f?.attributes?.get('gene')?.first() ?: '', product, origin,
                    hit?.metadata?.hitAccession ?: '', hit?.metadata?.pident ?: '',
                    hit?.metadata?.queryCoverage ?: '', hit?.metadata?.subjectCoverage ?: '', p.length]
                w.writeLine(row.collect { it.toString().replaceAll(/[\t\r\n]/, ' ') }.join('\t'))
            }
        }

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
        // The `source` column already records WHICH predictor produced each
        // annotation. When provenance is enabled, an extra `provenance` column
        // carries the trail (originating predictor + any enforcement actions),
        // so it is clear how each function was assigned.
        boolean prov = config.quality?.provenance
        output.withWriter { writer ->
            // `label` is the human-readable term name; `aspect` the GO sub-ontology
            // (MF/BP/CC) — both needed so a GO id is never shown bare, and for
            // standard-format (GAF) export / per-aspect summaries.
            def header = ['protein_id', 'type', 'value', 'label', 'score', 'source', 'evidence', 'aspect']
            if (prov) header << 'provenance'
            writer.writeLine(header.join('\t'))
            genome.proteins.each { protein ->
                protein.annotations.annotations.each { ann ->
                    def row = [
                        protein.id,
                        ann.type,
                        ann.value,
                        ann.goLabel ?: '',
                        String.format('%.4f', ann.score),
                        ann.source ?: '',
                        ann.evidence ?: '',
                        ann.goAspect ?: '',
                    ]
                    if (prov) {
                        def trail = ann.provenance ? ann.provenance.join(' | ') :
                            "predictor:${ann.source ?: 'unknown'}@${String.format('%.3f', ann.score)}"
                        row << trail
                    }
                    writer.writeLine(row.join('\t'))
                }
            }
        }
    }
}
