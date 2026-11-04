package gspa.integration.workflow

import com.fasterxml.jackson.annotation.JsonIgnore
import com.fasterxml.jackson.annotation.JsonAnySetter
import com.fasterxml.jackson.databind.DeserializationFeature
import com.fasterxml.jackson.databind.ObjectMapper

/**
 * Typed, serializable configuration of the shared integration workflow
 * ({@link IntegrationWorkflow}). One instance describes every route (CLI
 * {@code integrate}, batch, {@code annotate}, web, Nextflow): enabled stages
 * and priors, heavyweight reference files, per-genome inputs and per-genome
 * artifact paths. Unknown keys are rejected on read so a misspelled option can
 * never be silently ignored.
 *
 * <p>Reference files ({@code goOwlFile} ... {@code reactionGraphFile}) are loaded
 * once and shared across a batch; per-genome inputs ({@code operonsFile},
 * {@code gapsFile}, {@code genomeLayoutFile}, orthogroups, taxon, profile) and
 * artifacts ({@code outFile} and the {@code *Out} files) belong to one genome.</p>
 */
class IntegrationWorkflowConfig {
    static final Set<String> PRIOR_NAMES = ['essentiality', 'coherence', 'consistency', 'anatomical',
        'gap_filling', 'genomic_context', 'homology_transfer'].asImmutable()
    static final Set<String> ENFORCE_MODES = ['relax', 'remove', 'downrank', 'flag', 'minimal-flip'].asImmutable()
    static final Set<String> DARK_MATTER_MODES = ['separate', 'top-q', 'distributed'].asImmutable()
    static final String DEFAULT_PRIORS = 'essentiality,coherence,consistency,anatomical,gap_filling,genomic_context'
    private static final ObjectMapper MAPPER = new ObjectMapper()
        .configure(DeserializationFeature.FAIL_ON_UNKNOWN_PROPERTIES, true)

    // --- Evidence-model hyperparameters (theta.json equivalent) ---
    File thetaFile
    /** Inline theta overrides (merged over {@code thetaFile}); same keys as theta.json. */
    Map<String, Object> theta = [:]
    Integer maxIter
    Double epsilon
    Double damping

    // --- Shared heavyweight references ---
    File goOwlFile
    File essentialFunctionsFile
    String essentialProfile
    File ec2goFile
    File pathwaysFile
    File moduleDefsFile
    List<String> moduleExcludeOrganisms = []
    List<File> moduleFiles = []
    File taxonomyFile
    File taxonConstraintsFile
    File complexTermsFile
    boolean lite = false
    File reasonerCacheDir
    File reactionGraphFile
    File diffusionMetsFile
    File reactionEcAliasesFile
    double rlcCurrencyPercentile = 99.0d

    // --- Per-genome inputs ---
    String genomeId = 'genome'
    File claimsFile
    String organismTaxon
    /** Optional NCBI Taxonomy lineage IDs (comma/semicolon-separated) or TSV file mapping taxid -> lineage IDs. */
    String taxonLineage
    boolean inferTaxon = false
    File operonsFile
    File gapsFile
    File pathwayObservationsFile
    File orthogroupsFile
    File clusterConsensusFile
    File genomeLayoutFile
    File rxnLocusCatalogFile

    // --- Priors and stages ---
    String enabledPriors = DEFAULT_PRIORS
    boolean enforceStage2 = false
    String enforceConsistencyMode = 'relax'
    boolean darkMatter = false
    /** How dark-matter suggestions are applied to a genome: separate | top-q | distributed. */
    String darkMatterMode = 'separate'
    boolean iterateGapseq = false
    int maxGapseqIter = 5
    double gapseqTauCover = 0.5d
    double gapseqQBase = 0.5d
    double gapseqQStep = 0.05d
    double gapseqQCap = 0.75d
    boolean gapseqPinPromotions = true
    String promotionStrategy = 'default'
    int beamWidth = 5
    int beamCandidatesPerGap = 3
    double maxsatCoherenceBonus = 0.0d
    boolean refinedBf = true
    boolean rlcSuggester = false
    double rlcKernelBandwidth = 5000.0d
    int rlcRadiusK = 2
    double rlcAlpha = 0.5d
    double rlcAnchorThreshold = 0.3d
    double cgLambda = 1.0d
    int cgMinSupport = 3
    boolean cgRequireCredible = true
    boolean enableFbaRefiner = false
    File fbaDataDir
    File fbaSidecar = new File('scripts/stoichiometric_locus_refiner.py')
    String fbaPythonBin = 'python3'
    String fbaGram = 'neg'
    String fbaMedium = 'minimal'
    /** Stage 2B follow-up options (2026-10-31); all defaults keep the original Stage 2B behaviour. */
    boolean fbaSyntenyFromLayout = false
    File fbaKeggModules
    Double fbaCandidateMinScore
    Integer fbaCandidateTopK
    boolean enableAgent = false
    String agentModel = 'qwen3.8-27b'
    String agentPythonBin = 'python3'
    File agentSidecar = new File('benchmark/agent/run_gspa_agent.py')
    String agentApiBase
    boolean agentNoGuardrail = false
    boolean anomalyEcAware = false
    double anomalyFloor = 0.05d
    int anomalyMaxCandidates = 8
    int anomalyMaxLocus = 150

    // --- Learned stage, feature export and typed context ---
    /** Frozen production artifact of the learned context model. Naming it requests the learned stage. */
    File learnedModelFile
    /** Explicitly request the learned stage; without learnedModelFile this fails (no model is bundled, nothing is defaulted). */
    boolean learnedStage = false
    /** Remove a candidate whose learned probability is below this value (0..1); null never retracts. */
    Double learnedRetractBelow
    /** Add claimed functions that the preceding stages did not retain (use when the learned stage replaces evidence combination). */
    boolean learnedAddUnretainedClaims = false
    /** Export production feature rows (candidate_features.jsonl + contract + manifest) here; per genome in a batch. No model needed, no labels written. */
    File featureExportDir
    /** Strict JSON of typed ExternalContext (neighborhood, fba, roles, predictorStatus, assemblyCompleteness, candidates). */
    File externalContextFile

    // --- Joint cohort inference (two-pass in a batch) ---
    /** Request cohort context even without any of the cohort files below. */
    boolean cohortInference = false
    /** TSV of reference genomes that are peers only (genome_id, claims_path[, layout_path[, orthogroups_path]]). */
    File cohortManifestFile
    /** Cohort metadata JSON (clade, strain cluster, assembly, familySearchComplete, contigs, tree). Absent = unknown, treated conservatively. */
    File cohortMetadataFile
    /** Genome-aware orthogroup/subfamily TSV; without it each genome's --orthogroups file is used (cluster ids cohort-wide). */
    File cohortOrthogroupsFile
    File cohortRoleEvidenceFile
    File cohortInteractionEdgesFile
    File cohortReactionDirectionsFile
    File cohortPathwaysFile
    gspa.integration.cohort.CohortConfig cohortConfig = new gspa.integration.cohort.CohortConfig()

    // --- Metabolic experiments and community cross-feeding (conditional model results, never independent proof) ---
    static final Set<String> EXPERIMENT_MODES = ['manifest', 'evidence_guided'].asImmutable()
    static final Set<String> EXPERIMENT_STRATEGIES = ['fixed', 'random', 'evidence', 'qwen'].asImmutable()
    /**
     * Frozen experiment manifest (SBML-FBC model, baseline medium, explicit maintenance, optional context). Naming
     * it requests the stage. In mode {@code manifest} it runs exactly as supplied (a user-supplied analysis: the
     * ledger says so and claims no automated evidence guidance); in mode {@code evidence_guided} it is the TEMPLATE
     * (baseline, objectives, config) to which hypotheses and evidence records are added from this run's candidates.
     */
    File experimentManifestFile
    String experimentMode = 'manifest'
    /** Explicit transport/BGC/... evidence records referencing candidates by function key (evidence_guided only). */
    File experimentEvidenceFile
    String experimentStrategy = 'evidence'
    int experimentSeed = 0
    /** Experiments run after the baseline. */
    int experimentMaxExperiments = 8
    /** Hard cap on LP solves, baseline included. */
    int experimentMaxSolverCalls = 60
    /** Most candidate hypotheses added in evidence_guided mode; the rest are listed as truncated in the ledger. */
    int experimentMaxHypotheses = 20
    int experimentTimeoutSeconds = 600
    String experimentPythonBin = 'python3'
    File experimentSidecar = new File('scripts/metabolic_experiments.py')
    /** Defined community of separate organism models (batch-level: run once per batch). Never assigns an enzyme to a gene. */
    File communityManifestFile
    int communityMaxSolverCalls = 40
    int communityTimeoutSeconds = 600
    String communityPythonBin = 'python3'
    File communitySidecar = new File('scripts/community_metabolism.py')

    // --- Per-genome artifacts ---
    File outFile
    File provenanceFile
    File suggestionsOut
    File featuresOut
    File enforcementActionsOut
    File anomaliesOut
    File agentResolutionsOut

    /**
     * Any key that is not a property above. Rejected whichever mapper reads the file, including
     * the run-configuration loader that otherwise tolerates unknown keys.
     */
    @JsonAnySetter
    void unknownSetting(String name, Object value) {
        throw new IllegalArgumentException("Unknown integration workflow setting '${name}'")
    }

    /** Effective theta: file contents overlaid by inline entries and explicit scalar settings. */
    @JsonIgnore
    Map<String, Object> resolveTheta() {
        Map<String, Object> merged = [:]
        if (thetaFile != null && thetaFile.exists()) merged.putAll(MAPPER.readValue(thetaFile, Map))
        (theta ?: [:]).each { String key, Object value ->
            merged[key] = (merged[key] instanceof Map && value instanceof Map) ?
                ((merged[key] as Map) + (value as Map)) : value
        }
        if (maxIter != null) merged.max_iter = maxIter
        if (epsilon != null) merged.epsilon = epsilon
        if (damping != null) merged.damping = damping
        merged
    }

    @JsonIgnore
    Set<String> enabledPriorSet() {
        (enabledPriors ?: '').split(',').collect { it.trim() }.findAll { it } as Set<String>
    }

    /**
     * Reject inconsistent or unknown settings and explicitly requested files that
     * do not exist. Throws {@link IllegalArgumentException}; never downgrades a
     * requested stage to a silent no-op.
     */
    void validate() {
        if (iterateGapseq && !darkMatter) {
            throw new IllegalArgumentException(
                "--iterate-gapseq requires --dark-matter (the outer loop consumes suggester output)")
        }
        if (promotionStrategy != null && !(promotionStrategy in ['', 'default', 'greedy', 'maxsat', 'beam'])) {
            throw new IllegalArgumentException(
                "--promotion-strategy must be one of: default, greedy, maxsat, beam (got '${promotionStrategy}')")
        }
        if (beamWidth < 1) throw new IllegalArgumentException("--beam-width must be ≥ 1")
        if (beamCandidatesPerGap < 1) throw new IllegalArgumentException("--beam-candidates-per-gap must be ≥ 1")
        enabledPriorSet().findAll { it != 'none' && !PRIOR_NAMES.contains(it) }.each {
            throw new IllegalArgumentException("Unknown prior '${it}'. Known priors: ${PRIOR_NAMES.join(', ')}")
        }
        if (enforceConsistencyMode != null && !ENFORCE_MODES.contains(enforceConsistencyMode.toLowerCase(Locale.ROOT))) {
            throw new IllegalArgumentException("Unknown enforce-consistency-mode '${enforceConsistencyMode}'. Choose: ${ENFORCE_MODES.join(', ')}")
        }
        if (!DARK_MATTER_MODES.contains((darkMatterMode ?: '').toLowerCase(Locale.ROOT))) {
            throw new IllegalArgumentException("Unknown dark-matter mode '${darkMatterMode}'. Choose: ${DARK_MATTER_MODES.join(', ')}")
        }
        // An explicitly named reference must exist: a missing file is an error, not an absent prior.
        Map<String, File> named = [goOwlFile: goOwlFile, essentialFunctionsFile: essentialFunctionsFile,
            ec2goFile: ec2goFile, pathwaysFile: pathwaysFile, moduleDefsFile: moduleDefsFile,
            taxonomyFile: taxonomyFile, taxonConstraintsFile: taxonConstraintsFile,
            complexTermsFile: complexTermsFile, reactionGraphFile: reactionGraphFile,
            diffusionMetsFile: diffusionMetsFile, reactionEcAliasesFile: reactionEcAliasesFile,
            operonsFile: operonsFile, gapsFile: gapsFile, pathwayObservationsFile: pathwayObservationsFile, orthogroupsFile: orthogroupsFile,
            clusterConsensusFile: clusterConsensusFile, genomeLayoutFile: genomeLayoutFile,
            thetaFile: thetaFile, learnedModelFile: learnedModelFile, externalContextFile: externalContextFile,
            cohortManifestFile: cohortManifestFile, cohortMetadataFile: cohortMetadataFile,
            cohortOrthogroupsFile: cohortOrthogroupsFile, cohortRoleEvidenceFile: cohortRoleEvidenceFile,
            cohortInteractionEdgesFile: cohortInteractionEdgesFile, cohortReactionDirectionsFile: cohortReactionDirectionsFile,
            cohortPathwaysFile: cohortPathwaysFile,
            experimentManifestFile: experimentManifestFile, experimentEvidenceFile: experimentEvidenceFile,
            communityManifestFile: communityManifestFile]
        named.each { String name, File f ->
            if (f != null && !f.isFile()) throw new IllegalArgumentException("${name} does not exist: ${f}")
        }
        (moduleFiles ?: []).each { File f ->
            if (f != null && !f.isFile()) throw new IllegalArgumentException("moduleFiles entry does not exist: ${f}")
        }
        if (pathwaysFile != null && ec2goFile == null) {
            throw new IllegalArgumentException("pathwaysFile requires ec2goFile (EC to GO mapping)")
        }
        if (learnedStage && learnedModelFile == null) {
            throw new IllegalArgumentException("the learned stage was requested without a model (--learned-model): no trained model " +
                "is bundled and nothing is defaulted; omit the learned stage to run the legacy heuristic scoring")
        }
        if (learnedRetractBelow != null && !(learnedRetractBelow >= 0d && learnedRetractBelow <= 1d)) {
            throw new IllegalArgumentException("learnedRetractBelow must be in [0,1]")
        }
        if (learnedModelFile == null && (learnedRetractBelow != null || learnedAddUnretainedClaims)) {
            throw new IllegalArgumentException("learnedRetractBelow / learnedAddUnretainedClaims only apply to the learned stage (--learned-model)")
        }
        if (!EXPERIMENT_MODES.contains(experimentMode)) {
            throw new IllegalArgumentException("experimentMode must be one of ${EXPERIMENT_MODES.join(', ')} (got '${experimentMode}')")
        }
        if (!EXPERIMENT_STRATEGIES.contains(experimentStrategy)) {
            throw new IllegalArgumentException("experimentStrategy must be one of ${EXPERIMENT_STRATEGIES.join(', ')} (got '${experimentStrategy}')")
        }
        if (experimentMaxExperiments < 1 || experimentMaxSolverCalls < 1 || experimentMaxHypotheses < 1 ||
            experimentTimeoutSeconds < 1 || communityMaxSolverCalls < 1 || communityTimeoutSeconds < 1) {
            throw new IllegalArgumentException("experiment and community budgets and timeouts must be positive")
        }
        if (experimentManifestFile == null && (experimentEvidenceFile != null || experimentMode != 'manifest')) {
            throw new IllegalArgumentException("experimentMode / experimentEvidenceFile need experimentManifestFile (the manifest or template with baseline medium, objectives and maintenance)")
        }
        if (experimentEvidenceFile != null && experimentMode != 'evidence_guided') {
            throw new IllegalArgumentException("experimentEvidenceFile only applies to experimentMode evidence_guided")
        }
        cohortConfig.validate()
        if (cohortConfig.holdout == gspa.integration.cohort.HoldoutLevel.NONE) {
            throw new IllegalArgumentException("cohortConfig.holdout NONE would let a genome contribute to its own cohort features; use GENOME, STRAIN_CLUSTER or CLADE")
        }
        if (cohortRequested()) {
            if (cohortConfig.contextMode == gspa.integration.cohort.ContextMode.PATHWAY && cohortPathwaysFile == null) {
                throw new IllegalArgumentException('PATHWAY context requires --cohort-pathways')
            }
            if (cohortConfig.contextMode == gspa.integration.cohort.ContextMode.ORGANISM && cohortMetadataFile == null) {
                throw new IllegalArgumentException('ORGANISM context requires --cohort-metadata with plausibleReactions for every genome')
            }
            if (reactionGraphFile == null || reactionEcAliasesFile == null) {
                throw new IllegalArgumentException("cohort inference expands EC calls directly to reactions and needs both --reaction-graph " +
                    "and --reaction-ec-aliases (the EC bindings; the graph alone binds no EC to a reaction)")
            }
        } else if (cohortConfig.toMap() != new gspa.integration.cohort.CohortConfig().toMap()) {
            throw new IllegalArgumentException("cohortConfig was changed but cohort inference is not requested")
        }
    }

    /** Cohort context is requested explicitly, or implied by naming any cohort input. */
    @JsonIgnore
    boolean cohortRequested() {
        cohortInference || cohortManifestFile != null || cohortMetadataFile != null || cohortOrthogroupsFile != null ||
            cohortRoleEvidenceFile != null || cohortInteractionEdgesFile != null || cohortReactionDirectionsFile != null || cohortPathwaysFile != null
    }

    @JsonIgnore
    boolean experimentsRequested() { experimentManifestFile != null }

    @JsonIgnore
    boolean communityRequested() { communityManifestFile != null }

    /** The learned stage is requested by naming a model or asking for it. */
    @JsonIgnore
    boolean learnedRequested() { learnedStage || learnedModelFile != null }

    /**
     * The workflow configuration an annotation run implies: the explicit
     * {@code integration.workflow} block, with the scalar {@code integration.*}
     * settings and the run's own quality references filling whatever it leaves
     * unset. Equal inputs therefore give the same workflow as {@code gspa integrate}.
     */
    static IntegrationWorkflowConfig fromGspa(gspa.config.GspaConfig config) {
        gspa.config.GspaConfig.IntegrationConfig ic = config.integration
        IntegrationWorkflowConfig w = (ic.workflow ?: new IntegrationWorkflowConfig()).copy()
        if (w.maxIter == null) w.maxIter = ic.maxIter
        if (w.epsilon == null) w.epsilon = ic.epsilon
        if (w.damping == null) w.damping = ic.damping
        if (w.thetaFile == null && ic.thetaFile) w.thetaFile = new File(ic.thetaFile)
        if (ic.darkMatter?.enabled) {
            w.darkMatter = true
            w.darkMatterMode = ic.darkMatter.mode
            Map dm = [bf_min: ic.darkMatter.bfMin, gamma_in_p: ic.darkMatter.gammaInP,
                      coverage_threshold: ic.darkMatter.coverageThreshold]
            w.theta = (w.theta ?: [:]) + [dark_matter: ((w.theta?.dark_matter ?: [:]) as Map) + dm]
        }
        if (ic.outerLoop?.enabled) {
            w.iterateGapseq = true
            w.maxGapseqIter = ic.outerLoop.maxIter
            w.gapseqQBase = ic.outerLoop.qBase
            w.gapseqQStep = ic.outerLoop.qStep
            w.gapseqQCap = ic.outerLoop.qCap
            w.gapseqTauCover = ic.outerLoop.tauCover
            w.gapseqPinPromotions = ic.outerLoop.pinPromotions
        }
        if (w.goOwlFile == null && config.quality.goOwlFile) w.goOwlFile = new File(config.quality.goOwlFile)
        if (w.taxonomyFile == null && config.quality.consistency.taxonomyFile) {
            w.taxonomyFile = new File(config.quality.consistency.taxonomyFile)
        }
        if (w.taxonConstraintsFile == null && config.quality.consistency.constraintsFile) {
            w.taxonConstraintsFile = new File(config.quality.consistency.constraintsFile)
        }
        if (w.organismTaxon == null) w.organismTaxon = config.quality.consistency.organismTaxon
        if (w.essentialFunctionsFile == null && config.quality.completeness.customFile) {
            w.essentialFunctionsFile = new File(config.quality.completeness.customFile)
        }
        if (w.essentialProfile == null && !config.quality.completeness.customFile) {
            String profile = config.resolveCompletenessProfile()
            if (profile in ['bacteria', 'archaea', 'eukaryote', 'virus']) w.essentialProfile = profile
        }
        w
    }

    /**
     * Apply a settings overlay (the Nextflow route's {@code --integration_config}). Overlays carry
     * behavior only: a key that names a file or directory (those must be staged by the caller
     * into the task) or a per-genome binding ({@code genomeId}), and any unknown key, is rejected,
     * so an overlay can neither point at host paths the task cannot see nor redirect one sample's
     * inputs or outputs.
     */
    void applyOverlay(Map<String, Object> overlay) {
        (overlay ?: [:]).each { String key, Object value ->
            java.lang.reflect.Field property = IntegrationWorkflowConfig.declaredFields.find {
                it.name == key && !it.synthetic && !java.lang.reflect.Modifier.isStatic(it.modifiers)
            }
            if (property == null) throw new IllegalArgumentException("Unknown integration workflow setting '${key}'")
            if (File.isAssignableFrom(property.type) || key == 'moduleFiles' || key == 'genomeId') {
                throw new IllegalArgumentException(
                    "Overlay setting '${key}' names a file or per-genome binding; stage reference files into the task " +
                    "through the pipeline's own parameters instead of the overlay")
            }
            Object converted = MAPPER.convertValue(value, property.type)
            this[key] = key == 'theta' ? (((theta ?: [:]) as Map) + (converted as Map)) : converted
        }
    }

    IntegrationWorkflowConfig copy() {
        fromMap(toMap())
    }

    Map<String, Object> toMap() {
        MAPPER.convertValue(this, Map)
    }

    static IntegrationWorkflowConfig fromMap(Map<String, Object> values) {
        MAPPER.convertValue(values ?: [:], IntegrationWorkflowConfig)
    }

    static IntegrationWorkflowConfig fromJson(File file) {
        MAPPER.readValue(file, IntegrationWorkflowConfig)
    }

    void writeJson(File file) {
        file.parentFile?.mkdirs()
        MAPPER.writerWithDefaultPrettyPrinter().writeValue(file, toMap())
    }
}
