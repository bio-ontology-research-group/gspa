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
    boolean inferTaxon = false
    File operonsFile
    File gapsFile
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
    String gapseqTarget = 'genome'
    String intragenomeCluster = 'off'
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
    boolean enableAgent = false
    String agentModel = 'qwen3.8-27b'
    String agentPythonBin = 'python3'
    File agentSidecar = new File('benchmark/agent/run_gspa_agent.py')
    String agentApiBase
    boolean agentNoGuardrail = false

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
        if (gapseqTarget != 'genome' && gapseqTarget != 'proteome' && gapseqTarget != 'reps') {
            throw new IllegalArgumentException(
                "--gapseq-target must be one of: genome, proteome, reps (got '${gapseqTarget}')")
        }
        if (gapseqTarget == 'reps' && (intragenomeCluster == null || intragenomeCluster == 'off')) {
            throw new IllegalArgumentException(
                "--gapseq-target=reps requires --intragenome-cluster to be enabled (0.9 | 0.95 | 1.0)")
        }
        if (intragenomeCluster != 'off' && intragenomeCluster != '0.9' &&
                intragenomeCluster != '0.95' && intragenomeCluster != '1.0') {
            throw new IllegalArgumentException(
                "--intragenome-cluster must be one of: off, 0.9, 0.95, 1.0 (got '${intragenomeCluster}')")
        }
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
            operonsFile: operonsFile, gapsFile: gapsFile, orthogroupsFile: orthogroupsFile,
            clusterConsensusFile: clusterConsensusFile, genomeLayoutFile: genomeLayoutFile,
            thetaFile: thetaFile]
        named.each { String name, File f ->
            if (f != null && !f.isFile()) throw new IllegalArgumentException("${name} does not exist: ${f}")
        }
        (moduleFiles ?: []).each { File f ->
            if (f != null && !f.isFile()) throw new IllegalArgumentException("moduleFiles entry does not exist: ${f}")
        }
        if (pathwaysFile != null && ec2goFile == null) {
            throw new IllegalArgumentException("pathwaysFile requires ec2goFile (EC to GO mapping)")
        }
    }

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
