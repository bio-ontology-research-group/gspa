package gspa.cli

import com.fasterxml.jackson.databind.ObjectMapper
import gspa.integration.CalibrationTable
import gspa.integration.ClaimExtractor
import gspa.integration.EvidenceClaim
import gspa.integration.IntegrationState
import gspa.integration.workflow.IntegrationWorkflow
import gspa.integration.workflow.IntegrationWorkflowConfig
import picocli.CommandLine.Command
import picocli.CommandLine.Option

/**
 * {@code gspa integrate} — run the Phase 7 evidence integrator on a
 * pre-parsed claims file with a theta.json of hyperparameters.
 *
 * Designed to be the per-iteration evaluator invoked by the Phase 7.4
 * benchmark optimizer. Runs without the full predictor stack: reads claims.jsonl
 * and delegates every stage to the shared core
 * {@link gspa.integration.workflow.IntegrationWorkflow}; this class only parses
 * options into an {@link IntegrationWorkflowConfig}, reads claims and exports
 * (integrated TSV, provenance, stage ledger, batch statistics).
 */
@Command(
    name = 'integrate',
    description = 'Run the Phase 7 evidence integrator on pre-parsed claims.',
    mixinStandardHelpOptions = true
)
class IntegrateCommand implements Runnable {

    @picocli.CommandLine.Mixin PresetOptions preset = new PresetOptions()
    @picocli.CommandLine.Spec picocli.CommandLine.Model.CommandSpec spec

    void applyPreset() {
        if (!preset.name) return
        gspa.config.WorkflowPresets.integrationDefaults(preset.name).each { name, value ->
            def option = spec.findOption(name)
            // An explicit request to infer taxonomy must not be shadowed by
            // the broad domain taxon supplied by the workflow.
            if (name == '--taxon' && inferTaxon) return
            if (!spec.commandLine().parseResult.hasMatchedOption(option)) {
                option.setValue(value)
                presetApplied << option.longestName()
            }
        }
    }

    /** Long names of options whose value came from the preset rather than the command line. */
    private final Set<String> presetApplied = new LinkedHashSet<>()


    @Option(names = ['--claims'], required = false,
            description = 'JSONL file of evidence claims (one per line). Required unless --batch-manifest is given.')
    File claimsFile

    @Option(names = ['--theta'],
            description = 'JSON file with reliability + prior-weight hyperparameters.')
    File thetaFile

    @Option(names = ['--out'], required = false,
            description = 'Output TSV path. Required unless --batch-manifest is given.')
    File outFile

    @Option(names = ['--batch-manifest'],
            description = 'TSV manifest for multi-genome batch integration sharing loaded GO ontology, taxon/anatomical constraints, and pathway DB in memory across genomes.')
    File batchManifestFile

    @Option(names = ['--batch-threads'],
            description = 'Number of parallel worker threads when running --batch-manifest (default: 8).')
    int batchThreads = 8

    @Option(names = ['--batch-stats-out'],
            description = 'Optional TSV path to record per-genome wall-clock runtime (sec) and memory (MB) in --batch-manifest mode.')
    File batchStatsOut

    @Option(names = ['--provenance'],
            description = 'Optional provenance JSON output path.')
    File provenanceFile

    @Option(names = ['--go-owl'],
            description = 'Path to GO OWL file for ontology-driven priors.')
    File goOwlFile

    @Option(names = ['--essential-functions'],
            description = 'Essential function TSV (bacteria/archaea/eukaryote or custom).')
    File essentialFunctionsFile

    @Option(names = ['--essential-profile'],
            description = 'Built-in profile name: bacteria, archaea, eukaryote.')
    String essentialProfile

    @Option(names = ['--ec2go'], description = 'EC → GO mapping file.')
    File ec2goFile

    @Option(names = ['--pathways'], description = 'Pathway definitions TSV.')
    File pathwaysFile

    @Option(names = ['--module-defs'], description = 'KEGG module definitions JSON with step logic; enables the module-based pathway prior and writes <out>.modules.json.')
    File moduleDefsFile

    @Option(names = ['--module-exclude-organisms'], split = ',', description = 'KEGG organism codes left out of the module completion counts (organism context), for leave-one-out validation of a reference organism.')
    List<String> moduleExcludeOrganisms = []

    @Option(names = ['--modules'],
            description = 'Additional pathway TSV to stack on top of --pathways ' +
                          '(same schema). Designed for KEGG Modules — narrower units ' +
                          'than KEGG main pathways, so per-genome coverage is meaningful. ' +
                          'Repeatable: pass several files separated by commas to layer ' +
                          'KEGG modules + MetaCyc + BioCyc.', split = ',')
    List<File> moduleFiles = []

    @Option(names = ['--taxonomy'],
            description = 'NCBI taxonomy hierarchy file (for ConsistencyPrior).')
    File taxonomyFile

    @Option(names = ['--taxon-constraints'],
            description = 'GO taxon constraints OBO or TSV file.')
    File taxonConstraintsFile

    @Option(names = ['--taxon'],
            description = 'Assert organism NCBI taxon (e.g. NCBITaxon_2, 2, bacteria, eukaryote, NCBITaxon_4932).')
    String organismTaxon

    @Option(names = ['--infer-taxon'],
            description = 'Infer the organism taxon from high-confidence claims via GO + anatomical taxon constraints.')
    boolean inferTaxon = false

    @Option(names = ['--complex-terms'],
            description = 'Optional TSV of complex GO terms (term<TAB>h|n|a); defaults to bundled complex_terms.tsv.')
    File complexTermsFile

    @Option(names = ['--enforce-stage2'], negatable = true,
            description = 'Run Stage 2 Symbolic Constraint & Coherence Enforcement after Stage 1 refinement.')
    boolean enforceStage2 = false

    @Option(names = ['--enforce-consistency-mode'],
            description = 'Stage 2 consistency mode: relax (default), remove, downrank, flag, or minimal-flip.')
    String enforceConsistencyMode = 'relax'

    @Option(names = ['--enforcement-actions-out'],
            description = 'Optional TSV output path for Stage 2 enforcement actions.')
    File enforcementActionsOut

    @Option(names = ['--operons'],
            description = 'Operon assignments TSV (one operon per line, tab-separated protein IDs).')
    File operonsFile

    @Option(names = ['--gaps'],
            description = 'Metabolic gaps JSONL (one gap per line).')
    File gapsFile

    @Option(names = ['--orthogroups'],
            description = 'Orthogroup TSV (protein_id <TAB> cluster_id), from MMseqs2 easy-cluster across multiple genomes. Required with homology_transfer prior.')
    File orthogroupsFile

    @Option(names = ['--cluster-consensus'],
            description = 'Cluster consensus TSV (cluster_id <TAB> function_type <TAB> function_id <TAB> consensus_prob) built from a prior all-genome baseline pass.')
    File clusterConsensusFile

    @Option(names = ['--enable-priors'],
            description = 'Which priors to enable (comma-separated). Default: all.')
    String enabledPriors = 'essentiality,coherence,consistency,anatomical,gap_filling,genomic_context'

    @Option(names = ['--lite'], negatable = true, description = 'Skip ELK initialization (no process coherence).')
    boolean lite = false

    @Option(names = ['--reasoner-cache'],
            description = 'Directory for caching expensive reasoner results (has_part pairs).')
    File reasonerCacheDir

    @Option(names = ['--dark-matter'],
            description = 'Enable the Phase 8 dark-matter / contextual-gap suggester.')
    boolean darkMatter = false

    @Option(names = ['--suggestions-out'],
            description = 'Suggestions TSV output path (only with --dark-matter).')
    File suggestionsOut

    // --- Phase 10 outer-loop flags ---

    @Option(names = ['--iterate-gapseq'],
            description = 'Enable the Phase 10 outer fixed-point loop (requires --dark-matter).')
    boolean iterateGapseq = false

    @Option(names = ['--max-gapseq-iter'],
            description = 'Max outer-loop iterations (default 5).')
    int maxGapseqIter = 5

    @Option(names = ['--gapseq-tau-cover'],
            description = 'Posterior-probability threshold for coverage-based gap recomputation (default 0.5).')
    double gapseqTauCover = 0.5d

    @Option(names = ['--gapseq-q-base'],
            description = 'Base q threshold for DarkMatter promotion (default 0.5).')
    double gapseqQBase = 0.5d

    @Option(names = ['--gapseq-q-step'],
            description = 'Per-iteration q-threshold increment (default 0.05).')
    double gapseqQStep = 0.05d

    @Option(names = ['--gapseq-pin-promotions'], arity = '1',
            description = 'Pin promoted singletons as posterior floors (default true).')
    boolean gapseqPinPromotions = true

    @Option(names = ['--gapseq-target'],
            description = 'Unsupported in claims-only integration (only the default "genome" is accepted): gapseq runs ' +
                          'outside this command, so a proteome/reps search target cannot be applied here.')
    String gapseqTarget = 'genome'

    @Option(names = ['--intragenome-cluster'],
            description = 'Unsupported in claims-only integration (only "off" is accepted): clustering needs ' +
                          'sequences. Use gspa annotate --intragenome-cluster before inference.')
    String intragenomeCluster = 'off'

    @Option(names = ['--promotion-strategy'],
            description = 'How the outer loop picks promotions from DarkMatter suggestions: ' +
                          'default (all above q threshold), greedy (log-posterior rank, ' +
                          'conflict-free batch), maxsat (SAT4J weighted MaxSAT), ' +
                          'beam (top-k per gap with beam search over assignments).')
    String promotionStrategy = 'default'

    @Option(names = ['--beam-width'], description = 'BeamSearchStrategy: beam width W. Default 5.')
    int beamWidth = 5

    @Option(names = ['--beam-candidates-per-gap'],
            description = 'BeamSearchStrategy: top-k candidates per gap. Default 3.')
    int beamCandidatesPerGap = 3

    @Option(names = ['--gapseq-q-cap'],
            description = 'Upper cap on the rising q threshold across outer iterations. Default 0.75.')
    double gapseqQCap = 0.75d

    @Option(names = ['--maxsat-coherence-bonus'],
            description = 'MaxSatStrategy: weight on pairwise pathway-coherence bonus. ' +
                          'Zero (default) disables; positive values reward jointly committing ' +
                          'candidates in the same pathway (makes MaxSAT diverge from greedy).')
    double maxsatCoherenceBonus = 0.0d

    @Option(names = ['--refined-bf'], arity = '0..1', fallbackValue = 'true',
            description = 'Use the Phase 10.1 refined BF (Noisy-OR + IC + purity) in the ' +
                          'suggester. Accepts true/false; bare --refined-bf = true. Default true.')
    boolean refinedBf = true

    // --- Phase 12 RLGC flags ---

    @Option(names = ['--reaction-graph'],
            description = 'gapsmith seed_reactions.tsv — builds the reaction graph for the RLGC suggester.')
    File reactionGraphFile

    @Option(names = ['--diffusion-mets'],
            description = 'gapsmith diffusion_mets.tsv (currency metabolite list).')
    File diffusionMetsFile

    @Option(names = ['--reaction-ec-aliases'],
            description = 'gapsmith seed_Enzyme_Class_Reactions_Aliases_unique.tsv (EC → rxn binding).')
    File reactionEcAliasesFile

    @Option(names = ['--genome-layout'],
            description = 'Panel layout TSV (protein_id contig start end strand) built by make_layout.py.')
    File genomeLayoutFile

    @Option(names = ['--rlc-suggester'],
            description = 'Enable the Phase 12 Reaction-Local Context suggester (M1).')
    boolean rlcSuggester = false

    @Option(names = ['--rlc-kernel-bandwidth'],
            description = 'Gaussian kernel bandwidth (bp) for genomic density. Default 5000.')
    double rlcKernelBandwidth = 5000.0d

    @Option(names = ['--rlc-radius-k'],
            description = 'Reaction-graph BFS radius. Default 2.')
    int rlcRadiusK = 2

    @Option(names = ['--rlc-alpha'],
            description = 'α-decay across reaction-graph hops. Default 0.5.')
    double rlcAlpha = 0.5d

    @Option(names = ['--rlc-currency-pct'],
            description = 'Percentile threshold for degree-based currency detection. Default 99.0.')
    double rlcCurrencyPercentile = 99.0d

    @Option(names = ['--rlc-anchor-threshold'],
            description = 'Min posterior prob to qualify as an RLGC anchor. Default 0.3.')
    double rlcAnchorThreshold = 0.3d

    @Option(names = ['--features-out'],
            description = 'Optional TSV output: per-candidate feature vector from the RLGC suggester (M3 training data).')
    File featuresOut

    // --- Phase 12 M2 cross-genome LR flags ---

    @Option(names = ['--rxn-locus-catalog'],
            description = 'ReactionLocusCatalog TSV (from build_catalog.py). Enables M2 cross-genome re-scoring on top of RLGC.')
    File rxnLocusCatalogFile

    @Option(names = ['--cg-lambda'],
            description = 'Exponent on cross-genome LR in posterior update. Default 1.0.')
    double cgLambda = 1.0d

    @Option(names = ['--cg-min-support'],
            description = 'Min n_sig_total for a (C, R) LR to be trusted. Default 3.')
    int cgMinSupport = 3

    @Option(names = ['--cg-require-credible'], arity = '1',
            description = 'Drop LRs whose 90% CI overlaps 1.0 (log-CI excludes 0). Default true.')
    boolean cgRequireCredible = true

    // --- Stage 3 Optional Neuro-Symbolic Agentic Completion flags ---

    @Option(names = ['--enable-agent'],
            description = 'Enable the optional Stage 3 Neuro-Symbolic Agentic Annotation Completion Layer (gspa-agent).')
    boolean enableAgent = false

    @Option(names = ['--agent-model'],
            description = 'Model for Stage 3 Agent: qwen3.8-27b (default, on unimatrix01) or qwen-small.')
    String agentModel = 'qwen3.8-27b'

    @Option(names = ['--agent-python-bin'], description = 'Python interpreter with the agent dependencies installed.')
    String agentPythonBin = 'python3'

    @Option(names = ['--agent-sidecar'],
            description = 'Path to benchmark/agent/run_gspa_agent.py.')
    File agentSidecar = new File('benchmark/agent/run_gspa_agent.py')

    @Option(names = ['--agent-api-base'],
            description = 'Optional OpenAI-compatible API base URL for Stage 3 Agent.')
    String agentApiBase = null

    @Option(names = ['--agent-no-guardrail'],
            description = 'Disable the deterministic Symbolic Guardrail (for guarded vs. unguarded ablation only).')
    boolean agentNoGuardrail = false

    @Option(names = ['--anomalies-out'],
            description = 'Output JSONL path for extracted genome-scale anomalies.')
    File anomaliesOut

    @Option(names = ['--agent-resolutions-out'],
            description = 'Output JSONL path for Stage 3 agent resolutions.')
    File agentResolutionsOut

    // --- Stage 2B: StoichiometricLocusRefiner (FBA-coupled metabolic repair & locus attribution) ---

    @Option(names = ['--enable-fba-refiner'],
            description = 'Enable Stage 2B StoichiometricLocusRefiner (Positive-Cost Flux-Certified Repair + pFBA + Backward Locus Responsibility).')
    boolean enableFbaRefiner = false

    @Option(names = ['--fba-data-dir'],
            description = 'Directory containing universal_metadata.json, fba_matrices_v6RE.pkl, tight_bounds_v6_{neg,pos}.pkl, seedr2ec.pkl, seedec2r.pkl, medium.pkl.')
    File fbaDataDir

    @Option(names = ['--fba-sidecar'],
            description = 'Path to scripts/stoichiometric_locus_refiner.py.')
    File fbaSidecar = new File('scripts/stoichiometric_locus_refiner.py')

    @Option(names = ['--fba-python-bin'],
            description = 'Python interpreter for StoichiometricLocusRefiner sidecar.')
    String fbaPythonBin = 'python3'

    @Option(names = ['--fba-gram'],
            description = 'Gram stain biomass template for StoichiometricLocusRefiner: neg (default) or pos.')
    String fbaGram = 'neg'

    @Option(names = ['--fba-medium'],
            description = 'Medium name for StoichiometricLocusRefiner (default: minimal).')
    String fbaMedium = 'minimal'
    @Option(names = ['--genome-id'],
            description = 'Identifier of the genome in a single run (cohort metadata and reports are keyed by it). Default: genome.')
    String genomeId = 'genome'

    // --- Learned stage, feature export and typed context ---

    @Option(names = ['--learned-model'],
            description = 'Frozen production artifact of the learned context model. Runs the learned stage LAST over the original raw ' +
                          'claims; its calibrated probabilities REPLACE the heuristic scores. Fails if the file is missing, not a ' +
                          'production artifact, or built for another feature contract. Without it the scoring is the legacy heuristic ' +
                          '(marked as such in the ledger); no trained model is bundled.')
    File learnedModelFile

    @Option(names = ['--learned-stage'],
            description = 'Explicitly request the learned stage; fails without --learned-model instead of silently using the heuristic scoring.')
    boolean learnedStage = false

    @Option(names = ['--learned-retract-below'],
            description = 'Remove a candidate whose learned probability is below this value (0..1). Default: never retract.')
    Double learnedRetractBelow

    @Option(names = ['--learned-add-unretained-claims'],
            description = 'Also add claimed functions that the preceding stages did not retain (use when the learned stage replaces evidence combination).')
    boolean learnedAddUnretainedClaims = false

    @Option(names = ['--feature-export-dir'],
            description = 'Write production feature rows (candidate_features.jsonl, feature_contract.json, export_manifest.json) here, one ' +
                          'sub-directory per genome in a batch. Needs no model; writes no labels, splits or groups.')
    File featureExportDir

    @Option(names = ['--external-context'],
            description = 'Strict JSON of typed ExternalContext (neighborhood, fba, roles, predictorStatus, assemblyCompleteness, candidates). ' +
                          'Keys also supplied by cohort inference are an error.')
    File externalContextFile

    // --- Joint cohort inference ---

    @Option(names = ['--cohort'],
            description = 'Compute cohort recurrence context. In --batch-manifest mode this is a two-pass run: stage-0 raw baseline for ALL genomes, ' +
                          'one joint inference, then each genome is scored with its held-out features. Needs --reaction-graph (EC -> reaction).')
    boolean cohortInference = false

    @Option(names = ['--cohort-manifest'],
            description = 'TSV of reference genomes that are peers only: genome_id, claims_path[, layout_path[, orthogroups_path]]. Implies --cohort.')
    File cohortManifestFile

    @Option(names = ['--cohort-metadata'],
            description = 'Cohort metadata JSON (clade, strainCluster, assemblyStatus, familySearchComplete, contigs, tree). Absent facts are ' +
                          'treated as unknown, never as independent. Implies --cohort.')
    File cohortMetadataFile

    @Option(names = ['--cohort-orthogroups'],
            description = 'Genome-aware orthogroup TSV: genome_id, protein_id, orthogroup_id[, subfamily_id, weight]. Without it each genome\'s ' +
                          '--orthogroups file (protein_id, cluster_id; cohort-wide cluster ids) is used.')
    File cohortOrthogroupsFile

    @Option(names = ['--cohort-role-evidence'],
            description = 'TSV of protein-intrinsic raw role scores: genome_id, protein_id, source, role, score.')
    File cohortRoleEvidenceFile

    @Option(names = ['--cohort-interaction-edges'],
            description = 'TSV of independent curated reaction-reaction edges: reaction_a, reaction_b, source[, directed].')
    File cohortInteractionEdgesFile

    @Option(names = ['--cohort-reaction-directions'],
            description = 'TSV of reaction directions: reaction_id, direction (FORWARD, REVERSE, REVERSIBLE, UNKNOWN).')
    File cohortReactionDirectionsFile

    @Option(names = ['--cohort-holdout'],
            description = 'Which genomes are removed from the peers of a target: GENOME, STRAIN_CLUSTER (default) or CLADE.')
    String cohortHoldout

    @Option(names = ['--cohort-half-width-bp'], description = 'Half-width of the neighbourhood window (default 5000).')
    Integer cohortHalfWidthBp

    @Option(names = ['--cohort-permutations'], description = 'Permutations of the matched-window null (default 200; 0 disables).')
    Integer cohortPermutations

    @Option(names = ['--cohort-seed'], description = 'Seed of the permutation null (default 20261026).')
    Long cohortSeed

    @Option(names = ['--cohort-tree-cut-distance'], description = 'Clade cut height on a supplied tree (default 0.05).')
    Double cohortTreeCutDistance

    @Option(names = ['--workflow-config'],
            description = 'JSON file with the full typed integration workflow configuration ' +
                          '(IntegrationWorkflowConfig). Precedence: --preset defaults < file < explicit command-line options. ' +
                          'Unknown keys are rejected.')
    File workflowConfigFile

    @Option(names = ['--workflow-overlay'],
            description = 'JSON of behavior-only workflow settings applied over everything else (used by the Nextflow ' +
                          'route). Keys naming files/directories or the per-genome binding genomeId, and unknown keys, are rejected.')
    File workflowOverlayFile

    @Override
    void run() {
        applyPreset()
        if (preset.printConfig) {
            Map settings = spec.options().findAll {
                !it.usageHelp() && !it.versionHelp() && !(it.longestName() in ['--preset', '--print-config'])
            }.collectEntries { option -> [(option.longestName()): option.getValue()] }
            spec.commandLine().out.println(new ObjectMapper().writerWithDefaultPrettyPrinter().writeValueAsString(settings))
            return
        }
        if (preset.name) println "Preset: ${gspa.config.WorkflowPresets.canonicalName(preset.name)}"
        if (intragenomeCluster != 'off') {
            throw new IllegalArgumentException("--intragenome-cluster ${intragenomeCluster} is not supported by claims-only integration " +
                "(clustering needs the proteome sequences); cluster during gspa annotate with --intragenome-cluster, before inference.")
        }
        if (gapseqTarget != 'genome') {
            throw new IllegalArgumentException("--gapseq-target ${gapseqTarget} is not supported by claims-only integration " +
                "(gapseq runs outside this command); only the default 'genome' is accepted.")
        }
        // Validation belongs to the workflow's reference_preparation stage so that a failure
        // is recorded in the genome's ledger.
        IntegrationWorkflowConfig config = toWorkflowConfig()
        if (batchManifestFile != null) {
            runBatchManifest(config)
            return
        }
        if (config.claimsFile == null || config.outFile == null) {
            throw new IllegalArgumentException("Either --batch-manifest or both --claims and --out must be provided.")
        }
        runGenome(config, null)
    }

    /**
     * The typed workflow configuration this invocation describes. Every option
     * whose name matches a {@link IntegrationWorkflowConfig} property maps onto it;
     * with {@code --workflow-config} the file is the base and only options given
     * explicitly on the command line override it.
     */
    IntegrationWorkflowConfig toWorkflowConfig() {
        IntegrationWorkflowConfig config = workflowConfigFile != null ?
            IntegrationWorkflowConfig.fromJson(workflowConfigFile) : new IntegrationWorkflowConfig()
        // Precedence with a file: preset defaults < file < explicit command-line options.
        Set<String> fileKeys = workflowConfigFile == null ? null :
            new ObjectMapper().readValue(workflowConfigFile, Map).keySet() as Set<String>
        IntegrateCommand.declaredFields.findAll { it.getAnnotation(Option) != null }.each { field ->
            if (!IntegrationWorkflowConfig.metaClass.hasProperty(config, field.name)) return
            if (fileKeys != null) {
                String[] names = field.getAnnotation(Option).names()
                boolean explicit = names.any { spec.commandLine().parseResult.hasMatchedOption(it) }
                boolean fromPreset = names.any { presetApplied.contains(it) } && !fileKeys.contains(field.name)
                if (!explicit && !fromPreset) return
            }
            field.accessible = true
            def value = field.get(this)
            config[field.name] = value instanceof List ? new ArrayList(value) : value
        }
        // Cohort engine settings are nested in the workflow config; a value given on the command line wins.
        if (cohortHoldout != null) {
            try { config.cohortConfig.holdout = gspa.integration.cohort.HoldoutLevel.valueOf(cohortHoldout.toUpperCase(Locale.ROOT)) }
            catch (IllegalArgumentException e) {
                throw new IllegalArgumentException("--cohort-holdout must be one of GENOME, STRAIN_CLUSTER, CLADE (got '${cohortHoldout}')")
            }
            if (config.cohortConfig.holdout == gspa.integration.cohort.HoldoutLevel.NONE) {
                throw new IllegalArgumentException("--cohort-holdout NONE would let a genome score itself; use GENOME, STRAIN_CLUSTER or CLADE")
            }
        }
        if (cohortHalfWidthBp != null) config.cohortConfig.halfWidthBp = cohortHalfWidthBp
        if (cohortPermutations != null) config.cohortConfig.permutations = cohortPermutations
        if (cohortSeed != null) config.cohortConfig.seed = cohortSeed
        if (cohortTreeCutDistance != null) config.cohortConfig.treeCutDistance = cohortTreeCutDistance
        if (workflowOverlayFile != null) {
            config.applyOverlay(new ObjectMapper().readValue(workflowOverlayFile, Map))
        }
        config
    }

    /** Runs the shared workflow for one genome and exports its artifacts. */
    private Map runGenome(IntegrationWorkflowConfig config, IntegrationState sharedReferences,
                          gspa.integration.contextbridge.CohortContextBuilder.Result cohort = null) {
        println "GSPA integrate"
        println "  Claims: ${config.claimsFile}"
        println "  Theta:  ${config.thetaFile ?: '(defaults)'}"
        println "  Out:    ${config.outFile}"

        def extractor = new ClaimExtractor(new CalibrationTable())
        List<EvidenceClaim> claims = extractor.readClaimsJsonl(config.claimsFile)
        println "  Loaded ${claims.size()} claims"

        def workflow = new IntegrationWorkflow(config)
        workflow.sharedReferences = sharedReferences
        workflow.batchCohort = cohort
        workflow.info = { String m -> println m } as Closure<Void>
        IntegrationWorkflow.Result result = workflow.run(claims)
        workflow.writeArtifacts(result)

        Runtime rt = Runtime.getRuntime()
        [genome_id: result.genomeId, out_path: config.outFile.path,
         n_claims: claims.size(), n_integrated: result.integrated.annotations.size(),
         n_provenance: result.integrated.provenance.size(), max_iter: result.maxIter,
         enforce_stage2: config.enforceStage2, enforce_mode: config.enforceConsistencyMode,
         runtime_sec: result.elapsedNanos / 1e9d,
         used_heap_mb: (rt.totalMemory() - rt.freeMemory()) / (1024d * 1024d),
         total_heap_mb: rt.totalMemory() / (1024d * 1024d)]
    }

    /**
     * Batch execution runs the same workflow as a single genome. Immutable
     * ontology/pathway/reaction references are loaded once and shared; taxon
     * checkers, genome evidence, options and artifact paths belong to each worker.
     *
     * Manifest TSV columns (tab-separated, '-' or empty uses CLI default):
     *   0: genome_id (or tag)
     *   1: claims_path
     *   2: out_path
     *   3: taxon (optional, e.g. NCBITaxon_562, bacteria, archaea, virus, eukaryote)
     *   4: essential_profile (optional, e.g. bacteria, archaea, virus, eukaryote)
     *   5: operons_path (optional)
     *   6: gaps_path (optional)
     *   7: max_iter (optional integer override, e.g. 1, 2, 3)
     *   8: enforce_stage2 (optional true/false override)
     *   9: enforce_mode (optional relax/remove override)
     *  10: enable_priors (optional comma-separated prior list override)
     *  11: genome_layout_path (optional)
     *  12: orthogroups_path (optional)
     *  13: cluster_consensus_path (optional)
     *  14: fba_medium (optional)
     *  15: fba_gram (optional)
     */
    /** One genome's configuration: the batch settings, its manifest row, and artifact paths of its own. */
    private static IntegrationWorkflowConfig workerConfig(IntegrationWorkflowConfig base, String[] row) {
        IntegrationWorkflowConfig worker = base.copy()
        worker.genomeId = row[0].trim()
        worker.claimsFile = new File(row[1].trim())
        worker.outFile = new File(row[2].trim())
        Map<Integer, String> overrides = [3: 'organismTaxon', 4: 'essentialProfile',
            5: 'operonsFile', 6: 'gapsFile', 8: 'enforceStage2',
            9: 'enforceConsistencyMode', 10: 'enabledPriors', 11: 'genomeLayoutFile',
            12: 'orthogroupsFile', 13: 'clusterConsensusFile', 14: 'fbaMedium', 15: 'fbaGram']
        overrides.each { int column, String property ->
            if (row.length > column && row[column].trim() && row[column].trim() != '-') {
                String value = row[column].trim()
                if (property.endsWith('File')) worker[property] = new File(value)
                else if (property == 'enforceStage2') {
                    if (!(value in ['true', 'false'])) throw new IllegalArgumentException("Invalid enforce_stage2: ${value}")
                    worker[property] = Boolean.parseBoolean(value)
                } else worker[property] = value
            }
        }
        if (row.length > 7 && row[7].trim() && row[7].trim() != '-') worker.maxIter = Integer.parseInt(row[7].trim())
        if (worker.featureExportDir != null) {
            if (worker.genomeId.contains('/') || worker.genomeId.contains('\\') || worker.genomeId in ['.', '..']) {
                throw new IllegalArgumentException("genome id '${worker.genomeId}' cannot name a feature-export sub-directory")
            }
            worker.featureExportDir = new File(base.featureExportDir, worker.genomeId)
        }
        // Every worker gets separate files, even when outputs share a directory.
        [provenanceFile: '.provenance.json', suggestionsOut: '.suggestions.tsv',
         featuresOut: '.features.tsv', enforcementActionsOut: '.enforcement.tsv',
         anomaliesOut: '.anomalies.jsonl', agentResolutionsOut: '.agent_resolutions.jsonl'].each { property, suffix ->
            if (worker[property] != null) worker[property] = new File(worker.outFile.absolutePath + suffix)
        }
        worker
    }

    private void runBatchManifest(IntegrationWorkflowConfig base) {
        println "GSPA integrate --batch-manifest ${batchManifestFile} (threads=${batchThreads})"

        List<String[]> rows = []
        Set<String> genomeIds = new HashSet<>()
        Set<String> outputPaths = new HashSet<>()
        batchManifestFile.eachLine { line, idx ->
            if (!line.trim() || line.trim().startsWith('#')) return
            String[] parts = line.split('\t', -1)
            if (parts[0].trim() in ['genome_id', 'tag', 'id'] && parts.length > 1 &&
                parts[1].trim() in ['claims_path', 'claims']) return
            if (parts.length < 3 || !parts[0].trim() || !parts[1].trim() || !parts[2].trim()) {
                throw new IllegalArgumentException("Invalid batch manifest row ${idx}: genome_id, claims_path and out_path are required")
            }
            if (!genomeIds.add(parts[0].trim()) || !outputPaths.add(new File(parts[2].trim()).canonicalPath)) {
                throw new IllegalArgumentException("Duplicate genome ID or output path in batch manifest row ${idx}")
            }
            rows << parts
        }
        println "  Loaded ${rows.size()} batch tasks from ${batchManifestFile}"
        IntegrationState sharedState
        try {
            sharedState = IntegrationWorkflow.loadSharedReferences(base) { String m -> println m }
        } catch (Exception e) {
            // No genome can run without the shared references: give each its own failed ledger
            // (reference_preparation failed, requested stages not run) and fail the batch once.
            rows.each { String[] row ->
                IntegrationWorkflow.recordPreparationFailure(workerConfig(base, row), e)
            }
            throw new IllegalStateException("Batch not run: shared reference loading failed (${e.message}); " +
                "${rows.size()} genome ledgers record the failure", e)
        }
        // Pass 1 of a cohort run: the stage-0 raw baseline of EVERY genome, then ONE joint inference.
        // Each genome's learned stage afterwards consumes only its own held-out records.
        gspa.integration.contextbridge.CohortContextBuilder.Result cohort = null
        if (base.cohortRequested()) {
            try {
                def stage = new gspa.integration.contextbridge.CohortStage()
                stage.info = { String m -> println m } as Closure<Void>
                cohort = stage.run(base, gspa.integration.contextbridge.CohortStage.targetsOf(rows.collect { workerConfig(base, it) }),
                    sharedState.reactionGraph)
            } catch (Exception e) {
                rows.each { String[] row ->
                    IntegrationWorkflow.recordPreparationFailure(workerConfig(base, row), e, 'cohort_context')
                }
                throw new IllegalStateException("Batch not run: the cohort pass failed (${e.message}); " +
                    "${rows.size()} genome ledgers record the failure", e)
            }
        }
        def pool = java.util.concurrent.Executors.newFixedThreadPool(Math.max(1, batchThreads))
        List<Map> stats = []
        try {
            List<java.util.concurrent.Future<Map>> futures = rows.collect { String[] row ->
                pool.submit({ -> runGenome(workerConfig(base, row), sharedState, cohort) } as java.util.concurrent.Callable<Map>)
            }
            futures.each { stats << it.get() }
        } finally {
            pool.shutdownNow()
        }

        if (batchStatsOut != null) {
            batchStatsOut.parentFile?.mkdirs()
            batchStatsOut.withWriter { w ->
                w.writeLine(['genome_id', 'out_path', 'n_claims', 'n_integrated', 'n_provenance',
                             'max_iter', 'enforce_stage2', 'enforce_mode',
                             'runtime_sec', 'used_heap_mb', 'total_heap_mb'].join('\t'))
                for (Map s : stats) {
                    w.writeLine([
                        s.genome_id,
                        s.out_path,
                        s.n_claims,
                        s.n_integrated,
                        s.n_provenance,
                        s.max_iter,
                        s.enforce_stage2,
                        s.enforce_mode,
                        String.format(Locale.ROOT, '%.4f', s.runtime_sec as double),
                        String.format(Locale.ROOT, '%.1f', s.used_heap_mb as double),
                        String.format(Locale.ROOT, '%.1f', s.total_heap_mb as double)
                    ].join('\t'))
                }
            }
            println "  Wrote batch runtime/memory stats: ${batchStatsOut}"
        }
    }
}
