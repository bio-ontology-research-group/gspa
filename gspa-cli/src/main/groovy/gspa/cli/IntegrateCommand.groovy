package gspa.cli

import com.fasterxml.jackson.databind.ObjectMapper
import gspa.integration.CalibrationTable
import gspa.integration.ClaimExtractor
import gspa.integration.EvidenceClaim
import gspa.integration.EvidenceCombiner
import gspa.integration.EvidenceType
import gspa.integration.IntegratedAnnotationSet
import gspa.integration.IntegrationState
import gspa.integration.IterativeRefiner
import gspa.integration.OuterIterativeRefiner
import gspa.integration.PriorEngine
import gspa.integration.promotion.AllAboveThresholdStrategy
import gspa.integration.promotion.BeamSearchStrategy
import gspa.integration.promotion.GreedyStrategy
import gspa.integration.promotion.MaxSatStrategy
import gspa.integration.promotion.PromotionStrategy
import gspa.integration.prior.AnatomicalPrior
import gspa.integration.prior.ConsistencyPrior
import gspa.integration.prior.CoherencePrior
import gspa.integration.prior.EssentialityPrior
import gspa.integration.prior.GapFillingPrior
import gspa.integration.prior.GenomicContextPrior
import gspa.integration.crossgenome.CrossGenomeReScorer
import gspa.integration.crossgenome.ReactionLocusCatalog
import gspa.integration.suggester.DarkMatterSuggester
import gspa.integration.suggester.DisjunctiveSuggestion
import gspa.integration.suggester.ReactionLocalContextSuggester
import gspa.integration.suggester.SingletonSuggestion
import gspa.integration.suggester.Suggestion
import gspa.metrics.CoherenceEnforcer
import gspa.metrics.CompletenessEnforcer
import gspa.metrics.ConsistencyEnforcer
import gspa.metrics.EnforcementAction
import gspa.metrics.EnforcementReport
import gspa.metrics.TaxonInference
import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.model.Contig
import gspa.model.Genome
import gspa.model.GenomeLayout
import gspa.model.GenomeLayoutLoader
import gspa.model.Protein
import gspa.ontology.AnatomicalConstraints
import gspa.ontology.ReactionGraph
import gspa.ontology.ReactionGraphLoader
import picocli.CommandLine.Command
import picocli.CommandLine.Option

/**
 * {@code gspa integrate} — run the Phase 7 evidence integrator on a
 * pre-parsed claims file with a theta.json of hyperparameters.
 *
 * Designed to be the per-iteration evaluator invoked by the Phase 7.4
 * benchmark optimizer. Runs without the full predictor stack: just reads
 * claims.jsonl, optionally wires ontology handles, runs IterativeRefiner,
 * writes an integrated TSV.
 */
@Command(
    name = 'integrate',
    description = 'Run the Phase 7 evidence integrator on pre-parsed claims.',
    mixinStandardHelpOptions = true
)
class IntegrateCommand implements Runnable {

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

    @Option(names = ['--enforce-stage2'],
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

    @Option(names = ['--lite'], description = 'Skip ELK initialization (no process coherence).')
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
            description = 'Gapseq search target: genome (default) | proteome | reps (requires clustering).')
    String gapseqTarget = 'genome'

    @Option(names = ['--intragenome-cluster'],
            description = 'Intragenome protein clustering: off (default), 0.9, 0.95, or 1.0.')
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

    @Override
    void run() {
        validatePhase10Flags()
        if (batchManifestFile != null) {
            runBatchManifest()
            return
        }
        if (claimsFile == null || outFile == null) {
            throw new IllegalArgumentException("Either --batch-manifest or both --claims and --out must be provided.")
        }
        println "GSPA integrate"
        println "  Claims: ${claimsFile}"
        println "  Theta:  ${thetaFile ?: '(defaults)'}"
        println "  Out:    ${outFile}"

        // --- Load claims ---
        def calibration = new CalibrationTable()
        def extractor = new ClaimExtractor(calibration)
        List<EvidenceClaim> claims = extractor.readClaimsJsonl(claimsFile)
        println "  Loaded ${claims.size()} claims"

        // --- Load theta.json ---
        Map theta = thetaFile?.exists() ? new ObjectMapper().readValue(thetaFile, Map) : [:]

        // --- Build combiner with learned reliability ---
        Map<EvidenceType, Double> reliability = EvidenceCombiner.defaultReliability()
        if (theta.reliability instanceof Map) {
            (theta.reliability as Map).each { k, v ->
                try {
                    reliability[EvidenceType.valueOf(k as String)] = (v as Number).doubleValue()
                } catch (ignored) { /* unknown type */ }
            }
        }
        def combiner = new EvidenceCombiner(reliability)
        if (theta.second_opinion_bonus != null) {
            combiner.secondOpinionBonus = (theta.second_opinion_bonus as Number).doubleValue()
        }

        // --- Build state with optional ontology + pathway + operon + gap handles ---
        def state = new IntegrationState(new Genome(id: 'benchmark'))
        wireReferenceData(state, claims)

        // --- Build prior engine from enabled list + theta strengths ---
        Set<String> enabled = (enabledPriors ?: '').split(',').collect { it.trim() } as Set
        def engine = buildPriorEngine(theta, enabled)

        // --- Build refiner ---
        def refiner = new IterativeRefiner(combiner)
        refiner.priorEngine = engine
        if (theta.max_iter != null)  refiner.maxIter = (theta.max_iter as Number).intValue()
        if (theta.epsilon != null)   refiner.epsilon = (theta.epsilon as Number).doubleValue()
        if (theta.damping != null)   refiner.damping = (theta.damping as Number).doubleValue()

        // --- Refine ---
        IntegratedAnnotationSet integrated = refiner.refine(claims, state)
        println "  Produced ${integrated.annotations.size()} integrated annotations"

        // --- Phase 12: optional Reaction-Local Context suggester (M1) ---
        if (rlcSuggester) {
            def rlc = new ReactionLocalContextSuggester()
            rlc.radiusK = rlcRadiusK
            rlc.alpha = rlcAlpha
            rlc.kernelBandwidth = rlcKernelBandwidth
            rlc.anchorPosteriorThreshold = rlcAnchorThreshold
            if (featuresOut != null) rlc.featuresOut = featuresOut
            rlc.suggest(state, integrated)
            println "  RLGC suggester emitted ${integrated.suggestions.size()} suggestions"

            // --- Phase 12 M2: cross-genome LR re-scorer ---
            if (rxnLocusCatalogFile != null && rxnLocusCatalogFile.exists()) {
                def catalog = ReactionLocusCatalog.readFrom(rxnLocusCatalogFile)
                println "  ReactionLocusCatalog: ${catalog.size()} (C, R) entries; panel=${catalog.panelSize}"
                def cgr = new CrossGenomeReScorer()
                cgr.lambda = cgLambda
                cgr.minSupport = cgMinSupport
                cgr.requireCredible = cgRequireCredible
                cgr.rescore(state, integrated, catalog)
                println "  Cross-genome rescoring complete; ${integrated.suggestions.size()} suggestions after rescore"
            }

            if (suggestionsOut != null && !darkMatter) {
                suggestionsOut.parentFile?.mkdirs()
                writeSuggestionsTsv(integrated.suggestions, suggestionsOut)
                println "  Suggestions: ${suggestionsOut}"
            }
        }

        // --- Phase 8: optional dark-matter suggester ---
        if (darkMatter) {
            def suggester = buildSuggester(theta)
            suggester.useRefinedBayesFactor = refinedBf

            if (iterateGapseq) {
                // Phase 10 outer fixed-point loop over (refine → suggest → promote → pin).
                println "  Running Phase 10 outer loop (maxIter=${maxGapseqIter}, qBase=${gapseqQBase}, qStep=${gapseqQStep}, pin=${gapseqPinPromotions}, strategy=${promotionStrategy}, refinedBf=${refinedBf})"
                def outer = new OuterIterativeRefiner(refiner)
                outer.suggester = suggester
                outer.maxIter = maxGapseqIter
                outer.qBase = gapseqQBase
                outer.qStep = gapseqQStep
                outer.qCap = gapseqQCap
                outer.pinPromotions = gapseqPinPromotions
                outer.promotionStrategy = buildPromotionStrategy()
                def gs = new OuterIterativeRefiner.CoverageGapSource(tauCover: gapseqTauCover, pathwayDb: state.pathwayDatabase)
                outer.gapSource = gs
                def outerResult = outer.refine(claims, state)
                integrated = outerResult.integrated
                println "  Outer loop: iter=${outerResult.outerIterationsRun}, fixedPoint=${outerResult.fixedPointReached}, cascade=${outerResult.cascadeRolledBack}, promoted_per_iter=${outerResult.promotedPerIter}, gaps_per_iter=${outerResult.gapsPerIter}"
            } else {
                suggester.suggest(state, integrated)
                println "  Dark-matter suggester emitted ${integrated.suggestions.size()} suggestions"
            }

            if (suggestionsOut != null) {
                suggestionsOut.parentFile?.mkdirs()
                writeSuggestionsTsv(integrated.suggestions, suggestionsOut)
                println "  Suggestions: ${suggestionsOut}"
            }
        }

        // --- Optional Stage 2A: Symbolic Constraint & Coherence Enforcement ---
        if (enforceStage2) {
            runStage2Enforcement(state, integrated)
        }

        // --- Optional Stage 2B: StoichiometricLocusRefiner ---
        if (enableFbaRefiner && fbaDataDir != null) {
            def fbaRefiner = new gspa.integration.fba.StoichiometricLocusRefiner(
                fbaDataDir: fbaDataDir,
                sidecarScript: fbaSidecar,
                pythonBin: fbaPythonBin,
                gram: fbaGram,
                medium: fbaMedium
            )
            def fbaSum = fbaRefiner.refine(state, integrated)
            println "  Stage 2B StoichiometricLocusRefiner (gram=${fbaGram}, medium=${fbaMedium}): coreRxns=${fbaSum.nCoreRxns}, repairedRxns=${fbaSum.nRepairedRxns}, uncertified=${fbaSum.nUncertifiedAdditions}, updatedLoci=${fbaSum.nUpdatedLoci} (${String.format(Locale.ROOT, '%.2f', fbaSum.solveTimeSec)}s)"
        }

        // --- Optional Stage 3: Neuro-Symbolic Agentic Completion (gspa-agent) ---
        if (anomaliesOut != null && !enableAgent) {
            def extractorAnom = new gspa.agent.AnomalyExtractor()
            def anomalies = extractorAnom.extract(state, integrated)
            gspa.agent.AnomalyExtractor.writeJsonl(anomalies, anomaliesOut)
            println "  Extracted ${anomalies.size()} anomalies to ${anomaliesOut}"
        }
        if (enableAgent) {
            File anomFile = anomaliesOut ?: new File(outFile.parentFile ?: new File('.'), "${outFile.name}.anomalies.jsonl")
            File resFile = agentResolutionsOut ?: new File(outFile.parentFile ?: new File('.'), "${outFile.name}.agent_resolutions.jsonl")
            def runner = new gspa.agent.AgentCompletionRunner(
                sidecarScript: agentSidecar,
                agentModel: agentModel,
                apiBaseUrl: agentApiBase,
                guardrail: !agentNoGuardrail
            )
            def sum = runner.run(state, integrated, anomFile, resFile)
            println "  Stage 3 Agent (${agentModel}, guardrail=${!agentNoGuardrail}): anomalies=${sum.anomaliesExtracted}, proposed=${sum.resolutionsProposed}, promoted=${sum.promoted}, relaxed=${sum.relaxed}, demoted=${sum.demoted}, abstained=${sum.abstained}, rejected_by_guardrail=${sum.rejectedByGuardrail} (${String.format(Locale.ROOT, '%.2f', sum.elapsedSeconds)}s)"
        }

        // --- Write output TSV ---
        writeIntegratedOutputTsv(integrated, outFile)
        println "  Wrote: ${outFile}"

        // --- Optional provenance JSON ---
        if (provenanceFile != null) {
            def mapper = new ObjectMapper()
            def provList = integrated.provenance.values().collect { prov ->
                [
                    function_key: prov.functionKey,
                    protein_id: prov.proteinId,
                    likelihood_logodds: prov.likelihoodLogOdds,
                    final_logodds: prov.finalLogOdds,
                    final_probability: prov.finalProbability,
                    prior_contributions: prov.priorContributions,
                    n_supporting: prov.supportingClaims?.size() ?: 0,
                    supporting_sources: prov.supportingClaims ? prov.supportingClaims*.source : [],
                    convergence_iter: prov.convergenceIter,
                ]
            }
            provenanceFile.parentFile?.mkdirs()
            mapper.writerWithDefaultPrettyPrinter().writeValue(provenanceFile, provList)
            println "  Provenance: ${provenanceFile}"
        }
    }

    private static void writeIntegratedOutputTsv(IntegratedAnnotationSet integrated, File targetOutFile) {
        targetOutFile.parentFile?.mkdirs()
        targetOutFile.withWriter { w ->
            w.writeLine(['protein_id', 'type', 'function_id', 'go_aspect',
                         'posterior_prob', 'likelihood_logodds', 'final_logodds',
                         'n_supporting', 'priors_fired', 'convergence_iter'].join('\t'))
            integrated.provenance.values().each { prov ->
                def fkey = IntegrationState.splitFunctionKey(prov.functionKey)
                if (fkey == null) return
                def fired = prov.priorContributions.collect { k, v ->
                    "${k}:${String.format(Locale.ROOT, '%.3f', v)}"
                }.join(',')
                def firstClaim = prov.supportingClaims ? prov.supportingClaims.first() : null
                w.writeLine([
                    fkey[0],
                    fkey[1],
                    fkey[2],
                    firstClaim?.goAspect ?: '',
                    String.format(Locale.ROOT, '%.6f', prov.finalProbability),
                    String.format(Locale.ROOT, '%.4f', prov.likelihoodLogOdds),
                    String.format(Locale.ROOT, '%.4f', prov.finalLogOdds),
                    (prov.supportingClaims?.size() ?: 0).toString(),
                    fired,
                    prov.convergenceIter.toString(),
                ].join('\t'))
            }
        }
    }

    /**
     * Batch execution mode (--batch-manifest): loads GO ontology, anatomical/taxon
     * constraints, complex classifications, essential function presets, and pathway
     * databases ONCE into memory, then processes all manifest rows across a parallel
     * worker thread pool, recording per-genome wall-clock runtime and memory usage.
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
     */
    private void runBatchManifest() {
        long t0Batch = System.nanoTime()
        println "GSPA integrate --batch-manifest ${batchManifestFile} (threads=${batchThreads})"

        Map theta = thetaFile?.exists() ? new ObjectMapper().readValue(thetaFile, Map) : [:]
        Map<EvidenceType, Double> reliability = EvidenceCombiner.defaultReliability()
        if (theta.reliability instanceof Map) {
            (theta.reliability as Map).each { k, v ->
                try {
                    reliability[EvidenceType.valueOf(k as String)] = (v as Number).doubleValue()
                } catch (ignored) {}
            }
        }

        // 1. Load shared reference data ONCE
        def sharedState = new IntegrationState(new Genome(id: 'shared_reference'))
        wireReferenceData(sharedState, [])
        Map<String, gspa.config.EssentialFunctions> essentialPresets = [:]
        ['bacteria', 'archaea', 'eukaryote', 'virus'].each { prof ->
            try {
                essentialPresets[prof] = gspa.config.EssentialFunctions.loadPreset(prof)
            } catch (Exception ignored) {}
        }
        // Ensure base SAT consistency checker is loaded even if --taxon was not passed on CLI
        if (sharedState.satConsistencyChecker == null) {
            try {
                def tc = new gspa.ontology.TaxonConstraints()
                if (taxonConstraintsFile != null && taxonConstraintsFile.exists()) {
                    if (taxonConstraintsFile.name.toLowerCase().endsWith('.obo')) tc.loadFromObo(taxonConstraintsFile)
                    else tc.loadFromTsv(taxonConstraintsFile)
                } else {
                    File bundledTc = classpathToTempFile('/taxon-constraints/go-taxon-constraints.tsv', 'go-taxon-constraints', '.tsv')
                    if (bundledTc != null && bundledTc.exists()) tc.loadFromTsv(bundledTc)
                }
                boolean wantAnatShared = enablePriors == null || enablePriors.trim().toLowerCase() == 'all' ||
                    enablePriors.split(',').collect { it.trim().toLowerCase() }.contains('anatomical')
                if (sharedState.anatomicalConstraints != null && wantAnatShared) {
                    sharedState.anatomicalConstraints.augmentTaxonConstraints(tc)
                }
                def checker = new gspa.ontology.SatConsistencyChecker(tc)
                if (taxonomyFile != null && taxonomyFile.exists()) {
                    checker.loadTaxonomyTsv(taxonomyFile)
                } else {
                    File bundledHier = classpathToTempFile('/taxon-constraints/ncbi-taxon-hierarchy.tsv', 'ncbi-taxon-hierarchy', '.tsv')
                    if (bundledHier != null && bundledHier.exists()) checker.loadTaxonomyTsv(bundledHier)
                }
                sharedState.satConsistencyChecker = checker
            } catch (Exception e) {
                System.err.println "  [warn] Failed to initialize shared SAT checker: ${e.message}"
            }
        }
        double initSec = (System.nanoTime() - t0Batch) / 1e9d
        println String.format(Locale.ROOT, "  Shared reference initialization completed in %.2f s", initSec)

        // 2. Parse manifest rows
        List<String[]> rows = []
        batchManifestFile.eachLine { line, idx ->
            String trimmed = line.trim()
            if (trimmed.isEmpty() || trimmed.startsWith('#')) return
            String[] parts = trimmed.split('\t', -1)
            if (idx == 1 && (parts[0] in ['genome_id', 'tag', 'id'] || parts[1] in ['claims_path', 'claims'])) return
            if (parts.length >= 3) {
                rows << parts
            }
        }
        println "  Loaded ${rows.size()} batch tasks from ${batchManifestFile}"

        // 3. Process rows across thread pool
        int nThreads = Math.max(1, batchThreads)
        def pool = java.util.concurrent.Executors.newFixedThreadPool(nThreads)
        List<java.util.concurrent.Future<Map>> futures = []

        for (int rowIdx = 0; rowIdx < rows.size(); rowIdx++) {
            final String[] r = rows[rowIdx]
            futures << pool.submit({ ->
                long tStart = System.nanoTime()
                String tag = r[0].trim()
                File rowClaimsFile = new File(r[1].trim())
                File rowOutFile = new File(r[2].trim())
                String rowTaxonRaw = (r.length > 3 && r[3].trim() && r[3].trim() != '-') ? r[3].trim() : organismTaxon
                String rowEssProf = (r.length > 4 && r[4].trim() && r[4].trim() != '-') ? r[4].trim() : essentialProfile
                File rowOperonsFile = (r.length > 5 && r[5].trim() && r[5].trim() != '-') ? new File(r[5].trim()) : operonsFile
                File rowGapsFile = (r.length > 6 && r[6].trim() && r[6].trim() != '-') ? new File(r[6].trim()) : gapsFile
                Integer rowMaxIter = (r.length > 7 && r[7].trim() && r[7].trim() != '-') ? Integer.parseInt(r[7].trim()) : (theta.max_iter != null ? (theta.max_iter as Number).intValue() : null)
                boolean rowEnforceS2 = (r.length > 8 && r[8].trim() && r[8].trim() != '-') ? Boolean.parseBoolean(r[8].trim()) : enforceStage2
                String rowEnforceMode = (r.length > 9 && r[9].trim() && r[9].trim() != '-') ? r[9].trim() : enforceConsistencyMode
                String rowPriorsStr = (r.length > 10 && r[10].trim() && r[10].trim() != '-') ? r[10].trim() : enabledPriors

                def calibration = new CalibrationTable()
                def extractor = new ClaimExtractor(calibration)
                List<EvidenceClaim> claims = extractor.readClaimsJsonl(rowClaimsFile)

                def state = new IntegrationState(new Genome(id: tag))
                state.goOntology = sharedState.goOntology
                state.goReasoner = sharedState.goReasoner
                state.anatomicalConstraints = sharedState.anatomicalConstraints
                state.complexClassification = sharedState.complexClassification
                state.pathwayDatabase = sharedState.pathwayDatabase
                state.reactionGraph = sharedState.reactionGraph
                if (rowEssProf != null && essentialPresets.containsKey(rowEssProf.toLowerCase())) {
                    state.essentialFunctions = essentialPresets[rowEssProf.toLowerCase()]
                } else {
                    state.essentialFunctions = sharedState.essentialFunctions
                }
                String resolvedTax = resolveOrganismTaxon(rowTaxonRaw)
                if (sharedState.satConsistencyChecker != null && (resolvedTax != null || rowEnforceS2)) {
                    state.satConsistencyChecker = sharedState.satConsistencyChecker.copyForTaxon(resolvedTax)
                }
                if (rowOperonsFile != null && rowOperonsFile.exists()) {
                    List<List<String>> ops = []
                    rowOperonsFile.eachLine { ln ->
                        ln = ln.trim()
                        if (!ln.isEmpty() && !ln.startsWith('#')) ops << (ln.split('\t') as List)
                    }
                    state.operons = ops
                }
                if (rowGapsFile != null && rowGapsFile.exists()) {
                    def mapper = new ObjectMapper()
                    List gaps = []
                    rowGapsFile.eachLine { ln ->
                        ln = ln.trim()
                        if (!ln.isEmpty() && !ln.startsWith('#')) {
                            Map rec = mapper.readValue(ln, Map)
                            gaps << new gspa.integration.MetabolicGap(
                                pathwayId: rec.pathway_id as String,
                                reactionId: rec.reaction_id as String,
                                ecNumber: rec.ec_number as String,
                                goTerm: rec.go_term as String,
                                gapseqGuessed: (rec.gapseq_guessed ?: false) as boolean,
                            )
                        }
                    }
                    state.metabolicGaps = gaps
                }

                def combiner = new EvidenceCombiner(reliability)
                if (theta.second_opinion_bonus != null) {
                    combiner.secondOpinionBonus = (theta.second_opinion_bonus as Number).doubleValue()
                }
                Set<String> enabled = (rowPriorsStr ?: '').split(',').collect { it.trim() }.findAll { it && it != 'none' } as Set
                def engine = buildPriorEngine(theta, enabled)
                def refiner = new IterativeRefiner(combiner)
                refiner.priorEngine = engine
                if (rowMaxIter != null) refiner.maxIter = rowMaxIter
                if (theta.epsilon != null) refiner.epsilon = (theta.epsilon as Number).doubleValue()
                if (theta.damping != null) refiner.damping = (theta.damping as Number).doubleValue()

                IntegratedAnnotationSet integrated = refiner.refine(claims, state)
                if (rowEnforceS2) {
                    runStage2Enforcement(state, integrated, rowEnforceMode, null, false)
                }
                writeIntegratedOutputTsv(integrated, rowOutFile)

                double elapsedSec = (System.nanoTime() - tStart) / 1e9d
                Runtime rt = Runtime.getRuntime()
                double usedHeapMb = (rt.totalMemory() - rt.freeMemory()) / (1024.0d * 1024.0d)
                double totalHeapMb = rt.totalMemory() / (1024.0d * 1024.0d)
                return [
                    genome_id: tag,
                    out_path: rowOutFile.path,
                    n_claims: claims.size(),
                    n_integrated: integrated.annotations.size(),
                    n_provenance: integrated.provenance.size(),
                    max_iter: refiner.maxIter,
                    enforce_stage2: rowEnforceS2,
                    enforce_mode: rowEnforceMode,
                    runtime_sec: elapsedSec,
                    used_heap_mb: usedHeapMb,
                    total_heap_mb: totalHeapMb
                ]
            } as java.util.concurrent.Callable<Map>)
        }

        List<Map> stats = []
        int done = 0
        for (def f : futures) {
            stats << f.get()
            done++
            if (done % 50 == 0 || done == rows.size()) {
                double elapsedAll = (System.nanoTime() - t0Batch) / 1e9d
                println String.format(Locale.ROOT, "  Batch progress: %d / %d completed (%.2f s elapsed)", done, rows.size(), elapsedAll)
            }
        }
        pool.shutdown()

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

    /**
     * Wire optional reference data into the IntegrationState. Each piece is
     * independent: if a file isn't provided, the corresponding prior
     * no-ops.
     */
    private void wireReferenceData(IntegrationState state, List<EvidenceClaim> claims = []) {
        // GO ontology + ELK reasoner
        if (goOwlFile != null) {
            try {
                def goOntology = new gspa.ontology.GoOntology()
                goOntology.loadOwl(goOwlFile)
                state.goOntology = goOntology
                if (!lite) {
                    def reasoner = new gspa.ontology.GoReasoner(goOntology)
                    if (reasonerCacheDir != null) {
                        reasoner.cacheDir = reasonerCacheDir
                    }
                    reasoner.initialize()
                    state.goReasoner = reasoner
                }
                println "  Loaded GO ontology from ${goOwlFile}"
            } catch (Exception e) {
                System.err.println "  [warn] Failed to load GO ontology: ${e.message}"
            }
        }

        // Bundled anatomical constraints & complex terms
        try {
            state.anatomicalConstraints = AnatomicalConstraints.loadBundled()
            println "  Anatomical constraints: ${state.anatomicalConstraints.constrainedGoTermCount()} GO terms, ${state.anatomicalConstraints.constrainedAnatomyCount()} anatomical structures"
        } catch (Exception e) {
            System.err.println "  [warn] Failed to load bundled anatomical constraints: ${e.message}"
        }

        try {
            if (complexTermsFile != null && complexTermsFile.exists()) {
                complexTermsFile.eachLine { line ->
                    if (line.startsWith('GO_term') || line.startsWith('#') || line.trim().isEmpty()) return
                    String[] f = line.split('\t')
                    if (f.length >= 2) state.complexClassification[f[0].trim()] = f[1].trim()
                }
            } else {
                state.loadBundledComplexTerms()
            }
            println "  Complex classification: ${state.complexClassification.size()} GO complex terms"
        } catch (Exception e) {
            System.err.println "  [warn] Failed to load complex classification: ${e.message}"
        }

        // Essential functions
        if (essentialFunctionsFile != null) {
            try {
                state.essentialFunctions = gspa.config.EssentialFunctions
                    .loadFromTsv(essentialFunctionsFile, essentialProfile)
                println "  Essential functions: ${state.essentialFunctions.goTerms.size()} GO terms"
            } catch (Exception e) {
                System.err.println "  [warn] Failed to load essential functions: ${e.message}"
            }
        } else if (essentialProfile != null) {
            try {
                state.essentialFunctions = gspa.config.EssentialFunctions.loadPreset(essentialProfile)
                println "  Essential functions (preset=${essentialProfile}): ${state.essentialFunctions.goTerms.size()} GO terms"
            } catch (Exception e) {
                System.err.println "  [warn] Failed to load essential preset '${essentialProfile}': ${e.message}"
            }
        }

        // Pathway database — supports stacking multiple sources (KEGG main +
        // KEGG Modules + MetaCyc + BioCyc) so a bacterial genome's operons
        // get reasonable enrichment hits at each granularity.
        if ((pathwaysFile != null && ec2goFile != null) || (moduleFiles != null && !moduleFiles.isEmpty())) {
            try {
                if (pathwaysFile != null && ec2goFile != null) {
                    state.pathwayDatabase = gspa.ontology.PathwayLoader.loadPathways(
                        pathwaysFile,
                        gspa.ontology.PathwayLoader.loadEc2Go(ec2goFile),
                    )
                    println "  Pathways: ${state.pathwayDatabase.pathways.size()} loaded from ${pathwaysFile.name}"
                } else {
                    state.pathwayDatabase = new gspa.ontology.PathwayDatabase()
                }
                for (File mod : (moduleFiles ?: [])) {
                    if (mod != null && mod.exists()) {
                        gspa.ontology.PathwayLoader.loadPathwaysInto(state.pathwayDatabase, mod)
                        println "  Pathways: now ${state.pathwayDatabase.pathways.size()} (added ${mod.name})"
                    }
                }
            } catch (Exception e) {
                System.err.println "  [warn] Failed to load pathway database: ${e.message}"
            }
        }

        // SAT consistency checker: wire when --taxonomy, --taxon, --infer-taxon, or --enforce-stage2 is given.
        boolean wantSat = (taxonomyFile != null) || (organismTaxon != null && !organismTaxon.trim().isEmpty()) || inferTaxon || enforceStage2
        if (wantSat) {
            try {
                def taxonConstraints = new gspa.ontology.TaxonConstraints()
                if (taxonConstraintsFile != null && taxonConstraintsFile.exists()) {
                    if (taxonConstraintsFile.name.toLowerCase().endsWith('.obo')) {
                        taxonConstraints.loadFromObo(taxonConstraintsFile)
                    } else {
                        taxonConstraints.loadFromTsv(taxonConstraintsFile)
                    }
                } else {
                    File bundledTc = classpathToTempFile('/taxon-constraints/go-taxon-constraints.tsv', 'go-taxon-constraints', '.tsv')
                    if (bundledTc != null && bundledTc.exists()) {
                        taxonConstraints.loadFromTsv(bundledTc)
                    } else if (state.goOntology != null) {
                        taxonConstraints.loadFromGoOntology(state.goOntology)
                    }
                }
                boolean wantAnat = enablePriors == null || enablePriors.trim().toLowerCase() == 'all' ||
                    enablePriors.split(',').collect { it.trim().toLowerCase() }.contains('anatomical')
                if (state.anatomicalConstraints != null && wantAnat) {
                    state.anatomicalConstraints.augmentTaxonConstraints(taxonConstraints)
                }
                if (taxonConstraints.constrainedTermCount() > 0) {
                    def checker = new gspa.ontology.SatConsistencyChecker(taxonConstraints)
                    if (taxonomyFile != null && taxonomyFile.exists()) {
                        // Detect 3-column TSV (Term, Relationship, Other) vs 2-column legacy
                        String firstData = taxonomyFile.readLines().find { !it.trim().isEmpty() && !it.startsWith('#') } ?: ''
                        if (firstData.split('\t').length >= 3) {
                            checker.loadTaxonomyTsv(taxonomyFile)
                        } else {
                            checker.loadTaxonomyHierarchy(taxonomyFile)
                        }
                    } else {
                        File bundledHier = classpathToTempFile('/taxon-constraints/ncbi-taxon-hierarchy.tsv', 'ncbi-taxon-hierarchy', '.tsv')
                        if (bundledHier != null && bundledHier.exists()) {
                            checker.loadTaxonomyTsv(bundledHier)
                        }
                    }
                    String resolvedTax = resolveOrganismTaxon(organismTaxon)
                    if (resolvedTax) {
                        checker.organismTaxon = resolvedTax
                        println "  Asserted organism taxon: ${resolvedTax}"
                    } else if (inferTaxon && claims) {
                        Map<String, Double> termScore = [:].withDefault { 0.0d }
                        for (EvidenceClaim c : claims) {
                            if (c.functionType == 'GO' && c.calibratedProb > termScore[c.functionId]) {
                                termScore[c.functionId] = c.calibratedProb
                            }
                        }
                        def inf = new TaxonInference(checker: checker, goOntology: state.goOntology, minScore: 0.5d).infer(termScore)
                        if (inf?.taxon) {
                            checker.organismTaxon = inf.taxon
                            println "  Inferred organism taxon: ${inf.taxon} (${inf.label})"
                        }
                    }
                    state.satConsistencyChecker = checker
                    println "  SAT checker: ${taxonConstraints.onlyInTaxon.size()} only_in + ${taxonConstraints.neverInTaxon.size()} never_in (organismTaxon=${checker.organismTaxon ?: 'none'})"
                }
            } catch (Exception e) {
                System.err.println "  [warn] Failed to wire SAT consistency checker: ${e.message}"
            }
        } else {
            println "  SAT checker: disabled (pass --taxon, --infer-taxon, or --taxonomy to enable)"
        }

        // Operons TSV: one line per operon, tab-separated protein IDs.
        if (operonsFile != null && operonsFile.exists()) {
            List<List<String>> operons = []
            operonsFile.eachLine { line ->
                line = line.trim()
                if (line.isEmpty() || line.startsWith('#')) return
                operons << (line.split('\t') as List)
            }
            state.operons = operons
            println "  Operons: ${operons.size()}"
        }

        // Metabolic gaps JSONL
        if (gapsFile != null && gapsFile.exists()) {
            def mapper = new ObjectMapper()
            List gaps = []
            gapsFile.eachLine { line ->
                line = line.trim()
                if (line.isEmpty() || line.startsWith('#')) return
                Map rec = mapper.readValue(line, Map)
                gaps << new gspa.integration.MetabolicGap(
                    pathwayId: rec.pathway_id as String,
                    reactionId: rec.reaction_id as String,
                    ecNumber: rec.ec_number as String,
                    goTerm: rec.go_term as String,
                    gapseqGuessed: (rec.gapseq_guessed ?: false) as boolean,
                )
            }
            state.metabolicGaps = gaps
            println "  Metabolic gaps: ${gaps.size()}"
        }

        // Phase 11 cross-genome: orthogroup membership + consensus.
        if (orthogroupsFile != null && orthogroupsFile.exists()) {
            Map<String, String> orthogroups = new LinkedHashMap<>()
            orthogroupsFile.eachLine { line ->
                line = line.trim()
                if (line.isEmpty() || line.startsWith('#')) return
                String[] parts = line.split('\t')
                if (parts.length >= 2) {
                    // Strip optional "namespace:" prefix from the rep so the
                    // orthogroup ID matches the bare-accession format used by
                    // build_catalog.py and ReactionLocusCatalog.
                    String rep = parts[1]
                    int colon = rep.indexOf(':')
                    if (colon >= 0) rep = rep.substring(colon + 1)
                    orthogroups[parts[0]] = rep
                }
            }
            state.orthogroupMap = orthogroups
            println "  Orthogroups: ${orthogroups.size()} protein assignments"
        }
        if (clusterConsensusFile != null && clusterConsensusFile.exists()) {
            Map<String, Double> consensus = new LinkedHashMap<>()
            clusterConsensusFile.eachLine { line ->
                line = line.trim()
                if (line.isEmpty() || line.startsWith('#')) return
                String[] parts = line.split('\t')
                if (parts.length >= 4) {
                    String key = "${parts[0]}|${parts[1]}|${parts[2]}".toString()
                    try {
                        consensus[key] = Double.parseDouble(parts[3])
                    } catch (NumberFormatException ignored) {}
                }
            }
            state.orthogroupConsensus = consensus
            println "  Cluster consensus entries: ${consensus.size()}"
        }

        // Phase 12 RLGC: reaction graph + genome layout.
        if (reactionGraphFile != null && reactionGraphFile.exists()) {
            try {
                state.reactionGraph = ReactionGraphLoader.load(
                    reactionGraphFile, diffusionMetsFile, rlcCurrencyPercentile,
                    reactionEcAliasesFile)
                println "  Reaction graph: ${state.reactionGraph.reactions.size()} reactions, ${state.reactionGraph.currencyMetabolites.size()} currency metabolites, ${state.reactionGraph.ecToReactions.size()} ECs"
            } catch (Exception e) {
                System.err.println "  [warn] Failed to load reaction graph: ${e.message}"
            }
        }
        if (genomeLayoutFile != null && genomeLayoutFile.exists()) {
            try {
                state.genomeLayout = GenomeLayoutLoader.load(genomeLayoutFile)
                println "  Genome layout: ${state.genomeLayout.size()} loci across ${state.genomeLayout.byContig.size()} contig(s)"
            } catch (Exception e) {
                System.err.println "  [warn] Failed to load genome layout: ${e.message}"
            }
        }
    }

    /**
     * Run Stage 2 Symbolic Constraint & Coherence Enforcement on the Stage 1
     * IntegratedAnnotationSet, updating both {@code integrated.annotations} and
     * {@code integrated.provenance} so downstream CAFA and GAEF evaluators
     * observe the enforced/relaxed/imputed annotations.
     */
    private void runStage2Enforcement(IntegrationState state, IntegratedAnnotationSet integrated,
                                      String modeOverride = null, File actionsOutOverride = null,
                                      boolean verbose = true) {
        if (state.satConsistencyChecker == null) {
            if (verbose) println "  [warn] Stage 2 enforcement skipped: no SAT consistency checker available"
            return
        }
        String effMode = modeOverride ?: (enforceConsistencyMode ?: 'relax')
        File effActionsOut = actionsOutOverride != null ? actionsOutOverride : (verbose ? enforcementActionsOut : null)
        Genome genome = new Genome(id: 'stage2')
        Contig contig = new Contig(id: 'contig_1', genome: genome)
        genome.contigs << contig
        Map<String, Protein> byProtein = new LinkedHashMap<>()

        // Populate genome with all GO claims in integrated.provenance
        integrated.provenance.values().each { prov ->
            String[] fkey = IntegrationState.splitFunctionKey(prov.functionKey)
            if (fkey == null || fkey[1] != 'GO') return
            String pid = fkey[0]
            String goId = fkey[2]
            Protein p = byProtein.computeIfAbsent(pid) {
                def np = new Protein(id: pid, contig: contig)
                contig.proteins << np
                np
            }
            def firstClaim = prov.supportingClaims ? prov.supportingClaims.first() : null
            p.annotations.add(new Annotation(
                type: AnnotationType.GO,
                value: goId,
                score: prov.finalProbability,
                source: firstClaim?.source ?: 'integrator',
                evidence: 'IEA',
                goAspect: firstClaim?.goAspect
            ))
        }

        Set<String> enabledStage2Set = (enablePriors ?: '').split(',').collect { it.trim() }.findAll { it && it != 'none' } as Set
        boolean useAnatomicalInStage2 = enabledStage2Set.isEmpty() || enabledStage2Set.contains('anatomical')
        def report = new EnforcementReport(provenance: true)
        def consEnf = new ConsistencyEnforcer(
            goOntology: state.goOntology,
            satChecker: state.satConsistencyChecker,
            anatomicalConstraints: useAnatomicalInStage2 ? state.anatomicalConstraints : null,
            mode: effMode,
            report: report
        )
        def consRes = consEnf.enforce(genome)
        if (verbose) println "  Stage 2 Consistency (${effMode}): ${consRes}"

        // Rebuild GO entries in integrated.provenance directly from the enforced genome proteins
        // so removals, downranks, and ancestor relaxations are 100% reflected.
        int relaxedAdded = 0
        Set<String> goKeysToRemove = integrated.provenance.keySet().findAll { k ->
            String[] fkey = IntegrationState.splitFunctionKey(k)
            fkey != null && fkey[1] == 'GO'
        } as Set
        Map<String, gspa.integration.ClaimProvenance> oldGoProv = [:]
        for (String k : goKeysToRemove) {
            oldGoProv[k] = integrated.provenance.remove(k)
        }
        integrated.annotations.annotations.removeAll { it.type == AnnotationType.GO }

        for (Protein p : contig.proteins) {
            for (Annotation a : p.annotations.goAnnotations()) {
                String fkey = IntegrationState.splitFunctionKey(IntegrationState.functionKey(p.id, 'GO', a.value)) ? IntegrationState.functionKey(p.id, 'GO', a.value) : "${p.id}|GO|${a.value}"
                String fromGo = null
                if (a.evidence?.contains('relaxed-from:')) {
                    def m = (a.evidence =~ /relaxed-from:(GO:\d+)/)
                    if (m.find()) fromGo = m.group(1)
                }
                def origProv = fromGo ? oldGoProv[IntegrationState.functionKey(p.id, 'GO', fromGo)] : null
                if (origProv != null) {
                    double validPriorSum = (origProv.priorContributions?.findAll { k, v ->
                        !(k in ['consistency', 'anatomical']) || (v != null && v > 0.0d)
                    }?.values()?.sum() ?: 0.0d) as double
                    double unpenalizedLo = origProv.likelihoodLogOdds + validPriorSum
                    double unpenalizedProb = CalibrationTable.logOddsToProb(unpenalizedLo) * consEnf.relaxFactor
                    if (unpenalizedProb > a.score) {
                        a.score = unpenalizedProb
                    }
                }
                def existing = oldGoProv[fkey]
                if (existing != null) {
                    existing.finalProbability = Math.max(existing.finalProbability, a.score)
                    existing.finalLogOdds = CalibrationTable.probToLogOdds(existing.finalProbability)
                    a.score = existing.finalProbability
                    if (a.evidence?.contains('relaxed-from:')) {
                        existing.priorContributions['stage2_relax'] = 1.0d
                    }
                    integrated.provenance[fkey] = existing
                } else {
                    def dummyClaim = new EvidenceClaim(
                        proteinId: p.id,
                        functionId: a.value,
                        functionType: AnnotationType.GO,
                        goAspect: a.goAspect ?: (origProv?.supportingClaims ? origProv.supportingClaims.first().goAspect : (state.goOntology?.getAspect(a.value) ?: '')),
                        source: a.source ?: 'stage2_relax',
                        evidenceType: EvidenceType.SEQUENCE_SIMILARITY,
                        rawScore: a.score,
                        calibratedProb: a.score
                    )
                    double lo = CalibrationTable.probToLogOdds(a.score)
                    Map<String, Double> priors = [:]
                    if (origProv?.priorContributions) {
                        origProv.priorContributions.each { k, v ->
                            if (!(k in ['consistency', 'anatomical']) || (v != null && v > 0.0d)) {
                                priors[k] = v
                            }
                        }
                    }
                    priors['stage2_relax'] = 1.0d
                    integrated.provenance[fkey] = new gspa.integration.ClaimProvenance(
                        functionKey: fkey,
                        proteinId: p.id,
                        supportingClaims: origProv?.supportingClaims ?: [dummyClaim],
                        likelihoodLogOdds: origProv ? origProv.likelihoodLogOdds : lo,
                        priorContributions: priors,
                        finalLogOdds: lo,
                        finalProbability: a.score,
                        convergenceIter: origProv ? origProv.convergenceIter : 1
                    )
                    relaxedAdded++
                }
                if (a.score >= 0.5d) {
                    integrated.annotations.add(a)
                }
            }
        }
        // Obligate heteromeric complex singleton pruning: if an obligate heteromer ('n') has only a single
        // active protein (>= 0.50) and no partner candidate (>= 0.08) across the proteome, attenuate it to 0.45.
        if (state.complexClassification != null && !state.complexClassification.isEmpty()) {
            Map<String, List<String>> activeHeteromerByTerm = [:].withDefault { [] }
            Map<String, List<String>> candHeteromerByTerm = [:].withDefault { [] }
            for (def prov : integrated.provenance.values()) {
                String[] fkey = IntegrationState.splitFunctionKey(prov.functionKey)
                if (fkey == null || fkey[1] != 'GO') continue
                String goId = fkey[2]
                if (state.complexClassification[goId] != 'n') continue
                if (prov.finalProbability >= 0.5d) {
                    activeHeteromerByTerm[goId] << prov.functionKey
                } else if (prov.finalProbability >= 0.08d) {
                    candHeteromerByTerm[goId] << prov.functionKey
                }
            }
            for (def entry : activeHeteromerByTerm.entrySet()) {
                String goId = entry.key
                List<String> activeKeys = entry.value
                List<String> candKeys = candHeteromerByTerm[goId] ?: []
                if (activeKeys.size() == 1 && candKeys.isEmpty()) {
                    def singleProv = integrated.provenance[activeKeys[0]]
                    if (singleProv != null && singleProv.finalProbability > 0.45d) {
                        double beforeProb = singleProv.finalProbability
                        singleProv.finalProbability = 0.45d
                        singleProv.finalLogOdds = CalibrationTable.probToLogOdds(0.45d)
                        singleProv.priorContributions['stage2_singleton_prune'] = -0.5d
                        report.record(new EnforcementAction(
                            dimension: 'coherence',
                            action: 'PRUNED_SINGLETON',
                            contigId: 'contig_1',
                            proteinId: singleProv.proteinId,
                            term: goId,
                            reason: 'obligate heteromer singleton without partner candidate (>=0.08)',
                            basis: 'complexClassification=n',
                            scoreBefore: beforeProb,
                            scoreAfter: 0.45d
                        ))
                        integrated.annotations.annotations.removeAll {
                            it.type == AnnotationType.GO && it.value == goId && singleProv.proteinId == (it.metaClass.respondsTo(it, 'getProteinId') ? it.proteinId : singleProv.proteinId)
                        }
                    }
                }
            }
        }
        if (verbose) println "  Stage 2 synced to integrated output (relaxed_new=${relaxedAdded})"

        if (effActionsOut != null && report.count() > 0) {
            if (state.goOntology != null) {
                report.actions.each { it.termLabel = state.goOntology.getLabel(it.term) }
            }
            effActionsOut.parentFile?.mkdirs()
            report.writeTsv(effActionsOut)
            if (verbose) println "  Stage 2 enforcement actions: ${effActionsOut}"
        }
    }

    private static String resolveOrganismTaxon(String raw) {
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
        return raw
    }

    private File classpathToTempFile(String resourcePath, String prefix, String suffix) {
        def stream = getClass().getResourceAsStream(resourcePath)
        if (stream == null) return null
        try {
            File tmp = File.createTempFile(prefix, suffix)
            tmp.deleteOnExit()
            tmp.withOutputStream { out -> stream.transferTo(out) }
            return tmp
        } finally {
            stream.close()
        }
    }

    /**
     * Write suggestions to TSV. Columns:
     *   kind, pathway_id, reaction_id, function_id, go_aspect, operon_id,
     *   bayes_factor, suggestion_score, protein_ids, q_values, provenance
     */
    private static void writeSuggestionsTsv(List<Suggestion> suggestions, File out) {
        out.withWriter { w ->
            w.writeLine([
                'kind', 'pathway_id', 'reaction_id', 'function_id', 'go_aspect',
                'operon_id', 'bayes_factor', 'suggestion_score',
                'protein_ids', 'q_values', 'provenance'
            ].join('\t'))
            for (Suggestion s : suggestions) {
                String proteinList
                String qList
                if (s instanceof SingletonSuggestion) {
                    proteinList = s.proteinId
                    qList = String.format(Locale.ROOT, '%.4f', s.q)
                } else if (s instanceof DisjunctiveSuggestion) {
                    proteinList = s.proteinIds.join(',')
                    qList = s.qValues.collect {
                        String.format(Locale.ROOT, '%.4f', it)
                    }.join(',')
                } else {
                    proteinList = ''
                    qList = ''
                }
                w.writeLine([
                    s.kind(),
                    s.pathwayId ?: '',
                    s.reactionId ?: '',
                    s.functionId ?: '',
                    s.goAspect ?: '',
                    s.operonId ?: '',
                    String.format(Locale.ROOT, '%.2f', s.bayesFactor),
                    String.format(Locale.ROOT, '%.4f', s.suggestionScore),
                    proteinList,
                    qList,
                    s.provenance ?: '',
                ].join('\t'))
            }
        }
    }

    /** Validate cross-flag constraints for Phase 10 options. */
    private void validatePhase10Flags() {
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
        if (promotionStrategy != null && promotionStrategy != 'default' &&
                promotionStrategy != 'greedy' && promotionStrategy != 'maxsat' &&
                promotionStrategy != 'beam' && promotionStrategy != '') {
            throw new IllegalArgumentException(
                "--promotion-strategy must be one of: default, greedy, maxsat, beam (got '${promotionStrategy}')")
        }
        if (beamWidth < 1) throw new IllegalArgumentException("--beam-width must be ≥ 1")
        if (beamCandidatesPerGap < 1) throw new IllegalArgumentException("--beam-candidates-per-gap must be ≥ 1")
    }

    private PromotionStrategy buildPromotionStrategy() {
        switch (promotionStrategy) {
            case 'default':
            case null:
            case '':
                return new AllAboveThresholdStrategy()
            case 'greedy':
                return new GreedyStrategy()
            case 'maxsat':
                return new MaxSatStrategy(coherenceBonusWeight: maxsatCoherenceBonus)
            case 'beam':
                return new BeamSearchStrategy(beamWidth: beamWidth, candidatesPerGap: beamCandidatesPerGap)
            default:
                throw new IllegalArgumentException(
                    "Unknown --promotion-strategy: '${promotionStrategy}' (expected default|greedy|maxsat|beam)")
        }
    }

    private DarkMatterSuggester buildSuggester(Map theta) {
        def suggester = new DarkMatterSuggester()
        if (theta.dark_matter instanceof Map) {
            def dm = theta.dark_matter as Map
            if (dm.bf_min != null)
                suggester.bfMin = (dm.bf_min as Number).doubleValue()
            if (dm.gamma_in_p != null)
                suggester.gammaInP = (dm.gamma_in_p as Number).doubleValue()
            if (dm.coverage_threshold != null)
                suggester.coverageThreshold = (dm.coverage_threshold as Number).doubleValue()
        }
        suggester
    }

    private PriorEngine buildPriorEngine(Map theta, Set<String> enabled) {
        def engine = new PriorEngine()
        Map priorWeights = (theta.prior_weights ?: [:]) as Map

        if ('essentiality' in enabled) {
            def p = new EssentialityPrior()
            if (theta.alpha_ess != null) p.alphaEss = (theta.alpha_ess as Number).doubleValue()
            engine.register(p, (priorWeights.essentiality ?: 1.0) as double)
        }
        if ('coherence' in enabled) {
            def p = new CoherencePrior()
            if (theta.alpha_coh != null) p.alphaCoh = (theta.alpha_coh as Number).doubleValue()
            if (theta.alpha_complex_boost != null) p.complexPartnerWeight = (theta.alpha_complex_boost as Number).doubleValue()
            if (theta.alpha_complex_penalty != null) p.complexSingletonPenaltyWeight = (theta.alpha_complex_penalty as Number).doubleValue()
            engine.register(p, (priorWeights.coherence ?: 1.0) as double)
        }
        if ('consistency' in enabled) {
            def p = new ConsistencyPrior()
            if (theta.alpha_cons != null) p.alphaCons = (theta.alpha_cons as Number).doubleValue()
            if (theta.consistency_hard_filter != null) {
                p.hardFilter = theta.consistency_hard_filter as boolean
            }
            engine.register(p, (priorWeights.consistency ?: 1.0) as double)
        }
        if ('anatomical' in enabled) {
            def p = new AnatomicalPrior()
            if (theta.alpha_anat_penalty != null) p.alphaAnatPenalty = (theta.alpha_anat_penalty as Number).doubleValue()
            if (theta.alpha_anat_boost != null) p.alphaAnatBoost = (theta.alpha_anat_boost as Number).doubleValue()
            if (theta.anatomical_hard_filter != null) p.hardFilter = theta.anatomical_hard_filter as boolean
            engine.register(p, (priorWeights.anatomical ?: 1.0) as double)
        }
        if ('gap_filling' in enabled) {
            def p = new GapFillingPrior()
            if (theta.alpha_gap != null) p.alphaGap = (theta.alpha_gap as Number).doubleValue()
            engine.register(p, (priorWeights.gap_filling ?: 1.0) as double)
        }
        if ('genomic_context' in enabled) {
            def p = new GenomicContextPrior()
            if (theta.alpha_ctx != null) p.alphaCtx = (theta.alpha_ctx as Number).doubleValue()
            if (theta.alpha_gap_ctx != null) p.alphaGapCtx = (theta.alpha_gap_ctx as Number).doubleValue()
            engine.register(p, (priorWeights.genomic_context ?: 1.0) as double)
        }
        if ('homology_transfer' in enabled) {
            def p = new gspa.integration.prior.HomologyTransferPrior()
            if (theta.alpha_homology != null) p.alpha = (theta.alpha_homology as Number).doubleValue()
            if (theta.homology_min_consensus != null) p.minConsensus = (theta.homology_min_consensus as Number).doubleValue()
            if (theta.homology_min_delta != null) p.minDelta = (theta.homology_min_delta as Number).doubleValue()
            if (theta.homology_max_boost != null) p.maxBoost = (theta.homology_max_boost as Number).doubleValue()
            engine.register(p, (priorWeights.homology_transfer ?: 1.0) as double)
        }
        engine
    }
}
