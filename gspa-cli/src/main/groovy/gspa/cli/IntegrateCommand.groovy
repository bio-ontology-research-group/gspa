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
import gspa.model.Genome
import gspa.model.GenomeLayout
import gspa.model.GenomeLayoutLoader
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

    @Option(names = ['--claims'], required = true,
            description = 'JSONL file of evidence claims (one per line).')
    File claimsFile

    @Option(names = ['--theta'],
            description = 'JSON file with reliability + prior-weight hyperparameters.')
    File thetaFile

    @Option(names = ['--out'], required = true,
            description = 'Output TSV path.')
    File outFile

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
            description = 'GO taxon constraints OBO file (go-computed-taxon-constraints.obo).')
    File taxonConstraintsFile

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
    String enabledPriors = 'essentiality,coherence,consistency,gap_filling,genomic_context'

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

    @Override
    void run() {
        println "GSPA integrate"
        println "  Claims: ${claimsFile}"
        println "  Theta:  ${thetaFile ?: '(defaults)'}"
        println "  Out:    ${outFile}"

        // Phase 10 flag validation.
        validatePhase10Flags()

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
        wireReferenceData(state)

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

        // --- Write output TSV ---
        outFile.parentFile?.mkdirs()
        outFile.withWriter { w ->
            w.writeLine(['protein_id', 'type', 'function_id', 'go_aspect',
                         'posterior_prob', 'likelihood_logodds', 'final_logodds',
                         'n_supporting', 'priors_fired', 'convergence_iter'].join('\t'))
            integrated.provenance.values().each { prov ->
                def fkey = IntegrationState.splitFunctionKey(prov.functionKey)
                if (fkey == null) return
                def fired = prov.priorContributions.collect { k, v ->
                    "${k}:${String.format(Locale.ROOT, '%.3f', v)}"
                }.join(',')
                def firstClaim = prov.supportingClaims.first()
                w.writeLine([
                    fkey[0],
                    fkey[1],
                    fkey[2],
                    firstClaim.goAspect ?: '',
                    String.format(Locale.ROOT, '%.6f', prov.finalProbability),
                    String.format(Locale.ROOT, '%.4f', prov.likelihoodLogOdds),
                    String.format(Locale.ROOT, '%.4f', prov.finalLogOdds),
                    prov.supportingClaims.size().toString(),
                    fired,
                    prov.convergenceIter.toString(),
                ].join('\t'))
            }
        }
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
                    n_supporting: prov.supportingClaims.size(),
                    supporting_sources: prov.supportingClaims*.source,
                    convergence_iter: prov.convergenceIter,
                ]
            }
            provenanceFile.parentFile?.mkdirs()
            mapper.writerWithDefaultPrettyPrinter().writeValue(provenanceFile, provList)
            println "  Provenance: ${provenanceFile}"
        }
    }

    /**
     * Wire optional reference data into the IntegrationState. Each piece is
     * independent: if a file isn't provided, the corresponding prior
     * no-ops.
     */
    private void wireReferenceData(IntegrationState state) {
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
        if (pathwaysFile != null && ec2goFile != null) {
            try {
                state.pathwayDatabase = gspa.ontology.PathwayLoader.loadPathways(
                    pathwaysFile,
                    gspa.ontology.PathwayLoader.loadEc2Go(ec2goFile),
                )
                println "  Pathways: ${state.pathwayDatabase.pathways.size()} loaded from ${pathwaysFile.name}"
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

        // SAT consistency checker: only wire when a taxonomy hierarchy is
        // provided (--taxonomy), so the SAT solver knows which taxon the
        // genome belongs to. Without the hierarchy, the checker can't
        // distinguish "this is E. coli" from "this could be any organism",
        // which causes massive false-positive violations.
        if (taxonomyFile != null) {
            try {
                def taxonConstraints = new gspa.ontology.TaxonConstraints()
                if (taxonConstraintsFile != null && taxonConstraintsFile.exists()) {
                    taxonConstraints.loadFromObo(taxonConstraintsFile)
                } else if (state.goOntology != null) {
                    taxonConstraints.loadFromGoOntology(state.goOntology)
                }
                if (taxonConstraints.constrainedTermCount() > 0) {
                    def checker = new gspa.ontology.SatConsistencyChecker(taxonConstraints)
                    checker.loadTaxonomyHierarchy(taxonomyFile)
                    state.satConsistencyChecker = checker
                    println "  SAT checker: ${taxonConstraints.onlyInTaxon.size()} only_in + ${taxonConstraints.neverInTaxon.size()} never_in (taxonomy from ${taxonomyFile})"
                }
            } catch (Exception e) {
                System.err.println "  [warn] Failed to wire SAT consistency checker: ${e.message}"
            }
        } else {
            println "  SAT checker: disabled (no --taxonomy provided; ConsistencyPrior requires genome taxon context)"
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
