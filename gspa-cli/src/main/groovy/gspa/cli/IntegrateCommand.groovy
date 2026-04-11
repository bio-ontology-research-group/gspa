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
import gspa.integration.PriorEngine
import gspa.integration.prior.ConsistencyPrior
import gspa.integration.prior.CoherencePrior
import gspa.integration.prior.EssentialityPrior
import gspa.integration.prior.GapFillingPrior
import gspa.integration.prior.GenomicContextPrior
import gspa.integration.suggester.DarkMatterSuggester
import gspa.integration.suggester.DisjunctiveSuggestion
import gspa.integration.suggester.SingletonSuggestion
import gspa.integration.suggester.Suggestion
import gspa.model.Genome
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
    description = 'Run the Phase 7 evidence integrator on pre-parsed claims.'
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

    @Option(names = ['--taxonomy'],
            description = 'NCBI taxonomy hierarchy file (for ConsistencyPrior).')
    File taxonomyFile

    @Option(names = ['--operons'],
            description = 'Operon assignments TSV (one operon per line, tab-separated protein IDs).')
    File operonsFile

    @Option(names = ['--gaps'],
            description = 'Metabolic gaps JSONL (one gap per line).')
    File gapsFile

    @Option(names = ['--enable-priors'],
            description = 'Which priors to enable (comma-separated). Default: all.')
    String enabledPriors = 'essentiality,coherence,consistency,gap_filling,genomic_context'

    @Option(names = ['--lite'], description = 'Skip ELK initialization (no process coherence).')
    boolean lite = false

    @Option(names = ['--dark-matter'],
            description = 'Enable the Phase 8 dark-matter / contextual-gap suggester.')
    boolean darkMatter = false

    @Option(names = ['--suggestions-out'],
            description = 'Suggestions TSV output path (only with --dark-matter).')
    File suggestionsOut

    @Override
    void run() {
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

        // --- Phase 8: optional dark-matter suggester ---
        if (darkMatter) {
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
            suggester.suggest(state, integrated)
            println "  Dark-matter suggester emitted ${integrated.suggestions.size()} suggestions"

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

        // Pathway database
        if (pathwaysFile != null && ec2goFile != null) {
            try {
                state.pathwayDatabase = gspa.ontology.PathwayLoader.loadPathways(
                    pathwaysFile,
                    gspa.ontology.PathwayLoader.loadEc2Go(ec2goFile),
                )
                println "  Pathways: ${state.pathwayDatabase.pathways.size()} loaded"
            } catch (Exception e) {
                System.err.println "  [warn] Failed to load pathway database: ${e.message}"
            }
        }

        // Taxonomy + SAT consistency checker
        if (taxonomyFile != null) {
            try {
                def checker = new gspa.ontology.SatConsistencyChecker()
                checker.loadTaxonomyHierarchy(taxonomyFile)
                state.satConsistencyChecker = checker
                println "  Consistency checker wired with taxonomy from ${taxonomyFile}"
            } catch (Exception e) {
                System.err.println "  [warn] Failed to load taxonomy: ${e.message}"
            }
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
        engine
    }
}
