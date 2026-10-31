package gspa.integration.workflow

import com.fasterxml.jackson.databind.ObjectMapper
import gspa.integration.CalibrationTable
import gspa.integration.EvidenceClaim
import gspa.integration.EvidenceCombiner
import gspa.integration.EvidenceType
import gspa.integration.IntegratedAnnotationSet
import gspa.integration.IntegrationRunAudit
import gspa.integration.IntegrationState
import gspa.integration.IntegrationWriter
import gspa.integration.IterativeRefiner
import gspa.integration.OuterIterativeRefiner
import gspa.integration.PriorEngine
import gspa.integration.contextbridge.CohortContextBuilder
import gspa.integration.contextbridge.CohortStage
import gspa.integration.contextbridge.CommunityStage
import gspa.integration.contextbridge.ExperimentStage
import gspa.integration.contextbridge.ExternalContextSupport
import gspa.integration.crossgenome.CrossGenomeReScorer
import gspa.integration.crossgenome.ReactionLocusCatalog
import gspa.integration.learned.ExternalContext
import gspa.integration.learned.FeatureExporter
import gspa.integration.learned.LearnedContextModel
import gspa.integration.learned.LearnedIntegrationStage
import gspa.integration.learned.LearnedStageReport
import gspa.integration.learned.LearnedStageSettings
import gspa.integration.prior.AnatomicalPrior
import gspa.integration.prior.CoherencePrior
import gspa.integration.prior.ConsistencyPrior
import gspa.integration.prior.EssentialityPrior
import gspa.integration.prior.GapFillingPrior
import gspa.integration.prior.GenomicContextPrior
import gspa.integration.promotion.AllAboveThresholdStrategy
import gspa.integration.promotion.BeamSearchStrategy
import gspa.integration.promotion.GreedyStrategy
import gspa.integration.promotion.MaxSatStrategy
import gspa.integration.promotion.PromotionStrategy
import gspa.integration.suggester.DarkMatterSuggester
import gspa.integration.suggester.DarkMatterWriter
import gspa.integration.suggester.ReactionLocalContextSuggester
import gspa.model.Genome
import gspa.model.GenomeLayout
import gspa.model.ProteinLocus

/**
 * The one authoritative execution of the inference workflow: evidence
 * integration, dark-matter suggestion, reaction-local context, cross-genome
 * rescoring, symbolic enforcement, metabolic refinement and agent completion,
 * always in that order, always recorded in an {@link IntegrationRunAudit}
 * (executed / skipped / failed / not_run).
 *
 * <p>Used by {@code gspa integrate} (single and batch) and by the annotation
 * pipeline. Callers differ only in where claims come from and what they do with
 * the result: the CLI exports TSV/JSON artifacts, the pipeline writes the
 * integrated posteriors back to its {@link Genome}. Core depends on neither.</p>
 *
 * <p>For a batch, load heavyweight references once with
 * {@link #loadSharedReferences} and pass them to every workflow; the per-genome
 * state, taxon checker and artifacts stay private to each run.</p>
 */
class IntegrationWorkflow {

    /** Everything one run produced. */
    static class Result {
        String genomeId
        IntegratedAnnotationSet integrated
        IntegrationState state
        IntegrationRunAudit audit
        List<EvidenceClaim> claims
        Map qualityReport
        int maxIter
        long elapsedNanos
        /** Typed context the learned stage and the feature export consumed (null: none was available). */
        ExternalContext context
        /** Cohort-context report of this genome (null: cohort inference was not requested). */
        Map<String, Object> cohortReport
        /** Learned-stage audit (null: the stage was not requested, scores are the legacy heuristic ones). */
        LearnedStageReport learned
        /** Summary of the feature export (null: not requested). */
        FeatureExporter.Summary featureExport
    }

    /** Everything the stages produce besides the annotation set. */
    private static class Outputs {
        ExternalContext context
        Map<String, Object> cohortReport
        LearnedStageReport learned
        FeatureExporter.Summary featureExport
    }

    final IntegrationWorkflowConfig c
    /** Heavyweight references loaded once for a batch (null: load them for this run). */
    IntegrationState sharedReferences
    /** Progress sink; the CLI prints, library callers log. */
    Closure<Void> info = { String m -> }
    /** Ledger of the most recent {@link #run}, available even when the run threw. */
    IntegrationRunAudit lastAudit
    /**
     * Result of the batch cohort pass (one joint inference over every genome). When set, this genome's
     * held-out features are taken from it instead of running a cohort of its own.
     */
    CohortContextBuilder.Result batchCohort
    /** Result of a batch-level community run (one community per batch); when set, this genome's ledger refers to it. */
    CommunityStage.Outcome batchCommunity

    private ExternalContext explicitContext
    private LearnedContextModel learnedModel
    /** Actual AntiFam outcome supplied by the annotation route, never inferred from absent hits. */
    String antiFamScanStatus

    IntegrationWorkflow(IntegrationWorkflowConfig config) {
        this.c = config
    }

    def propertyMissing(String name) { c.getProperty(name) }

    private void println(Object message) { info.call(String.valueOf(message)) }

    /**
     * Artifact loader and stage factory. The defaults are the production ones: they refuse missing,
     * mechanism-test and contract-mismatched artifacts. Only test sources replace them (to run a
     * hand-set mechanism model through the route); no main-source code does.
     */
    Closure<LearnedContextModel> learnedModelLoader = { File file -> LearnedContextModel.load(file) }
    Closure<LearnedIntegrationStage> learnedStageFactory = { LearnedStageSettings s -> new LearnedIntegrationStage(s) }

    static IntegrationState loadSharedReferences(IntegrationWorkflowConfig config, Closure<Void> info = { String m -> }) {
        config.validate()
        new IntegrationReferenceLoader(config).with { it.info = info; it.loadShared() }
    }

    /**
     * Run all requested stages on {@code claims}. With a {@code genome}, state is
     * keyed to that genome and the final posteriors replace its annotations;
     * without one the result is returned for export only.
     *
     * <p>The ledger exists from the first step: configuration validation and reference
     * loading are the {@code reference_preparation} stage, so a failure there leaves a
     * readable ledger (failed, with the reason; every requested stage not run) next to
     * the configured output and never a stale file from an earlier run.</p>
     */
    Result run(List<EvidenceClaim> claims, Genome genome = null) {
        long started = System.nanoTime()
        String genomeId = genome?.id ?: c.genomeId
        clearStaleArtifacts(c)
        def audit = new IntegrationRunAudit(genomeId)
        lastAudit = audit

        Map theta = null
        EvidenceCombiner combiner = null
        IntegrationState state = null
        IterativeRefiner refiner = null
        try {
            audit.runPreparation('reference_preparation') {
                c.validate()
                if (genome == null && c.darkMatterMode != DarkMatterWriter.MODE_SEPARATE) {
                    throw new IllegalArgumentException(
                        "darkMatterMode '${c.darkMatterMode}' modifies a genome's proteins; run it through annotate (a claims-only run keeps suggestions on the separate channel)")
                }
                theta = c.resolveTheta()
                // A requested model or context that cannot be loaded fails here, before any stage runs.
                explicitContext = c.externalContextFile != null ? ExternalContextSupport.read(c.externalContextFile) : null
                learnedModel = c.learnedRequested() ? learnedModelLoader.call(c.learnedModelFile) : null

                // --- Build combiner with learned reliability ---
                Map<EvidenceType, Double> reliability = EvidenceCombiner.defaultReliability()
                if (theta.reliability instanceof Map) {
                    (theta.reliability as Map).each { k, v ->
                        try {
                            reliability[EvidenceType.valueOf(k as String)] = (v as Number).doubleValue()
                        } catch (ignored) { /* unknown type */ }
                    }
                }
                combiner = new EvidenceCombiner(reliability)
                if (theta.second_opinion_bonus != null) {
                    combiner.secondOpinionBonus = (theta.second_opinion_bonus as Number).doubleValue()
                }

                // --- State with ontology, pathway, operon, gap and layout handles ---
                state = new IntegrationState(genome ?: new Genome(id: genomeId))
                def loader = new IntegrationReferenceLoader(c)
                loader.info = info
                loader.wire(state, claims, sharedReferences)
                if (state.genomeLayout == null && genome != null) state.genomeLayout = layoutFromGenome(genome)

                // --- Prior engine and refiner ---
                refiner = new IterativeRefiner(combiner)
                refiner.priorEngine = buildPriorEngine(theta, c.enabledPriorSet())
                if (theta.max_iter != null)  refiner.maxIter = (theta.max_iter as Number).intValue()
                if (theta.epsilon != null)   refiner.epsilon = (theta.epsilon as Number).doubleValue()
                if (theta.damping != null)   refiner.damping = (theta.damping as Number).doubleValue()
            }
        } catch (Exception ex) {
            audit.abandon(requestedStages(c), 'not run: reference_preparation failed')
            writeLedger(c, audit)
            throw ex
        }

        IntegratedAnnotationSet integrated
        Outputs outputs = new Outputs()
        try {
            integrated = runStages(audit, claims, state, refiner, theta, genomeId, outputs)
        } finally {
            writeLedger(c, audit)
        }

        if (genome != null) {
            IntegrationWriter.applyIntegratedAnnotations(genome, integrated)
            // Only dark-matter suggestions follow the dark-matter mode; other sources stay on the channel.
            if (c.darkMatter) {
                int touched = DarkMatterWriter.apply(genome, integrated, c.darkMatterMode, ['dark_matter'] as Set)
                println "  Dark-matter mode ${c.darkMatterMode}: ${touched} proteins updated"
            }
        }
        new Result(genomeId: genomeId, integrated: integrated, state: state, audit: audit, claims: claims,
            qualityReport: gspa.integration.GeneModelQuality.report(claims, antiFamScanStatus ?: outputs.context?.predictorStatus?.get('antifam')),
            maxIter: refiner.maxIter, elapsedNanos: System.nanoTime() - started, context: outputs.context,
            cohortReport: outputs.cohortReport, learned: outputs.learned, featureExport: outputs.featureExport)
    }

    /**
     * Ledger for a run that failed before it could start (batch-wide reference loading): every
     * requested stage not run, preparation failed with the reason; same file name as a normal run.
     */
    static IntegrationRunAudit recordPreparationFailure(IntegrationWorkflowConfig config, Throwable cause,
                                                        String stage = 'reference_preparation') {
        clearStaleArtifacts(config)
        def audit = new IntegrationRunAudit(config.genomeId)
        try {
            if (stage == 'reference_preparation') audit.runPreparation('reference_preparation') { throw cause }
            else {
                // A batch-wide stage (the cohort pass) failed after the references loaded.
                audit.runPreparation('reference_preparation') { }
                audit.runSupport(stage, true, null) { throw cause }
            }
        } catch (Throwable expected) { /* recorded in the ledger */ }
        audit.abandon(requestedStages(config), "not run: ${stage} failed".toString())
        writeLedger(config, audit)
        audit
    }

    /** Stage names this configuration asks for, in execution order. */
    static List<String> requestedStages(IntegrationWorkflowConfig c) {
        List<String> stages = []
        if (c.cohortRequested()) stages << 'cohort_context'
        stages << 'evidence_integration'
        if (c.darkMatter) stages << 'dark_matter'
        if (c.rlcSuggester) stages << 'reaction_context'
        if (c.rxnLocusCatalogFile != null) stages << 'cross_genome'
        if (c.enforceStage2) stages << 'symbolic_enforcement'
        if (c.enableFbaRefiner) stages << 'metabolic_refinement'
        if (c.enableAgent) stages << 'agent_completion'
        if (c.experimentsRequested()) stages << 'metabolic_experiments'
        if (c.communityRequested()) stages << 'community_metabolism'
        if (c.learnedRequested()) stages << 'learned_integration'
        if (c.featureExportDir != null) stages << 'feature_export'
        stages
    }

    private static void writeLedger(IntegrationWorkflowConfig c, IntegrationRunAudit audit) {
        if (c.outFile != null) audit.write(new File(c.outFile.absolutePath + '.stages.json'))
    }

    /** Outputs of an earlier run at the configured paths must not survive a run that fails. */
    private static void clearStaleArtifacts(IntegrationWorkflowConfig c) {
        List<File> files = [c.provenanceFile, c.suggestionsOut, c.featuresOut, c.enforcementActionsOut,
                            c.anomaliesOut, c.agentResolutionsOut]
        if (c.outFile != null) {
            files << c.outFile
            ['.stages.json', '.rescues.tsv', '.modules.json', '.learned.json', '.cohort.json', '.quality.json'].each { files << new File(c.outFile.absolutePath + it) }
            // Only files a stage writes are removed (and the directory if that leaves it empty): a directory
            // that happens to carry the name but holds anything else is left alone.
            ['.experiments', '.community'].each { String suffix ->
                File d = new File(c.outFile.absolutePath + suffix)
                ['manifest.json', 'result.json', 'sidecar.log'].each { new File(d, it).delete() }
                if (suffix == '.experiments') new File(d, 'construction.json').delete()
                d.delete()
            }
        }
        if (c.featureExportDir != null) {
            [FeatureExporter.ROWS_FILE, FeatureExporter.CONTRACT_FILE, FeatureExporter.MANIFEST_FILE].each { files << new File(c.featureExportDir, it) }
        }
        files.findAll { it != null }.each { it.delete() }
    }


    /** Writes the configured per-genome artifacts of a finished run (no-op without paths). */
    void writeArtifacts(Result result) {
        if (c.outFile == null) return
        IntegrationArtifacts.writeIntegratedOutputTsv(result.integrated, c.outFile)
        println "  Wrote: ${c.outFile}"
        if (result.qualityReport != null) {
            new ObjectMapper().writerWithDefaultPrettyPrinter().writeValue(new File(c.outFile.absolutePath + '.quality.json'), result.qualityReport)
        }
        [IntegrationArtifacts.writeRescueLog(result.state, c.outFile),
         IntegrationArtifacts.writeModuleReport(result.state, result.integrated, c.outFile)].each {
            if (it != null) println it
        }
        if (c.provenanceFile != null) {
            IntegrationArtifacts.writeProvenanceJson(result.integrated, c.provenanceFile)
            println "  Provenance: ${c.provenanceFile}"
        }
        if (result.cohortReport != null) {
            File f = new File(c.outFile.absolutePath + '.cohort.json')
            new ObjectMapper().writerWithDefaultPrettyPrinter().writeValue(f, result.cohortReport)
            println "  Cohort context report: ${f}"
        }
        if (result.learned != null) {
            File f = IntegrationArtifacts.writeLearnedReport(result.learned, new File(c.outFile.absolutePath + '.learned.json'))
            println "  Learned stage report: ${f} (${result.learned})"
        }
    }

    private IntegratedAnnotationSet runStages(IntegrationRunAudit audit, List<EvidenceClaim> claims,
                                              IntegrationState state, IterativeRefiner refiner, Map theta,
                                              String genomeId, Outputs outputs) {
        IntegratedAnnotationSet integrated

        // --- Cohort context from the pre-cohort Stage 0 snapshot (raw claims only, independent of all later stages) ---
        ExternalContext cohortContext = null
        audit.runSupport('cohort_context', c.cohortRequested(),
            state.reactionGraph == null ? 'reaction graph is required (--reaction-graph)' : null) {
            CohortContextBuilder.Result cohort = batchCohort
            String route = 'batch'
            if (cohort == null) {
                route = 'single_run'
                def stage = new CohortStage()
                stage.info = info
                cohort = stage.run(c, [new CohortStage.Target(genomeId: genomeId, claims: { claims }, layout: state.genomeLayout,
                    orthogroupsFile: c.orthogroupsFile)], state.reactionGraph)
            }
            cohortContext = cohort.contexts.get(genomeId)
            outputs.cohortReport = cohort.reports.get(genomeId)
            if (outputs.cohortReport == null) {
                throw new IllegalStateException("the cohort result has no report for genome '${genomeId}' (is it a target of the batch?)")
            }
            Map report = outputs.cohortReport
            [route: route, baseline: 'stage0_raw_max_score', baseline_snapshot_id: cohort.snapshotId,
             cohort_genomes: (report.cohort_genomes as List).size(), n_peer_genomes: report.n_peer_genomes,
             keys_with_cohort_features: (report.candidates as Map).keys_with_cohort_features,
             keys_with_local_features: (report.candidates as Map).keys_with_local_features,
             no_peers: report.no_peers,
             consumed_by: [c.learnedRequested() ? 'learned_integration' : null, c.featureExportDir != null ? 'feature_export' : null].findAll { it }]
        }
        outputs.context = ExternalContextSupport.merge(explicitContext, cohortContext)

        // --- Refine ---
        integrated = audit.runStage('evidence_integration', true, null, null, state) {
            refiner.refine(claims, state)
        }
        println "  Produced ${integrated.annotations.size()} integrated annotations"

        // --- Phase 8: optional dark-matter suggester ---
        integrated = audit.runStage('dark_matter', c.darkMatter, null, integrated, state) {
            def suggester = buildSuggester(theta)
            suggester.useRefinedBayesFactor = c.refinedBf

            if (c.iterateGapseq) {
                // Phase 10 outer fixed-point loop over (refine → suggest → promote → pin).
                println "  Running Phase 10 outer loop (maxIter=${c.maxGapseqIter}, qBase=${c.gapseqQBase}, qStep=${c.gapseqQStep}, pin=${c.gapseqPinPromotions}, strategy=${c.promotionStrategy}, refinedBf=${c.refinedBf})"
                def outer = new OuterIterativeRefiner(refiner)
                outer.suggester = suggester
                outer.maxIter = c.maxGapseqIter
                outer.qBase = c.gapseqQBase
                outer.qStep = c.gapseqQStep
                outer.qCap = c.gapseqQCap
                outer.pinPromotions = c.gapseqPinPromotions
                outer.promotionStrategy = buildPromotionStrategy()
                def gs = new OuterIterativeRefiner.CoverageGapSource(tauCover: c.gapseqTauCover, pathwayDb: state.pathwayDatabase)
                outer.gapSource = gs
                def outerResult = outer.refine(claims, state)
                integrated = outerResult.integrated
                println "  Outer loop: iter=${outerResult.outerIterationsRun}, fixedPoint=${outerResult.fixedPointReached}, cascade=${outerResult.cascadeRolledBack}, promoted_per_iter=${outerResult.promotedPerIter}, gaps_per_iter=${outerResult.gapsPerIter}"
            } else {
                suggester.suggest(state, integrated)
                println "  Dark-matter suggester emitted ${integrated.suggestions.size()} suggestions"
            }

            return integrated
        }

        // --- Phase 12: optional Reaction-Local Context suggester (M1) ---
        integrated = audit.runStage('reaction_context', c.rlcSuggester,
            state.reactionGraph == null ? 'reaction graph is required' :
                (state.genomeLayout == null ? 'genome layout is required' : null), integrated, state, {
            def rlc = new ReactionLocalContextSuggester()
            rlc.radiusK = c.rlcRadiusK
            rlc.alpha = c.rlcAlpha
            rlc.kernelBandwidth = c.rlcKernelBandwidth
            rlc.anchorPosteriorThreshold = c.rlcAnchorThreshold
            if (c.featuresOut != null) rlc.featuresOut = c.featuresOut
            rlc.suggest(state, integrated)
            println "  RLGC suggester emitted ${integrated.suggestions.size()} suggestions"

            return integrated
        }, state.metabolicGaps ? null : 'no metabolic gaps to evaluate')
        integrated = audit.runStage('cross_genome', c.rxnLocusCatalogFile != null,
            c.rxnLocusCatalogFile != null && !c.rxnLocusCatalogFile.isFile() ? 'catalog file does not exist' :
                (!c.rlcSuggester && !c.darkMatter ? 'enable a context suggester to rescore candidates' : null), integrated, state) {
            def catalog = ReactionLocusCatalog.readFrom(c.rxnLocusCatalogFile)
            def cgr = new CrossGenomeReScorer(lambda: c.cgLambda, minSupport: c.cgMinSupport,
                requireCredible: c.cgRequireCredible)
            cgr.rescore(state, integrated, catalog)
            return integrated
        }

        if (c.suggestionsOut != null && (c.rlcSuggester || c.darkMatter)) {
            c.suggestionsOut.parentFile?.mkdirs()
            IntegrationArtifacts.writeSuggestionsTsv(integrated.suggestions, c.suggestionsOut)
        }

        // --- Optional Stage 2A: Symbolic Constraint & Coherence Enforcement ---
        integrated = audit.runStage('symbolic_enforcement', c.enforceStage2,
            state.satConsistencyChecker == null ? 'taxon consistency checker is unavailable' : null, integrated, state) {
            def stage = new SymbolicEnforcementStage(c)
            stage.info = info
            stage.run(state, integrated)
            return integrated
        }

        // --- Optional Stage 2B: StoichiometricLocusRefiner ---
        integrated = audit.runStage('metabolic_refinement', c.enableFbaRefiner,
            c.fbaDataDir == null || !c.fbaDataDir.isDirectory() ? 'FBA reference data directory is required' :
                (c.outFile == null ? 'an output path is required for model artifacts' :
                    (!hasEcCandidates(integrated) ? 'no EC predictions to reconstruct a metabolic model from' : null)),
            integrated, state) {
            def fbaRefiner = new gspa.integration.fba.StoichiometricLocusRefiner(
                fbaDataDir: c.fbaDataDir,
                sidecarScript: c.fbaSidecar,
                pythonBin: c.fbaPythonBin,
                gram: c.fbaGram,
                medium: c.fbaMedium
            )
            def fbaSum = fbaRefiner.refine(state, integrated, new File(c.outFile.absolutePath + ".metabolism"))
            // The medium and biomass template are assumptions of the model, not observations.
            audit.note('metabolic_refinement', [medium: c.fbaMedium, gram: c.fbaGram,
                core_reactions: fbaSum.nCoreRxns, repaired_reactions: fbaSum.nRepairedRxns,
                uncertified_additions: fbaSum.nUncertifiedAdditions, updated_loci: fbaSum.nUpdatedLoci])
            println "  Stage 2B StoichiometricLocusRefiner (gram=${c.fbaGram}, medium=${c.fbaMedium}): coreRxns=${fbaSum.nCoreRxns}, repairedRxns=${fbaSum.nRepairedRxns}, uncertified=${fbaSum.nUncertifiedAdditions}, updatedLoci=${fbaSum.nUpdatedLoci} (${String.format(Locale.ROOT, '%.2f', fbaSum.solveTimeSec)}s)"
            return integrated
        }

        // --- Optional Stage 3: Neuro-Symbolic Agentic Completion (gspa-agent) ---
        if (c.anomaliesOut != null && !c.enableAgent) {
            def extractorAnom = new gspa.agent.AnomalyExtractor()
            def anomalies = extractorAnom.extract(state, integrated)
            gspa.agent.AnomalyExtractor.writeJsonl(anomalies, c.anomaliesOut)
            println "  Extracted ${anomalies.size()} anomalies to ${c.anomaliesOut}"
        }
        integrated = audit.runStage('agent_completion', c.enableAgent,
            c.outFile == null && (c.anomaliesOut == null || c.agentResolutionsOut == null) ?
                'an output path (or explicit anomaly and resolution paths) is required' : null, integrated, state) {
            File anomFile = c.anomaliesOut ?: new File(c.outFile.parentFile ?: new File('.'), "${c.outFile.name}.anomalies.jsonl")
            File resFile = c.agentResolutionsOut ?: new File(c.outFile.parentFile ?: new File('.'), "${c.outFile.name}.agent_resolutions.jsonl")
            def runner = new gspa.agent.AgentCompletionRunner(
                sidecarScript: c.agentSidecar,
                agentModel: c.agentModel,
                pythonBin: c.agentPythonBin,
                apiBaseUrl: c.agentApiBase,
                guardrail: !c.agentNoGuardrail
            )
            def sum = runner.run(state, integrated, anomFile, resFile)
            println "  Stage 3 Agent (${c.agentModel}, guardrail=${!c.agentNoGuardrail}): anomalies=${sum.anomaliesExtracted}, proposed=${sum.resolutionsProposed}, promoted=${sum.promoted}, relaxed=${sum.relaxed}, demoted=${sum.demoted}, abstained=${sum.abstained}, rejected_by_guardrail=${sum.rejectedByGuardrail} (${String.format(Locale.ROOT, '%.2f', sum.elapsedSeconds)}s)"
            return integrated
        }

        // --- Metabolic experiments: conditional model results for this run's candidates (never a probability) ---
        audit.runSupport('metabolic_experiments', c.experimentsRequested(),
            c.outFile == null ? 'an output path is required for experiment artifacts' : null) {
            File refinerModel = audit.stages['metabolic_refinement'].status == 'executed' ?
                new File(c.outFile.absolutePath + '.metabolism/metabolic_model.xml') : null
            def stage = new ExperimentStage()
            stage.info = info
            ExperimentStage.Outcome outcome = stage.run(c, new File(c.outFile.absolutePath + '.experiments'),
                refinerModel != null && refinerModel.isFile() ? refinerModel : null, integrated, claims, outputs.context,
                { String ec -> (state.reactionGraph?.reactionsForEc(ec) ?: []) as Collection<String> })
            if (!outcome.fba.isEmpty()) {
                outputs.context = ExternalContextSupport.merge(outputs.context, new ExternalContext(fba: outcome.fba), 'metabolic experiments')
            }
            outcome.facts
        }

        // --- Community cross-feeding: organism-specific, no gene assignment, no candidate feature ---
        audit.runSupport('community_metabolism', c.communityRequested(),
            c.outFile == null && batchCommunity == null ? 'an output path is required for community artifacts' : null) {
            if (batchCommunity != null) {
                return batchCommunity.facts + [route: 'batch', note: 'one community run per batch; its result is shared by every genome ledger']
            }
            def stage = new CommunityStage()
            stage.info = info
            stage.run(c, new File(c.outFile.absolutePath + '.community')).facts + [route: 'single_run']
        }

        // --- Learned stage: LAST, over the ORIGINAL raw claims. Its calibrated probabilities replace the
        //     heuristic scores of every key it decides (never added to them). ---
        integrated = audit.runStage('learned_integration', c.learnedRequested(), null, integrated, state) {
            Set<String> forbidden = state.satConsistencyChecker?.organismTaxon ?
                (state.satConsistencyChecker.forbiddenConstrainedTerms ?: ([] as Set)) : ([] as Set)
            def settings = new LearnedStageSettings(model: learnedModel, retractBelow: c.learnedRetractBelow,
                addUnretainedClaims: c.learnedAddUnretainedClaims, forbiddenGoTerms: new LinkedHashSet<String>(forbidden))
            LearnedStageReport report = learnedStageFactory.call(settings).run(claims, state, integrated, outputs.context)
            outputs.learned = report
            audit.scoringMode = report.candidates && report.unsupportedCount() == 0 ? 'learned' :
                (report.candidates.size() > report.unsupportedCount() ? 'mixed' : 'legacy_heuristic')
            audit.note('learned_integration', [
                model_id: report.modelId, model_sha256: report.modelSha256, model_kind: report.modelKind,
                contract_id: report.contractId, contract_fingerprint: report.contractFingerprint,
                candidates_scored: report.candidates.size() - report.unsupportedCount(),
                candidates_unsupported: report.unsupportedCount(), unchanged_unsupported: report.count('unchanged_unsupported'),
                rescored: report.count('rescored'), added: report.count('added'),
                retracted: report.count('retracted'), removed_no_projection: report.count('removed_no_projection'),
                scored_only: report.count('scored_only'), pinned_floors_cleared: report.pinnedFloorsCleared.size(),
                retract_below: c.learnedRetractBelow, add_unretained_claims: c.learnedAddUnretainedClaims,
                structurally_forbidden_go_terms: forbidden.size(),
                context_supplied: outputs.context != null,
                replaced_baseline: 'legacy_prior_baseline',
                note: 'Supported candidates were scored by the frozen model. Projected learned scores replace the heuristic posterior; scored-only candidates remain in the report. A display threshold, if any, is advisory.'])
            println "  Learned stage (${report.modelId}, sha256 ${report.modelSha256.take(12)}): ${report}"
            return integrated
        }
        if (outputs.learned != null && outputs.learned.unsupportedCount() == outputs.learned.candidates.size()) {
            audit.note('learned_integration', [status: 'skipped', requested: true,
                reason: 'no candidates supported by the model input-coverage manifest; baseline scores retained'])
        }
        if (!c.learnedRequested()) {
            audit.note('learned_integration', [reason: 'learned stage not requested: scores are the legacy heuristic ones (evidence combiner and hand-set priors); no trained model produced them'])
        }

        // --- Feature export: production rows from the same extractor, with or without a model; never any labels ---
        audit.runSupport('feature_export', c.featureExportDir != null, null) {
            def exporter = new FeatureExporter(learnedModel)
            FeatureExporter.Summary s = exporter.exportToDirectory(claims, state, outputs.context, c.featureExportDir, integrated, true, genomeId)
            outputs.featureExport = s
            println "  Feature export: ${s.rows} rows over ${s.functions} functions -> ${c.featureExportDir}"
            [rows: s.rows, functions: s.functions, rows_with_sequence_evidence: s.rowsWithSequenceEvidence,
             rows_sha256: s.rowsSha256, contract_id: s.contractId, contract_fingerprint: s.contractFingerprint,
             model_sha256: s.modelSha256, labels_included: false, directory: c.featureExportDir.path,
             note: 'Rows carry no label, split or group; join independently curated truth and declared splits by candidate_key.']
        }
        integrated
    }

    private static boolean hasEcCandidates(IntegratedAnnotationSet integrated) {
        integrated.provenance.keySet().any { String key ->
            String[] parts = IntegrationState.splitFunctionKey(key)
            parts != null && parts[1] == 'EC'
        }
    }

    /** Protein loci from an in-memory genome, for genomes that carry CDS coordinates. */
    static GenomeLayout layoutFromGenome(Genome genome) {
        GenomeLayout layout = new GenomeLayout(genome.id)
        genome.proteins.each { p ->
            if (p.sourceFeature == null || p.contig == null || p.end < p.start || p.end == 0) return
            layout.add(new ProteinLocus(genomeId: genome.id, proteinId: p.id, contig: p.contig.id,
                start: p.start, end: p.end, strand: p.strand))
        }
        if (layout.size() == 0) return null
        layout.finishLoading()
        layout
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
            if (theta.ess_rescue_mode != null) p.rescueMode = theta.ess_rescue_mode as String
            if (theta.ess_rescue_target != null) p.rescueTarget = (theta.ess_rescue_target as Number).doubleValue()
            if (theta.ess_min_candidate != null) p.minCandidate = (theta.ess_min_candidate as Number).doubleValue()
            if (theta.ess_max_rescue != null) p.maxRescuePerCategory = (theta.ess_max_rescue as Number).intValue()
            if (theta.ess_w_support != null) p.wSupport = (theta.ess_w_support as Number).doubleValue()
            if (theta.ess_w_free != null) p.wFree = (theta.ess_w_free as Number).doubleValue()
            if (theta.ess_w_competing != null) p.wCompeting = (theta.ess_w_competing as Number).doubleValue()
            if (theta.ess_w_context != null) p.wContext = (theta.ess_w_context as Number).doubleValue()
            engine.register(p, (priorWeights.essentiality ?: 1.0) as double)
        }
        if ('coherence' in enabled) {
            def p = new CoherencePrior()
            if (theta.alpha_coh != null) p.alphaCoh = (theta.alpha_coh as Number).doubleValue()
            if (theta.alpha_complex_boost != null) p.complexPartnerWeight = (theta.alpha_complex_boost as Number).doubleValue()
            if (theta.alpha_complex_penalty != null) p.complexSingletonPenaltyWeight = (theta.alpha_complex_penalty as Number).doubleValue()
            if (theta.coh_pathway_mode != null) p.pathwayMode = theta.coh_pathway_mode as String
            if (theta.coh_rescue_target != null) p.rescueTarget = (theta.coh_rescue_target as Number).doubleValue()
            if (theta.coh_refine_min_completeness != null) p.refineMinCompleteness = (theta.coh_refine_min_completeness as Number).doubleValue()
            if (theta.coh_max_rescue_per_module != null) p.maxRescuePerModule = (theta.coh_max_rescue_per_module as Number).intValue()
            if (theta.coh_min_candidate != null) p.minCandidate = (theta.coh_min_candidate as Number).doubleValue()
            if (theta.coh_min_score != null) p.minScore = (theta.coh_min_score as Number).doubleValue()
            if (theta.coh_w_support != null) p.wSupport = (theta.coh_w_support as Number).doubleValue()
            if (theta.coh_w_free != null) p.wFree = (theta.coh_w_free as Number).doubleValue()
            if (theta.coh_w_competing != null) p.wCompeting = (theta.coh_w_competing as Number).doubleValue()
            if (theta.coh_w_operon != null) p.wOperon = (theta.coh_w_operon as Number).doubleValue()
            if (theta.coh_w_module != null) p.wModule = (theta.coh_w_module as Number).doubleValue()
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
