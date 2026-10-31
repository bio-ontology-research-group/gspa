package gspa.integration.contextbridge

import com.fasterxml.jackson.databind.ObjectMapper
import gspa.integration.EvidenceClaim
import gspa.integration.IntegratedAnnotationSet
import gspa.integration.MetabolicGap
import gspa.integration.learned.ExternalContext
import gspa.integration.workflow.IntegrationWorkflowConfig

/** Run individually verified absence comparisons with one shared budget and deadline. */
class GuidedAbsenceStage {
    static ExperimentStage.Outcome run(IntegrationWorkflowConfig config, File workDir, Map template,
                                      SbmlModelInfo model, IntegratedAnnotationSet integrated,
                                      Collection<EvidenceClaim> claims, ExternalContext context,
                                      Collection<MetabolicGap> gaps, boolean gapsAvailable,
                                      Closure<Collection<String>> graphLookup, Closure info) {
        if (config.experimentEvidenceFile != null) {
            throw new IllegalArgumentException('guided absence uses this run\'s gaps/claims/context and explicit per-reaction supply; experimentEvidenceFile is not supported')
        }
        long start = System.nanoTime()
        ObjectMapper mapper = new ObjectMapper()
        Map built = GuidedAbsenceBuilder.build(template, model, claims, context, gaps, gapsAvailable,
            config.experimentMaxHypotheses, graphLookup)
        workDir.mkdirs()
        mapper.writerWithDefaultPrettyPrinter().writeValue(new File(workDir, 'construction.json'), built.report)
        Map result = [schema: 'gspa.guided_absence_comparisons/1', input_model_sha256: model.sha256,
            max_solver_calls: config.experimentMaxSolverCalls, max_experiments: config.experimentMaxExperiments,
            comparisons: [], unvisited: [], solver_calls_used: 0, conditions_visited: 0,
            promotion: 'none', calibrated_probability: null, selected_explanation: null, gene_attribution: 'unverified']
        File resultFile = new File(workDir, 'result.json')
        Closure save = { mapper.writerWithDefaultPrettyPrinter().writeValue(resultFile, result) }
        int calls = 0, conditions = 0
        built.comparisons.eachWithIndex { Map item, int index ->
            int seconds = config.experimentTimeoutSeconds - (int) Math.ceil((System.nanoTime() - start) / 1e9d)
            String unavailable = calls >= config.experimentMaxSolverCalls ? 'solver_budget_exhausted' :
                conditions >= config.experimentMaxExperiments ? 'experiment_budget_exhausted' :
                seconds < 1 ? 'deadline_exhausted' : null
            if (unavailable != null) {
                result.unvisited << [id: item.id, reason: unavailable, conditions: item.request.conditions*.id]
                return
            }
            File childDir = new File(workDir, String.format(Locale.ROOT, 'comparison-%04d', index))
            childDir.mkdirs()
            Map manifest = new LinkedHashMap(template)
            manifest.model = model.model.absolutePath
            manifest.absence_comparison = item.request
            File manifestFile = new File(childDir, 'request.json')
            mapper.writerWithDefaultPrettyPrinter().writeValue(manifestFile, manifest)
            IntegrationWorkflowConfig child = config.copy()
            child.experimentMode = 'manifest'
            child.experimentManifestFile = manifestFile
            child.experimentMaxSolverCalls = config.experimentMaxSolverCalls - calls
            child.experimentMaxExperiments = config.experimentMaxExperiments - conditions
            child.experimentTimeoutSeconds = seconds
            try {
                ExperimentStage.Outcome checked = new ExperimentStage(info: info).run(child, childDir, null,
                    integrated, claims, context, graphLookup)
                calls += checked.facts.solver_calls_used as int
                conditions += checked.facts.conditions_visited.size()
                result.comparisons << [id: item.id, status: 'verified', result_file: checked.resultFile.path,
                    result_sha256: SidecarRun.sha256(checked.resultFile), facts: checked.facts]
                result.solver_calls_used = calls
                result.conditions_visited = conditions
                save.call()
            } catch (Exception failure) {
                result.comparisons << [id: item.id, status: 'failed', reason: "${failure.class.simpleName}: ${failure.message}".toString()]
                result.outcome = 'failed'
                result.solver_budget_status = 'last comparison unverified; no subsequent comparison started'
                result.unvisited.addAll(built.comparisons.drop(index + 1).collect { [id: it.id, reason: 'previous_comparison_failed'] })
                save.call()
                throw failure
            }
        }
        result.outcome = !built.comparisons ? 'no_admissible_experiment' :
            !result.unvisited && result.comparisons.every { it.facts.outcome == 'completed' } && !built.report.truncated ?
                'completed' : 'completed_with_skips_or_inconclusive'
        save.call()
        new ExperimentStage.Outcome(solverRun: calls > 0, resultFile: resultFile, facts: [
            status: calls > 0 ? 'executed' : 'skipped', reason: calls > 0 ? null : 'no complete comparison condition ran; inspect unavailable or budget-limited comparisons',
            analysis_type: 'guided_absence_comparisons', context_source: built.report.context_source,
            context: built.report, mode: config.experimentMode, strategy: config.experimentStrategy,
            model_file: model.model.path, model_sha256: model.sha256, model_origin: 'supplied model; not this run\'s refiner export',
            solver_calls_used: calls, solver_calls_max: config.experimentMaxSolverCalls,
            conditions_visited: conditions, max_experiments: config.experimentMaxExperiments,
            comparisons: result.comparisons, unvisited_comparisons: result.unvisited,
            outcome: result.outcome, candidates_with_conditional_features: 0, promotion: 'none',
            gene_attribution: 'unverified', selected_explanation: null,
            result_file: resultFile.path, result_sha256: SidecarRun.sha256(resultFile),
            interpretation: built.report.interpretation])
    }
}
