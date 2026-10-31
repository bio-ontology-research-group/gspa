package gspa.integration.contextbridge

import com.fasterxml.jackson.databind.ObjectMapper
import gspa.integration.EvidenceClaim
import gspa.integration.IntegratedAnnotationSet
import gspa.integration.learned.ExternalContext
import gspa.integration.workflow.IntegrationWorkflowConfig

/**
 * The production entry point of {@code scripts/metabolic_experiments.py}: builds (or takes) a frozen
 * manifest, runs the sidecar under a hard solver-call budget and wall-clock limit, validates what came back
 * and maps the conditional records onto {@link ExternalContext.FbaConditional}.
 *
 * <p><b>Two honestly different modes.</b> {@code manifest} runs a user-supplied manifest exactly as written;
 * the ledger records {@code context_source: user_manifest} and claims no automated evidence guidance.
 * {@code evidence_guided} treats the file as a template and adds hypotheses and evidence from this run's
 * candidates ({@link ExperimentManifestBuilder}); with nothing admissible it records
 * {@code no_admissible_experiment} and starts no solver.</p>
 *
 * <p>Nothing returned is a probability or a promotion. Results of the sidecar are checked to say so
 * ({@code calibrated_probability: null}, {@code promotion: none}, never independent proof) and to belong to the
 * model that was supplied (SHA-256) and the budget that was requested.</p>
 *
 * <p>The optional Qwen strategy is only forwarded as sidecar options ({@code --strategy qwen --qwen-model
 * --api-base}); the sidecar reads its own credentials. The agent can only choose a menu id.</p>
 */
class ExperimentStage {

    static final String SCHEMA = 'gspa.metabolic_experiments/1'
    static final String MASKED_SCHEMA = 'gspa.masked_assignment_experiments/1'
    private static final Set<String> KNOWN_OUTCOMES = ['completed', 'completed_with_skips_or_inconclusive', 'no_candidate_experiments'] as Set

    static class Outcome {
        /** Facts for the stage ledger. */
        Map<String, Object> facts = [:]
        Map<String, ExternalContext.FbaConditional> fba = [:]
        File resultFile
        boolean solverRun
    }

    Closure<Void> info = { String m -> }

    /**
     * @param workDir receives manifest.json, result.json (the sidecar's output, verbatim) and sidecar.log
     * @param refinerModel the model exported by the metabolic_refinement stage of this run, or null
     */
    Outcome run(IntegrationWorkflowConfig c, File workDir, File refinerModel, IntegratedAnnotationSet integrated,
                Collection<EvidenceClaim> claims, ExternalContext context, Closure<Collection<String>> graphEcToReactions = null) {
        File templateFile = c.experimentManifestFile
        Map template
        try {
            template = new ObjectMapper().readValue(templateFile, Map)
        } catch (Exception e) {
            throw new IllegalArgumentException("experiment manifest ${templateFile} is not valid JSON: ${e.message}", e)
        }
        boolean guided = c.experimentMode == 'evidence_guided'
        boolean maskedAssignments = template.containsKey('assignment_recovery')
        File model = resolveModel(template, templateFile, refinerModel, guided)
        if (maskedAssignments && guided && model == refinerModel) {
            throw new IllegalArgumentException('guided masked assignments require a supplied model whose topology was not selected from the target annotations; a same-run refiner export is not independent')
        }
        File refinerTable = (refinerModel != null && model == refinerModel) ? new File(refinerModel.absoluteFile.parentFile, 'reactions.tsv') : null
        SbmlModelInfo modelInfo = SbmlModelInfo.read(model, refinerTable?.isFile() ? refinerTable : null)

        workDir.mkdirs()
        Outcome out = new Outcome()
        Map<String, Object> report
        List<ExperimentManifestBuilder.Binding> bindings
        Map manifest
        if (guided && maskedAssignments) {
            List<Map> evidence = c.experimentEvidenceFile != null ? ExperimentManifestBuilder.readEvidenceFile(c.experimentEvidenceFile) : []
            def built = new MaskedAssignmentBuilder().build(template, modelInfo, claims, context, evidence,
                c.experimentMaxHypotheses, graphEcToReactions)
            manifest = built.manifest
            report = built.report
            bindings = []
            new ObjectMapper().writerWithDefaultPrettyPrinter().writeValue(new File(workDir, 'construction.json'), report)
            if (!built.admissible()) {
                out.facts = baseFacts(c, modelInfo, report, c.experimentMode) + [
                    outcome: 'no_admissible_experiment', solver_calls_used: 0, analysis_type: 'masked_assignment',
                    promotion: 'none', reason: 'no measured candidate maps to an admissible assignment; null alone is not a comparison']
                return out
            }
        } else if (guided) {
            List<Map> evidence = c.experimentEvidenceFile != null ? ExperimentManifestBuilder.readEvidenceFile(c.experimentEvidenceFile) : []
            def built = new ExperimentManifestBuilder().build(template, modelInfo, integrated, claims, context, evidence, c.experimentMaxHypotheses, graphEcToReactions)
            manifest = built.manifest
            report = built.report
            bindings = built.bindings
            if (!built.admissible() && ((manifest.context.evidence ?: []) as List).isEmpty()) {
                out.facts = baseFacts(c, modelInfo, report, 'evidence_guided') + [
                    outcome: 'no_admissible_experiment', solver_calls_used: 0,
                    reason: 'no candidate of this run maps to a gene and reaction of the model, and no evidence record was supplied; nothing was invented and no solver was started']
                return out
            }
        } else {
            manifest = new LinkedHashMap(template)
            manifest.model = model.absolutePath
            bindings = bindingsFromManifest(manifest)
            report = [context_source: 'user_manifest',
                      note: 'the context was written by the user; this is not an evidence-guided workflow',
                      hypotheses: ((manifest.context ?: [:]).hypotheses ?: []).size(),
                      hypotheses_with_function_key: bindings.size(),
                      model_path_rewritten_to_absolute: true]
        }

        File manifestFile = new File(workDir, 'manifest.json')
        new ObjectMapper().writerWithDefaultPrettyPrinter().writeValue(manifestFile, manifest)
        File result = new File(workDir, 'result.json')
        List<String> cmd = [c.experimentPythonBin, c.experimentSidecar.absolutePath, '--manifest', manifestFile.absolutePath,
            '--out', result.absolutePath, '--strategy', c.experimentStrategy, '--seed', c.experimentSeed.toString(),
            '--max-experiments', c.experimentMaxExperiments.toString(), '--max-solver-calls', c.experimentMaxSolverCalls.toString()]
        if (c.experimentStrategy == 'qwen') {
            cmd.addAll(['--qwen-model', c.agentModel])
            if (c.agentApiBase) cmd.addAll(['--api-base', c.agentApiBase])
        }
        if (!c.experimentSidecar.isFile()) throw new IllegalArgumentException("experiment sidecar does not exist: ${c.experimentSidecar}")
        info.call("  Metabolic experiments (${c.experimentMode}, strategy ${c.experimentStrategy}, <= ${c.experimentMaxSolverCalls} solver calls)".toString())
        SidecarRun.Result run = SidecarRun.run(cmd, result, new File(workDir, 'sidecar.log'), c.experimentTimeoutSeconds,
            maskedAssignments ? MASKED_SCHEMA : SCHEMA,
            { Map parsed -> if (maskedAssignments) verifyMasked(parsed, modelInfo, c, manifest) else verify(parsed, modelInfo, c) })
        out.solverRun = true
        out.resultFile = result

        if (maskedAssignments) {
            List<Map> cells = run.parsed.cells as List<Map>
            out.solverRun = !(run.parsed.solver_calls as List).isEmpty()
            out.facts = baseFacts(c, modelInfo, guided ? report : [context_source: 'explicit_masked_assignments'], c.experimentMode) + [
                analysis_type: 'masked_assignment', masked_model_sha256: run.parsed.model.sha256,
                public_problem_sha256: run.parsed.public_problem_sha256,
                solver_calls_used: (run.parsed.solver_calls as List).size(), solver_calls_max: run.parsed.max_solver_calls,
                outcome: cells && cells.every { it.status == 'completed' } && !run.parsed.unvisited_conditions ?
                    'completed' : 'completed_with_skips_or_inconclusive',
                hypotheses: run.parsed.public_problem.hypotheses.size(),
                conditions_visited: run.parsed.visited_conditions, conditions_unvisited: run.parsed.unvisited_conditions,
                comparisons_completed: cells.count { it.status == 'completed' },
                comparisons_inconclusive: cells.count { it.status != 'completed' },
                candidates_with_conditional_features: 0, promotion: 'none', seconds: run.seconds,
                result_file: result.path, result_sha256: run.outputSha256,
                interpretation: 'Target gene assignments were masked before scenario construction. These are conditional phenotype signatures, not verified function recovery or learned evidence features.']
            return out
        }

        Map<String, String> basis = modelInfo.table.collectEntries { rid, row -> [rid, row.basis] }
        ExperimentContextMapper.Mapped mapped = ExperimentContextMapper.map(run.parsed, bindings, basis)
        out.fba = mapped.features
        List<Map> experiments = ((run.parsed.experiments ?: []) as List<Map>).findAll { it.experiment_id != 'baseline' }
        Map baseRec = ((run.parsed.experiments ?: []) as List<Map>).find { it.experiment_id == 'baseline' }
        Map baseMeasurements = (baseRec?.measurements ?: [:]) as Map
        Map<String, Object> baselineFacts = [
            status: baseRec == null ? 'not_run' : baseRec.status, outcome: baseRec?.outcome,
            reactions_checked: ((baseMeasurements.focus_checks ?: []) as List<Map>)*.reaction,
            reactions_not_checked: baseMeasurements.focus_reactions_not_checked ?: [],
            note: baseMeasurements.focus_checks ? 'matched control: FVA and a reaction knockout in the unchanged medium; model feasibility only' :
                'no baseline check was run (no growth, no budget or no mapped reaction): unavailable, not negative']
        out.facts = baseFacts(c, modelInfo, report, c.experimentMode) + [
            baseline_focus: baselineFacts,
            outcome: run.parsed.outcome,
            solver_calls_used: (run.parsed.solver as Map).used,
            solver_calls_max: (run.parsed.solver as Map).max_calls,
            seconds: run.seconds,
            experiments_completed: experiments.count { it.status == 'completed' },
            experiments_skipped_or_inconclusive: experiments.count { it.status != 'completed' },
            menu_rejections: (run.parsed.menu_rejections ?: []).size(),
            candidates_with_conditional_features: mapped.features.size(),
            candidates_flagged_objective_requires_candidate: mapped.features.count { k, f -> f.objectiveRequiresCandidate },
            conditions_by_candidate: mapped.conditions,
            own_knockouts_not_used_as_features: mapped.ownKnockouts,
            result_sha256: run.outputSha256,
            result_file: result.path,
            promotion: 'none',
            interpretation: 'Conditional model predictions under the stated medium, bounds and maintenance constraint; not observed phenotypes and not independent proof of gene function. A product or growth objective that needs the candidate is flagged and kept out of the evidence features.']
        out
    }

    private Map<String, Object> baseFacts(IntegrationWorkflowConfig c, SbmlModelInfo model, Map report, String mode) {
        [mode: mode, strategy: c.experimentStrategy, seed: c.experimentSeed, model_file: model.model.path, model_sha256: model.sha256,
         max_experiments: c.experimentMaxExperiments, max_solver_calls: c.experimentMaxSolverCalls,
         context_source: report.context_source, context: report,
         model_origin: model.hasTable() ? 'refiner export with reactions.tsv (gene rules inferred from this run\'s EC evidence)' : 'supplied model without a refiner table'] as Map<String, Object>
    }

    private static File resolveModel(Map template, File templateFile, File refinerModel, boolean guided) {
        File named = template.model ? resolve(templateFile, template.model as String) : null
        if (named != null && refinerModel != null && guided) {
            throw new IllegalArgumentException("the template names a model and this run also exported one (metabolic_refinement): " +
                "remove 'model' from the template to use the run's model, or disable the refiner")
        }
        File model = named ?: (guided ? refinerModel : null)
        if (model == null) {
            throw new IllegalArgumentException(guided ?
                "evidence_guided mode needs a model: name one in the template or run metabolic_refinement in this run" :
                "the experiment manifest names no 'model'")
        }
        if (!model.isFile()) throw new IllegalArgumentException("experiment model does not exist: ${model}")
        model
    }

    private static File resolve(File base, String path) {
        File f = new File(path)
        f.isAbsolute() ? f : new File(base.absoluteFile.parentFile, path)
    }

    /** Hypotheses of a user manifest that name the function they test (function_key) can be mapped to features. */
    private static List<ExperimentManifestBuilder.Binding> bindingsFromManifest(Map manifest) {
        ((manifest.context ?: [:]).hypotheses ?: []).findAll { it.function_key }.collect { Map h ->
            new ExperimentManifestBuilder.Binding(hypothesisId: h.id as String, functionKey: h.function_key as String,
                reactions: ((h.reactions ?: []) as List<String>))
        }
    }

    /** Content checks: the result belongs to this model and budget and carries no confidence or promotion. */
    static void verifyMasked(Map parsed, SbmlModelInfo model, IntegrationWorkflowConfig c, Map requested) {
        if (parsed.input_model_sha256 != model.sha256 || !(parsed.model?.sha256 ==~ /[0-9a-f]{64}/) ||
            parsed.public_problem?.model != "masked:${parsed.model.sha256}".toString()) {
            throw new SidecarRun.SidecarFailure('masked result has missing or mismatched model provenance')
        }
        if (parsed.public_problem.hypotheses != requested.assignment_recovery.hypotheses ||
            parsed.public_problem.conditions != requested.assignment_recovery.conditions ||
            parsed.public_problem.baseline != requested.baseline || parsed.public_problem.objectives != requested.objectives ||
            (parsed.public_problem.condition_constraints ?: []) != (requested.condition_constraints ?: [])) {
            throw new SidecarRun.SidecarFailure('masked result changed the candidate assignments, conditions or model constraints')
        }
        if (parsed.strategy != c.experimentStrategy || parsed.max_solver_calls != c.experimentMaxSolverCalls ||
            !(parsed.solver_calls instanceof List) || parsed.solver_calls.size() > c.experimentMaxSolverCalls ||
            !(parsed.visited_conditions instanceof List) || parsed.visited_conditions.size() > c.experimentMaxExperiments ||
            parsed.promotion != 'none' || parsed.calibrated_probability != null) {
            throw new SidecarRun.SidecarFailure('masked result violates the requested budget, strategy or interpretation')
        }
        List<Map> hypotheses = parsed.public_problem.hypotheses as List<Map>
        List<Map> cells = parsed.cells as List<Map>
        Set pairs = cells.collect { [it.condition_id, it.hypothesis_id] } as Set
        Set expected = [] as Set
        parsed.visited_conditions.each { condition -> hypotheses.each { h -> expected.add([condition, h.id]) } }
        if (hypotheses.count { !it.assignments } != 1 || cells.size() != pairs.size() || pairs != expected) {
            throw new SidecarRun.SidecarFailure('masked result lacks a complete comparison including the null assignment')
        }
        cells.each { cell ->
            if (!(cell.status in ['completed', 'inconclusive']) ||
                cell.phenotypes?.keySet() != ['wild_type', 'target_knockout'] as Set ||
                !((cell.solver_calls ?: []) as List).every { it instanceof Number && it.intValue() >= 0 && it.intValue() < parsed.solver_calls.size() }) {
                throw new SidecarRun.SidecarFailure('invalid masked comparison or solver-call references')
            }
        }
    }

    static void verify(Map parsed, SbmlModelInfo model, IntegrationWorkflowConfig c) {
        String sha = (parsed.hashes as Map)?.model_sha256
        if (sha != model.sha256) {
            throw new SidecarRun.SidecarFailure("the result was computed on a different model (sha256 ${sha} != supplied ${model.sha256}); refusing it")
        }
        Map solver = parsed.solver as Map
        if (solver == null || (solver.max_calls as Number)?.intValue() != c.experimentMaxSolverCalls ||
            (solver.used as Number).intValue() > (solver.max_calls as Number).intValue()) {
            throw new SidecarRun.SidecarFailure("the result does not honour the requested solver-call budget of ${c.experimentMaxSolverCalls}: ${solver}")
        }
        if (parsed.strategy != c.experimentStrategy) {
            throw new SidecarRun.SidecarFailure("the result used strategy '${parsed.strategy}', not the requested '${c.experimentStrategy}'")
        }
        if (parsed.promotion != 'none') throw new SidecarRun.SidecarFailure("the result carries promotion '${parsed.promotion}'; experiments never promote annotations")
        if (!(parsed.outcome in KNOWN_OUTCOMES)) throw new SidecarRun.SidecarFailure("unknown experiment outcome '${parsed.outcome}'")
        ((parsed.conditional_evidence ?: []) as List<Map>).each { Map r ->
            if (r.calibrated_probability != null || r.independent_proof_of_function != false || r.promotion != 'none') {
                throw new SidecarRun.SidecarFailure("conditional record ${r.experiment_id} claims a probability, a promotion or independent proof")
            }
        }
    }
}
