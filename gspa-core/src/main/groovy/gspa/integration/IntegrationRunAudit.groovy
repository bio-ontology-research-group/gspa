package gspa.integration

import com.fasterxml.jackson.databind.ObjectMapper

/**
 * Per-run evidence ledger shared by integration entry points. Changes describe
 * predictions, never biological correctness; accuracy needs independent labels.
 */
class IntegrationRunAudit {
    static final List<String> STAGES = ['reference_preparation', 'cohort_context', 'evidence_integration', 'dark_matter',
        'reaction_context', 'cross_genome', 'symbolic_enforcement', 'metabolic_refinement',
        'agent_completion', 'metabolic_experiments', 'community_metabolism', 'learned_integration', 'feature_export'].asImmutable()

    /** Meaning of the two baselines a ledger can mention; never to be confused. */
    static final Map<String, String> BASELINE_DEFINITIONS = [
        stage0_raw_max_score : 'Stage 0: strongest raw predictor score per function; no priors, no combination, no posterior.',
        legacy_prior_baseline: 'Posterior after the evidence_integration stage (heuristic evidence combiner and priors). The ' +
                               'baseline_candidates and final_changes of this ledger are measured from it.',
    ].asImmutable()

    /** legacy_heuristic until the learned stage has replaced the scores; then learned. */
    String scoringMode = 'legacy_heuristic'

    String genomeId
    final Map<String, Map> stages = new LinkedHashMap<>()
    private Map<String, Double> baseline = [:]
    private Map<String, Double> current = [:]

    IntegrationRunAudit(String genomeId) {
        this.genomeId = genomeId
        STAGES.each { stages[it] = [name: it, status: 'not_run', changes: []] }
    }

    IntegratedAnnotationSet runStage(String name, boolean enabled, String unavailable,
                                    IntegratedAnnotationSet input, IntegrationState state,
                                    Closure<IntegratedAnnotationSet> operation, String noWorkReason = null) {
        if (!stages.containsKey(name)) throw new IllegalArgumentException("Unknown integration stage: ${name}")
        Map record = stages[name]
        if (record.status != 'not_run') throw new IllegalStateException("Stage already recorded: ${name}")
        if (!enabled) {
            record.status = 'skipped'
            record.reason = 'disabled'
            return input
        }
        long start = System.nanoTime()
        Map<String, Double> before = scores(input)
        try {
            if (unavailable) throw new IllegalStateException("${name}: ${unavailable}")
            // Eligibility is distinct from execution. Validate required resources first: an empty
            // workload must never hide a broken configuration. Preserve the preceding prediction
            // state and record why this requested stage had nothing to evaluate.
            if (noWorkReason) {
                record.status = 'skipped'
                record.reason = noWorkReason.toString()
                record.requested = true
                record.input_candidates = before.size()
                record.output_candidates = before.size()
                record.suggestions = input?.suggestions?.size() ?: 0
                return input
            }
            IntegratedAnnotationSet output = operation.call()
            if (output == null) throw new IllegalStateException("${name} returned no integration result")
            output.synchronize(state)
            current = scores(output)
            record.status = 'executed'
            record.input_candidates = before.size()
            record.output_candidates = current.size()
            record.suggestions = output.suggestions.size()
            record.changes = changes(before, current)
            if (name == 'evidence_integration') baseline = new LinkedHashMap<>(current)
            return output
        } catch (Exception ex) {
            record.status = 'failed'
            record.reason = "${ex.class.simpleName}: ${ex.message}".toString()
            throw ex
        } finally {
            record.elapsed_seconds = (System.nanoTime() - start) / 1e9d
        }
    }

    /**
     * Run a stage that produces context or artifacts but no annotation set (cohort inference,
     * metabolic experiments, feature export): skipped when disabled, failed with the reason when it
     * is unavailable or throws, executed otherwise. The closure may return a map of facts that is
     * attached to the stage's record.
     */
    void runSupport(String name, boolean enabled, String unavailable, Closure operation) {
        if (!stages.containsKey(name)) throw new IllegalArgumentException("Unknown integration stage: ${name}")
        Map record = stages[name]
        if (record.status != 'not_run') throw new IllegalStateException("Stage already recorded: ${name}")
        if (!enabled) {
            record.status = 'skipped'
            record.reason = 'disabled'
            return
        }
        long start = System.nanoTime()
        try {
            if (unavailable) throw new IllegalStateException("${name}: ${unavailable}")
            Object facts = operation.call()
            record.status = 'executed'
            if (facts instanceof Map) record.putAll(facts as Map)
        } catch (Exception ex) {
            record.status = 'failed'
            record.reason = "${ex.class.simpleName}: ${ex.message}".toString()
            throw ex
        } finally {
            record.elapsed_seconds = (System.nanoTime() - start) / 1e9d
        }
    }

    /**
     * Run a stage that produces no annotation set (configuration validation, reference loading):
     * recorded executed, or failed with the reason when it throws.
     */
    void runPreparation(String name, Closure operation) {
        Map record = stages[name]
        long start = System.nanoTime()
        try {
            operation.call()
            record.status = 'executed'
        } catch (Exception ex) {
            record.status = 'failed'
            record.reason = "${ex.class.simpleName}: ${ex.message}".toString()    // a GString would serialize as a bean
            throw ex
        } finally {
            record.elapsed_seconds = (System.nanoTime() - start) / 1e9d
        }
    }

    /** After an early failure, say why every requested stage that never started did not run. */
    void abandon(Collection<String> requested, String reason) {
        requested.each { String name ->
            Map record = stages[name]
            if (record != null && record.status == 'not_run') record.reason = reason
        }
    }

    /** Attach stage-specific facts (assumptions, inputs, counts) to a stage's record. */
    void note(String name, Map<String, Object> details) {
        if (!stages.containsKey(name)) throw new IllegalArgumentException("Unknown integration stage: ${name}")
        stages[name].putAll(details)
    }

    Map report() {
        [schema_version: 2, genome_id: genomeId,
         interpretation: 'Prediction changes are not measurements of accuracy.',
         scoring_mode: scoringMode,
         scoring_mode_note: scoringMode == 'learned' ?
             'Annotation scores were replaced by a frozen learned model (see learned_integration).' :
             scoringMode == 'mixed' ? 'Supported candidates use the frozen model; unsupported candidates retain legacy scores (see learned_integration).' :
             'Legacy heuristic scoring (evidence combiner and hand-set priors). No trained model produced these scores.',
         baseline_definitions: BASELINE_DEFINITIONS,
         changes_measured_from: 'legacy_prior_baseline',
         stages: new ArrayList<>(stages.values()),
         baseline_candidates: baseline.size(), final_candidates: current.size(),
         final_changes: changes(baseline, current)]
    }

    void write(File target) {
        target.parentFile?.mkdirs()
        new ObjectMapper().writerWithDefaultPrettyPrinter().writeValue(target, report())
    }

    private static Map<String, Double> scores(IntegratedAnnotationSet integrated) {
        Map<String, Double> out = new LinkedHashMap<>()
        integrated?.provenance?.each { key, p -> out[key] = p.finalProbability }
        out
    }

    private static List<Map> changes(Map<String, Double> before, Map<String, Double> after) {
        AuditChanges.between(before, after)
    }
}
