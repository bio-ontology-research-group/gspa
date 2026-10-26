package gspa.integration

import com.fasterxml.jackson.databind.ObjectMapper

/**
 * Per-run evidence ledger shared by integration entry points. Changes describe
 * predictions, never biological correctness; accuracy needs independent labels.
 */
class IntegrationRunAudit {
    static final List<String> STAGES = ['evidence_integration', 'dark_matter',
        'reaction_context', 'cross_genome', 'symbolic_enforcement', 'metabolic_refinement',
        'agent_completion'].asImmutable()

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
                                    Closure<IntegratedAnnotationSet> operation) {
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
            record.reason = "${ex.class.simpleName}: ${ex.message}"
            throw ex
        } finally {
            record.elapsed_seconds = (System.nanoTime() - start) / 1e9d
        }
    }

    Map report() {
        [schema_version: 1, genome_id: genomeId,
         interpretation: 'Prediction changes are not measurements of accuracy.',
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
        Set<String> keys = new LinkedHashSet<>(before.keySet())
        keys.addAll(after.keySet())
        List<Map> out = []
        keys.each { key ->
            Double old = before[key], now = after[key]
            if (old != null && now != null && Math.abs(old - now) <= 1e-12d) return
            String[] parts = IntegrationState.splitFunctionKey(key)
            out << [function_key: key, protein_id: parts[0], function_type: parts[1], function_id: parts[2],
                    change: old == null ? 'added' : (now == null ? 'removed' : (now > old ? 'increased' : 'decreased')),
                    before_probability: old, after_probability: now]
        }
        out
    }
}
