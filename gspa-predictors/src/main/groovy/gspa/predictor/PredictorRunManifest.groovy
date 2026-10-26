package gspa.predictor

import groovy.json.JsonOutput

/**
 * Machine-readable record of every predictor's completion status in a run.
 * Written next to the annotation outputs as {@code <genome>_predictor_status.json}.
 */
class PredictorRunManifest {

    static final String SCHEMA = 'gspa.predictor-run-manifest/1'

    static Map<String, Object> toMap(String genomeId, Collection<PredictorRunStatus> statuses) {
        Map<String, Integer> counts = [:]
        PredictorOutcome.values().each { counts[it.name()] = 0 }
        statuses.each { counts[it.outcome.name()] += 1 }
        [
            schema    : SCHEMA,
            genome_id : genomeId,
            summary   : counts,
            predictors: statuses.collect { it.toMap() },
        ] as Map<String, Object>
    }

    static void write(String genomeId, Collection<PredictorRunStatus> statuses, File output) {
        output.parentFile?.mkdirs()
        output.text = JsonOutput.prettyPrint(JsonOutput.toJson(toMap(genomeId, statuses))) + '\n'
    }
}
