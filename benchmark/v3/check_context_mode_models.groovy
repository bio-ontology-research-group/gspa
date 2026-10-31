import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import gspa.integration.cohort.ContextMode
import gspa.integration.learned.CandidateFeatures
import gspa.integration.learned.CandidateRole
import gspa.integration.learned.LearnedContextModel
import gspa.model.AnnotationType

// Actual production model API, including mode-aware coverage and calibrated scores.
def model = LearnedContextModel.load(new File(args[0]))
def names = model.contract.names()
new File(args[2]).withWriter('UTF-8') { writer ->
    new File(args[1]).eachLine('UTF-8') { line ->
        if (!line.trim()) return
        def row = new JsonSlurper().parseText(line)
        double[] vector = names.collect { row.features[it] as double } as double[]
        def f = new CandidateFeatures(row.protein_id, AnnotationType.valueOf(row.function_type),
            row.function_id, CandidateRole.valueOf(row.role), model.contract, vector, [], row.sequence_sources,
            row.quality_evidence ?: [], ContextMode.valueOf(row.cohort_context_mode ?: 'UNIVERSAL'), row.sequence_contracts)
        String reason = model.unsupportedReason(f)
        writer.write(JsonOutput.toJson([example_id:row.example_id, mode:f.cohortContextMode.name(),
            unsupported_reason:reason, probability:reason ? null : model.score(f).probability]) + '\n')
    }
}
