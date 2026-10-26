package gspa.agent

import gspa.integration.ClaimProvenance
import gspa.integration.IntegratedAnnotationSet
import gspa.integration.IntegrationState
import spock.lang.Specification
import spock.lang.TempDir

class AgentExecutionSpec extends Specification {
    @TempDir File dir

    def 'failed agent request cannot promote a protein using a stale resolution from a previous run'() {
        given:
        def state = new IntegrationState()
        state.complexClassification = ['GO:0000124': 'n']
        def integrated = new IntegratedAnnotationSet()
        ['subunit': 0.9d, 'candidate': 0.3d].each { gene, probability ->
            String key = "${gene}|GO|GO:0000124"
            integrated.provenance[key] = new ClaimProvenance(functionKey: key, proteinId: gene,
                finalProbability: probability, finalLogOdds: Math.log(probability / (1d - probability)))
        }
        integrated.synchronize(state)
        assert !new AnomalyExtractor().extract(state, integrated).isEmpty()
        File stale = new File(dir, 'resolutions.jsonl')
        stale.text = '{"action":"PROMOTE","protein_id":"candidate","function_type":"GO","function_id":"GO:0000124","confidence":0.99}\n'
        File sidecar = new File(dir, 'failed.py')
        sidecar.text = 'raise SystemExit(7)\n'
        def runner = new AgentCompletionRunner(sidecarScript: sidecar)

        when:
        runner.run(state, integrated, new File(dir, 'anomalies.jsonl'), stale)

        then:
        def error = thrown(IllegalStateException)
        error.message.contains('exit=7')
        integrated.provenance['candidate|GO|GO:0000124'].finalProbability == 0.3d
        integrated.annotationsByFunction['candidate|GO|GO:0000124'].score == 0.3d
    }
}
