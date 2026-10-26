package gspa.integration.workflow

import com.fasterxml.jackson.databind.ObjectMapper
import gspa.integration.CalibrationTable
import gspa.integration.ClaimExtractor
import gspa.integration.learned.LearnedTestSupport
import spock.lang.Specification
import spock.lang.TempDir

/**
 * Writes the real per-genome artifacts (`integrated.tsv`, `.stages.json`, `.learned.json`, provenance) of a
 * learned-stage run through the shared workflow, for the results-browser and portal tests
 * (`scripts/tests/test_visualize_learned_e2e.py`, which sets GSPA_LEARNED_E2E_DIR and reads the result).
 *
 * The model is a frozen, hand-set, TEST ONLY mechanism artifact (`intended_use = mechanism_test`). It was never
 * trained and nothing here is evidence about prediction quality; the point is that the JVM writers and the
 * Python readers agree on the files. The production loader refuses such a model, so the run goes through the
 * test-source seam `LearnedWorkflowSpec.mechanismWorkflow`, not through `gspa integrate`.
 */
class LearnedArtifactsWriterSpec extends Specification {
    @TempDir File tmp
    final ObjectMapper mapper = new ObjectMapper()

    def 'the shared workflow writes a ledger and a learned report that name the same frozen test model'() {
        given:
        String env = System.getenv('GSPA_LEARNED_E2E_DIR')
        File dir = env ? new File(env) : tmp
        dir.mkdirs()
        def fx = new ContextFixtures(dir)
        def files = fx.genome('g1')
        def model = LearnedTestSupport.logisticModel(['seq.group.homology.max_raw': 2.0d], -1.0d)
        File placeholder = new File(dir, 'model.json')
        placeholder.text = '{"placeholder": true}'
        def config = fx.worker(fx.config(learnedModelFile: placeholder), 'g1', files, new File(dir, 'integrated.tsv'))
        config.provenanceFile = new File(dir, 'provenance.json')
        def claims = new ClaimExtractor(new CalibrationTable()).readClaimsJsonl(files.claims)

        when:
        def wf = LearnedWorkflowSpec.mechanismWorkflow(config, model)
        def result = wf.run(claims)
        wf.writeArtifacts(result)
        def ledger = mapper.readValue(new File(dir, 'integrated.tsv.stages.json'), Map)
        def report = mapper.readValue(new File(dir, 'integrated.tsv.learned.json'), Map)

        then:
        ledger.scoring_mode == 'learned'
        ledger.stages.find { it.name == 'learned_integration' }.status == 'executed'
        ledger.stages.find { it.name == 'learned_integration' }.model_sha256 == model.sha256
        report.schema == 'gspa.learned-stage-report/1'
        report.model.sha256 == model.sha256
        report.model.intended_use == 'mechanism_test'
        report.candidates.size() > 0
        new File(dir, 'provenance.json').isFile()
    }
}
