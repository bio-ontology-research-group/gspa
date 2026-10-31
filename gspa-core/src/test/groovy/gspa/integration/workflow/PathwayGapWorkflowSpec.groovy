package gspa.integration.workflow

import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import gspa.integration.learned.ContextFeatureExtractor
import spock.lang.Specification
import spock.lang.TempDir

import static gspa.integration.learned.LearnedTestSupport.claim
import static gspa.integration.workflow.PathwayGapBuilderSpec.observation
import static gspa.integration.workflow.PathwayGapBuilderSpec.row

class PathwayGapWorkflowSpec extends Specification {
    @TempDir File dir
    def claims = [claim('p', 'EC:1.1.1.1', 0.2d)]

    private IntegrationWorkflowConfig config() {
        new IntegrationWorkflowConfig(lite:true, enabledPriors:'none', genomeId:'g', outFile:new File(dir,'out.tsv'))
    }
    private double available(IntegrationWorkflow.Result r) {
        new ContextFeatureExtractor().extractAll(claims, r.state)[0].get('ctx.gap.available')
    }

    def 'in-memory and file observations produce the same gaps, features and replayable artifacts'() {
        given:
        def obs = observation()
        def memory = new IntegrationWorkflow(config())
        memory.pathwaySearchObservations = [obs]
        def file = new File(dir, 'observations.json')
        file.text = JsonOutput.toJson(obs)
        def cfg = config(); cfg.pathwayObservationsFile = file
        def fromFile = new IntegrationWorkflow(cfg)

        when:
        def a = memory.run(claims)
        def b = fromFile.run(claims)
        fromFile.writeArtifacts(b)
        def replayCfg = config(); replayCfg.gapsFile = new File(cfg.outFile.path + '.gaps.jsonl')
        // Reusing the same output prefix must not delete the explicitly supplied gap input.
        def replayWorkflow = new IntegrationWorkflow(replayCfg)
        def replay = replayWorkflow.run(claims)
        replayWorkflow.writeArtifacts(replay)

        then:
        a.state.metabolicGaps == b.state.metabolicGaps
        replay.state.metabolicGaps == b.state.metabolicGaps
        available(a) == available(b) && available(b) == available(replay) && available(b) == 1d
        a.state.metabolicGapReport == b.state.metabolicGapReport
        new JsonSlurper().parse(new File(cfg.outFile.path+'.gap-context.json')).available
        b.audit.report().stages.find { it.name=='reference_preparation' }.metabolic_gap_records == 1
    }

    def 'observed zero and absent gap inputs have distinct feature availability'() {
        given:
        def wf = new IntegrationWorkflow(config())
        wf.pathwaySearchObservations = [observation([row(status:'good_blast')])]

        when:
        def observed = wf.run(claims)
        def absent = new IntegrationWorkflow(config()).run(claims)

        then:
        observed.state.metabolicGaps.isEmpty()
        absent.state.metabolicGaps.isEmpty()
        available(observed) == 1d
        available(absent) == 0d
    }

    def 'explicit gap set takes precedence over automatic gaps including explicit empty set'() {
        given:
        def cfg = config(); cfg.gapsFile = new File(dir,'empty.jsonl'); cfg.gapsFile.text = ''
        def wf = new IntegrationWorkflow(cfg)
        wf.pathwaySearchObservations = [observation()]

        when:
        def result = wf.run(claims)

        then:
        result.state.metabolicGaps.isEmpty()
        available(result) == 1d
        result.state.metabolicGapReport.source == 'explicit_gaps'
        result.state.metabolicGapReport.observed_searches_superseded == 1
    }

    def 'partial observations remain unavailable and obsolete gap artifacts are removed'() {
        given:
        def cfg = config()
        def wf = new IntegrationWorkflow(cfg)
        wf.pathwaySearchObservations = [observation()]
        wf.writeArtifacts(wf.run(claims))
        assert new File(cfg.outFile.path+'.gaps.jsonl').exists()
        wf.pathwaySearchObservations = [observation([], [outcome:'PARTIAL'])]

        when:
        def partial = wf.run(claims)
        wf.writeArtifacts(partial)

        then:
        available(partial) == 0d
        partial.state.metabolicGaps.isEmpty()
        !new File(cfg.outFile.path+'.gaps.jsonl').exists()
        !new JsonSlurper().parse(new File(cfg.outFile.path+'.gap-context.json')).available
    }

    def 'cross-genome observations fail preparation and leave a failed ledger'() {
        given:
        def cfg = config()
        def wf = new IntegrationWorkflow(cfg)
        wf.pathwaySearchObservations = [observation([], [genome_id:'wrong'])]

        when:
        wf.run(claims)

        then:
        thrown(IllegalArgumentException)
        new JsonSlurper().parse(new File(cfg.outFile.path+'.stages.json')).stages.find {
            it.name=='reference_preparation'
        }.status == 'failed'
    }
}
